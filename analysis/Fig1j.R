#! usr/bin/env Rscript

# GSEA for DEGs

if(!("pacman" %in% rownames(installed.packages()))){
  install.packages("pacman")
}

pacman::p_load(data.table,
               ggplot2,
               stringr,
               dplyr,
               fgsea,
               )

# Enter paths here
output_dir <- ""

# Creating directories needed for outputs
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

plot_dir <- paste0(output_dir, "/path_to_output_directory")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# GSEA
 
curr_gene_list <- read.csv("de_results.csv") # import de results from fig 1i

rownames(curr_gene_list) <- curr_gene_list$gene

gene_set <- gmtPathways("/cluster/projects/gaitigroup/Users/Benson/Parsebio/gene_lists/c5.go.v2023.2.Hs.symbols.gmt")

ranked_genes_response <- curr_gene_list %>%
    mutate(rank = log2FoldChange) %>%
    filter(rank != Inf) %>% filter(rank != -Inf) %>% arrange(-rank)

ranked_genes_response$ranking <- rank(-ranked_genes_response$rank, ties.method = "first")
ggplot(ranked_genes_response, aes(x = ranking, y = rank)) +
    geom_bar(stat = "identity") +
    ylab("Ranked List Metric (log2FoldChange)") +
    xlab("Rank in Ordered Dataset")

write.csv(ranked_genes_response, file.path(plot_dir, "ranked_genes_response.csv"))
response_level_stats <- ranked_genes_response$rank
names(response_level_stats) <- rownames(ranked_genes_response)

# Run gsea and filter the pathways with a p-val less than 0.05
fgseaRes_ctrl_multilevel <- fgseaMultilevel(pathways = gene_set,
                                            stats = response_level_stats,
                                            minSize=10,
                                            maxSize=2000,
                                            nPermSimple = 10000) %>% 
filter(padj < 0.5) %>% arrange(-abs(NES))


# Format for Cytoscape EnrichmentMap

current_ranks <- ranked_genes_response$ranking
names(current_ranks) <- ranked_genes_response$gene

# function borrowed from Ruth Isserlin from Bader lab 
write_patient_fgsea_results<- function(current_fgsea_results, current_results_dir, 
                                       file_header){
      
    directory_fullpath <- file.path(current_results_dir, 
                                                    file_header)
    
    if(!dir.exists(directory_fullpath)){
      dir.create(directory_fullpath)
    }
  
    #calculate the rank at max
    #fgsea returns the leading edge.  Just need to extract the highest rank from 
    # set to get the rank at max
    calculated_rank_at_max <- apply(current_fgsea_results,1,
                                    FUN=function(x){ 
                                      max(which(names(current_ranks)
                                                %in% unlist(x[8])))})
    
    
    enr <- cbind(current_fgsea_results$pathway,
                                     current_fgsea_results$pathway,
                                     "Details",
                                     current_fgsea_results$size,
                                     current_fgsea_results$ES,
                                     current_fgsea_results$NES,
                                     current_fgsea_results$pval,
                                     current_fgsea_results$padj,
                                     0,
                                     calculated_rank_at_max,
                                     apply(current_fgsea_results,1,
                                           FUN=function(x){paste(unlist(x[8]),collapse=",")})) 
    
    colnames(enr) <- c("name","description",
                                           "GS details","SIZE","ES",
                                           "NES","pval","padj","FWER",
                                           "Rank at Max","leading edge genes")
    
    enr_filename_positive <- paste0(file_header, "fgsea_enr_results_pos.txt",sep="")
    enr_filename_docker_positive <- file.path(directory_fullpath,
                                                  enr_filename_positive)
   
    write.table(enr[which(enr[,"NES"]>=0),] ,
                enr_filename_docker_positive,
                col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE,fileEncoding="latin1")
    
    enr_filename_negative <- paste0(file_header, "fgsea_enr_results_neg.txt",sep="")
    enr_filename_docker_negative <- file.path(directory_fullpath,
                                                  enr_filename_negative)
   
    write.table(enr[which(enr[,"NES"]<0),] ,
                enr_filename_docker_negative,
                col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE,fileEncoding="latin1")
    
      
}

write_patient_fgsea_results(fgseaRes_ctrl_multilevel, plot_dir, "opcnpc_pt_vs_tum_C5_GO")

