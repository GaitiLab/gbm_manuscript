import argparse
import logging
import os
import time
from argparse import ArgumentParser as AP
from os.path import abspath
from pathlib import Path

import anndata
import geopandas as gpd
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def write_10X_h5(adata, file):
    """Writes adata to a 10X-formatted h5 file.
    Ref: https://github.com/scverse/anndata/issues/595#issuecomment-1824376236

    Note that this function is not fully tested and may not work for all cases.
    It will not write the following keys to the h5 file compared to 10X:
    '_all_tag_keys', 'pattern', 'read', 'sequence'

    Args:
        adata (AnnData object): AnnData object to be written.
        file (str): File name to be written to. If no extension is given, '.h5' is appended.

    Raises:
        FileExistsError: If file already exists.

    Returns:
        None
    """

    if ".h5" not in file:
        file = f"{file}.h5"
    if Path(file).exists():
        raise FileExistsError(f"There already is a file `{file}`.")

    def int_max(x):
        return int(max(np.floor(len(str(int(max(x)))) / 4), 1) * 4)

    def str_max(x):
        return max([len(i) for i in x])

    w = h5py.File(file, "w")
    grp = w.create_group("matrix")
    grp.create_dataset(
        "barcodes",
        data=np.array(adata.obs_names, dtype=f"|S{str_max(adata.obs_names)}"),
    )
    grp.create_dataset(
        "data", data=np.array(adata.X.data, dtype=f"<i{int_max(adata.X.data)}")
    )
    ftrs = grp.create_group("features")
    # this group will lack the following keys:
    # '_all_tag_keys', 'feature_type', 'genome', 'id', 'name', 'pattern', 'read', 'sequence'
    ftrs.create_dataset(
        "feature_type",
        data=np.array(
            adata.var.feature_types, dtype=f"|S{str_max(adata.var.feature_types)}"
        ),
    )
    ftrs.create_dataset(
        "genome",
        data=np.array(adata.var.genome, dtype=f"|S{str_max(adata.var.genome)}"),
    )
    ftrs.create_dataset(
        "id",
        data=np.array(adata.var.gene_ids, dtype=f"|S{str_max(adata.var.gene_ids)}"),
    )
    ftrs.create_dataset(
        "name", data=np.array(adata.var.index, dtype=f"|S{str_max(adata.var.index)}")
    )
    grp.create_dataset(
        "indices", data=np.array(adata.X.indices, dtype=f"<i{int_max(adata.X.indices)}")
    )
    grp.create_dataset(
        "indptr", data=np.array(adata.X.indptr, dtype=f"<i{int_max(adata.X.indptr)}")
    )
    grp.create_dataset(
        "shape",
        data=np.array(list(adata.X.shape)[::-1], dtype=f"<i{int_max(adata.X.shape)}"),
    )


def get_args():
    # Script description
    description = """Mapping barcodes + generate Anndata object"""

    # Add parser
    parser = AP(
        description=description, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--sample_dir",
        type=str,
        help="Path to spatial outs directory for `square_002um`",
    )
    parser.add_argument("-id", "--sample_id", type=str, help="Sample id")
    parser.add_argument(
        "-o", "--output_dir", type=str, default="", help="Output directory"
    )
    parser.add_argument(
        "--gdf_barcodes_path",
        type=str,
        help="Path to GeoDataFrame with barcodes polygons, i.e. <sample_id>__tissue_positions__poly.parquet",
    )
    parser.add_argument(
        "--gdf_segmentation_path",
        type=str,
        help="Path to GeoDataFrame with polygons of segmented cells, i.e. <sample_id>__stardist_masks__expanded__poly.parquet",
    )
    arg = parser.parse_args()
    arg.output_dir = Path(abspath(arg.output_dir))

    if (arg.output_dir != "") & (not os.path.isdir(arg.output_dir)):
        arg.output_dir.mkdir(parents=True, exist_ok=True)

    arg.gdf_barcodes_path = Path(arg.gdf_barcodes_path)
    arg.gdf_segmentation_path = Path(arg.gdf_segmentation_path)
    arg.sample_dir = Path(arg.sample_dir)
    return arg


def map_barcodes(
    sample_dir: str,
    gdf_barcodes_path: str,
    gdf_segmentation_path: str,
    output_dir: str,
    sample_id: str,
):
    # Load gene-expression
    logging.info("Load feature matrix...")
    adata = sc.read_10x_h5(Path(sample_dir, "filtered_feature_bc_matrix.h5"))

    logging.info("Load tissue positions...")
    df_tissue_positions = pd.read_parquet(
        Path(sample_dir, "spatial", "tissue_positions.parquet")
    )
    # Set the index of the dataframe to the barcodes
    df_tissue_positions = df_tissue_positions.set_index("barcode")

    # Create an index in the dataframe to check joins
    df_tissue_positions["index"] = df_tissue_positions.index
    # Adding the tissue positions to the meta data
    adata.obs = pd.merge(
        adata.obs, df_tissue_positions, left_index=True, right_index=True
    )

    # Loading polygons
    logging.info("Load polygons (squares) for bins/barcodes...")
    gdf_barcodes = gpd.read_parquet(gdf_barcodes_path)
    logging.info("Load polygons (segmentations)...")
    gdf_segmentation = gpd.read_parquet(gdf_segmentation_path)

    logging.info("Determine whether bin is within a cell/nuclei masks...")
    # NOTE this is different from previous approach, comparing the bin polygon and the segmentation polygon
    result_spatial_join = gpd.sjoin(
        gdf_barcodes, gdf_segmentation, how="left", predicate="covered_by"
    )

    # Identify nuclei associated barcodes and find barcodes that are in more than one nucleus
    result_spatial_join["is_within_polygon"] = ~result_spatial_join[
        "index_right"
    ].isna()
    barcodes_in_overlaping_polygons = pd.unique(
        result_spatial_join[result_spatial_join.duplicated(subset=["index"])]["index"]
    )
    result_spatial_join["is_not_in_an_polygon_overlap"] = ~result_spatial_join[
        "index"
    ].isin(barcodes_in_overlaping_polygons)

    # Remove barcodes based on filters
    barcodes_in_one_polygon = result_spatial_join[
        result_spatial_join["is_within_polygon"]
        & result_spatial_join["is_not_in_an_polygon_overlap"]
    ]
    logging.info("Filter Anndata object...")
    # The AnnData object is filtered to only contain the barcodes that are in non-overlapping polygon regions
    filtered_obs_mask = adata.obs_names.isin(barcodes_in_one_polygon["index"])
    filtered_adata = adata[filtered_obs_mask, :]

    # Add the results of the point spatial join to the Anndata object
    filtered_adata.obs = pd.merge(
        filtered_adata.obs,
        barcodes_in_one_polygon[
            [
                "index",
                "geometry",
                "id",
                "is_within_polygon",
                "is_not_in_an_polygon_overlap",
            ]
        ],
        left_index=True,
        right_index=True,
    )

    logging.info("Performing a gene-wise count summation for binned data...")
    # Group the data by unique nucleous IDs
    groupby_object = filtered_adata.obs.groupby(["id"], observed=True)

    # Extract the gene expression counts from the AnnData object
    counts = filtered_adata.X

    # Obtain the number of unique nuclei and the number of genes in the expression data
    N_groups = groupby_object.ngroups
    N_genes = counts.shape[1]

    # Initialize a sparse matrix to store the summed gene counts for each nucleus
    summed_counts = sparse.lil_matrix((N_groups, N_genes))

    # Lists to store the IDs of polygons and the current row index
    polygon_id = []
    row = 0

    # Iterate over each unique polygon to calculate the sum of gene counts.
    for polygons, idx_ in groupby_object.indices.items():
        summed_counts[row] = counts[idx_].sum(0)
        row += 1
        polygon_id.append(polygons)

    # Create and AnnData object from the summed count matrix
    summed_counts = summed_counts.tocsr()
    grouped_filtered_adata = anndata.AnnData(
        X=summed_counts,
        obs=pd.DataFrame(polygon_id, columns=["id"], index=polygon_id),
        var=filtered_adata.var,
    )

    # Store the area of each nucleus in the GeoDataframe
    gdf_segmentation["area"] = gdf_segmentation["geometry"].area

    logging.info("Saving masks to json file for use in R...")
    gdf_segmentation.to_file(
        os.path.join(output_dir, f"{sample_id}__masks.json"), driver="GeoJSON"
    )

    logging.info(
        "Saving summed count matrix based on nuclei segmentation as hdf5 file..."
    )
    write_10X_h5(
        grouped_filtered_adata,
        str(
            Path(
                output_dir,
                f"{sample_id}__filtered_feature_bc_matrix_w_seg.h5",
            )
        ),
    )


def main(args):
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(message)s")
    logging.info("Mapping barcodes")
    map_barcodes(
        sample_dir=args.sample_dir,
        gdf_barcodes_path=args.gdf_barcodes_path,
        gdf_segmentation_path=args.gdf_segmentation_path,
        output_dir=args.output_dir,
        sample_id=args.sample_id,
    )
    logging.info("Finished!")


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
