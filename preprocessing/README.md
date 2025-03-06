# Instructions

## Environment

The following package versions are used to run `map_barcodes.py` with Python version `3.11.10`

- anndata==0.10.9
- geopandas==1.0.1
- h5py==3.11.0
- numpy==1.26.4
- pandas==2.2.3
- scanpy=1.10.3
- scipy=1.14.1

## Run

```sh

# Polygons from bins (can be downloaded online, see publication)
gdf_barcodes_path="${segmentation_outs}/${sample_id}__tissue_positions__poly.parquet"

# Polygons from expanded nuclei (cell segmentation) segmentation (can be downloaded online, see publication)
gdf_segmentation_path="${segmentation_outs}/${sample_id}__stardist_masks__expanded__poly.parquet"

python3 "map_barcodes.py" \
    --output_dir ${output_dir} \
    --sample_dir ${sample_dir} \
    --gdf_barcodes_path ${gdf_barcodes_path} \
    --gdf_segmentation_path ${gdf_segmentation_path} \
    --sample_id ${sample_id}
```
