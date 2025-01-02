import sys

import numpy as np
import tifffile
from scipy.sparse import csr_matrix

sys.path.append("libs")
import argparse
import logging
import os
import time
from argparse import ArgumentParser as AP
from multiprocessing import Pool
from os.path import abspath
from pathlib import Path

import geopandas as gpd
import pandas as pd
from libs.utils.helpers import create_polygon
from skimage.measure import label


def get_args():
    # Script description
    description = """Convert labeled image (TIF) to polygons (GeoDataFrame)"""

    # Add parser
    parser = AP(
        description=description, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--n_cores", type=int, help="Number of cores to use (default=1)", default=1
    )

    parser.add_argument(
        "--img_path",
        type=str,
        help="Path to segmentation results of cellpose (.npy)",
        required=True,
    )

    parser.add_argument(
        "--output_dir",
        type=str,
        help="Path to output folder to store generated files",
        required=False,
        default="",
    )

    parser.add_argument("--version", action="version", version="0.1.0")
    arg = parser.parse_args()
    arg.output_dir = abspath(arg.output_dir)

    if (arg.output_dir != "") & (not os.path.isdir(arg.output_dir)):
        os.mkdir(arg.output_dir)

    # Make Path
    arg.img_path = Path(arg.img_path)
    return arg


def compute_M(data):
    cols = np.arange(data.size)
    return csr_matrix((cols, (data.ravel(), cols)), shape=(data.max() + 1, data.size))


def get_indices_sparse(data):
    M = compute_M(data)
    return [np.unravel_index(row.data, data.shape) for row in M]


def convert_tif_to_poly(img_path: str, n_cores: int = 1):
    logging.info("Load labeled image...")
    # relabel to remove negative values
    labeled_img = label(tifffile.imread(img_path), background=0)
    print(len(np.unique(labeled_img)))

    # Remove background aka first item
    coords_of_masks = list(enumerate(get_indices_sparse(labeled_img)[1:]))

    logging.info("Create polygons for each mask...")
    pool = Pool(processes=n_cores)
    out = pool.map(create_polygon, coords_of_masks)

    logging.info("Combine dicts of all masks...")
    dict_items = {k: v for d in out for k, v in d.items()}

    print(len(dict_items))

    # Create geodataframe
    logging.info("Convert dict to pd.DataFrame...")
    df = pd.DataFrame.from_dict(dict_items, orient="index", columns=["geometry"])

    logging.info("Convert DataFrame to GeoDataFrame...")
    gdf = (
        gpd.GeoDataFrame(data=df, geometry="geometry")
        .reset_index()
        .rename(columns={"index": "id"})
    )

    logging.info("Set cell IDs...")
    # start ID with 1
    gdf.id = "ID_" + (gdf.id + 1).astype(str)

    return gdf


def main(args):
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(message)s")

    logging.info("Assign labels (cell or nuclei IDs) to barcodes...")
    poly_df = convert_tif_to_poly(img_path=args.img_path, n_cores=args.n_cores)

    logging.info("Save dataframe with barcode labels...")
    poly_df.to_parquet(Path(args.output_dir, f"{args.img_path.stem}__poly.parquet"))
    logging.info("Finished!")


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
