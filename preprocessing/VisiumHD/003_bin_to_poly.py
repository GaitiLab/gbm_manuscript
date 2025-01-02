# # Binning barcodes
import argparse
import json
import logging
import os
import time
from argparse import ArgumentParser as AP
from os.path import abspath
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from pandarallel import pandarallel
from shapely import Polygon


def get_args():
    # Script description
    description = """Labeling barcodes"""

    # Add parser
    parser = AP(
        description=description, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--tissue_positions_path",
        type=str,
        help="Tissue positions",
        required=True,
    )

    parser.add_argument(
        "--scale_factors_path", type=str, help="path to 'scalefactors_json.json'"
    )

    parser.add_argument("--radius_bin", type=int,
                        help="Radius of bin in um", default=1)

    parser.add_argument("--show_pb", type=bool,
                        help="Show progressbar", default=False)
    parser.add_argument(
        "--output_dir",
        type=str,
        help="Path to output folder to store generated files",
        required=False,
        default="",
    )
    parser.add_argument("--sample_id", type=str, help="Sample ID")
    parser.add_argument("--version", action="version", version="0.1.0")
    arg = parser.parse_args()
    arg.output_dir = abspath(arg.output_dir)

    if (arg.output_dir != "") & (not os.path.isdir(arg.output_dir)):
        os.mkdir(arg.output_dir)

    # Make Path
    arg.tissue_positions_path = Path(arg.tissue_positions_path)
    return arg


def bin_to_poly(
    tissue_positions_path: str,
    scale_factors_path: str,
    radius_bin,
    show_pb: bool = False,
):
    pandarallel.initialize(progress_bar=show_pb)

    # Convert scalefactor
    with open(scale_factors_path, "r") as file:
        scalefactors = json.load(file)
    um_per_px = scalefactors["microns_per_pixel"]

    radius_in_px = radius_bin / um_per_px

    # Load positions
    df_tissue_positions = pd.read_parquet(tissue_positions_path)

    # Set the index of the dataframe to the barcodes
    df_tissue_positions = df_tissue_positions.set_index("barcode")

    # Create an index in the dataframe to check joins
    df_tissue_positions["index"] = df_tissue_positions.index

    df_tissue_positions.loc[:, "row_min"] = np.floor(
        df_tissue_positions.loc[:, "pxl_row_in_fullres"] - radius_in_px
    )
    df_tissue_positions.loc[:, "row_max"] = np.ceil(
        df_tissue_positions.loc[:, "pxl_row_in_fullres"] + radius_in_px
    )

    df_tissue_positions.loc[:, "col_min"] = np.floor(
        df_tissue_positions.loc[:, "pxl_col_in_fullres"] - radius_in_px
    )
    df_tissue_positions.loc[:, "col_max"] = np.ceil(
        df_tissue_positions.loc[:, "pxl_col_in_fullres"] + radius_in_px
    )

    res = df_tissue_positions.parallel_apply(
        lambda row: Polygon(
            [
                (row.col_min, row.row_min),
                (row.col_max, row.row_min),
                (row.col_max, row.row_max),
                (row.col_min, row.row_max),
            ]
        ),
        axis=1,
    )
    return gpd.GeoDataFrame(df_tissue_positions, geometry=res)


def main(args):
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(message)s")

    logging.info("Assign labels (cell or nuclei IDs) to barcodes...")
    poly_df = bin_to_poly(
        tissue_positions_path=args.tissue_positions_path,
        show_pb=args.show_pb,
        scale_factors_path=args.scale_factors_path,
        radius_bin=args.radius_bin,
    )
    logging.info("Save dataframe with barcode labels...")

    poly_df.to_parquet(
        Path(
            args.output_dir,
            f"{args.sample_id}__tissue_positions__poly.parquet")
    )

    logging.info("Finished!")


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
