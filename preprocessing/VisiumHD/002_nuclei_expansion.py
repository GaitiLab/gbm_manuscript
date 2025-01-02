import argparse
import json
import logging
import os
import time
from argparse import ArgumentParser as AP
from os.path import abspath
from pathlib import Path

import tifffile
from csbdeep.io import save_tiff_imagej_compatible
from PIL import Image
from skimage.segmentation import expand_labels


def get_args():
    # Script description
    description = """Nuclear Expansion"""

    # Add parser
    parser = AP(
        description=description, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--seg_path",
        type=str,
        help="Path to segmentation results of cellpose (.npy)",
        required=True,
    )
    parser.add_argument(
        "--distance",
        type=int,
        required=False,
        help="Distance to expand in µm",
        default=5,
    )

    parser.add_argument(
        "--scale_factors_path", type=str, help="path to 'scalefactors_json.json'"
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

    if (arg.output_dir != "") & (
        not os.path.isdir(arg.output_dir) & (Path(arg.output_dir).exists())
    ):
        os.mkdir(arg.output_dir)

    # Make Path
    arg.seg_path = Path(arg.seg_path)
    return arg


def nuclear_expansion(seg_path: str, scale_factors_path: str, distance: int = 5):
    """Nuclear Expansion

    Args:
        seg_path (str): Path to segmentation mask
        spot_diameter (float): Side length of the VisiumHD 2um squared in px, can be found in `scalefactors_json.json` as property 'spot_diameter_fullres'
        distance (int, optional): Distance to expand in um. Defaults to 5.

    Returns:
        arr: Expanded segmentation mask
    """

    # Load segmentation results from CellPose
    logging.info("Load TIFF-file...")
    segmentation_mask = tifffile.imread(Path(seg_path))

    # Convert scalefactor
    logging.info("Load scalefactors...")
    with open(scale_factors_path, "r") as file:
        scalefactors = json.load(file)
    um_per_px = scalefactors["microns_per_pixel"]

    logging.info("Compute distance to expand...")
    dist_in_px = distance / um_per_px

    # Expand masks by distance
    logging.info("Perform nuclei expansion...")
    expanded_segmentation_mask = expand_labels(
        segmentation_mask, distance=dist_in_px)
    return expanded_segmentation_mask


def main(args):
    Image.MAX_IMAGE_PIXELS = None

    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(message)s")

    logging.info("Expand masks...")
    expanded_segmentation_mask = nuclear_expansion(
        seg_path=args.seg_path,
        distance=args.distance,
        scale_factors_path=args.scale_factors_path,
    )

    logging.info("Save expanded masks...")
    save_tiff_imagej_compatible(
        Path(args.output_dir, f"{args.seg_path.stem}__expanded.tif"),
        expanded_segmentation_mask,
        axes="YX",
    )

    logging.info("Finished!")


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
