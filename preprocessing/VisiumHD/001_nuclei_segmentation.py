import argparse
import logging
import os
import time
from argparse import ArgumentParser as AP
from os.path import abspath
from pathlib import Path

import geopandas as gpd

# import libs.vis.plot as vis
from csbdeep.io import save_tiff_imagej_compatible
from csbdeep.utils import normalize

# from matplotlib.colors import ListedColormap
from shapely.geometry import Polygon
from stardist.models import StarDist2D
from tifffile import imread


def get_args():
    # Script description
    description = """Nuclei segmentation with StarDist"""

    # Add parser
    parser = AP(
        description=description, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-f", "--image_path", type=str, help="Path to high-res microscope image"
    )
    parser.add_argument(
        "-m",
        "--model_dir",
        type=str,
        help="Path directory with models, models are in subdirectories",
    )
    parser.add_argument(
        "-mn", "--model_name", type=str, help="Name of model, {model_dir}/{model_name}"
    )
    parser.add_argument(
        "-of",
        "--output_filename",
        type=str,
        help="Filename for created output, image and geodataframe (nuclei segmentation mask) (default='nuclei_segmentation')",
        default="nuclei_segmentation",
    )
    parser.add_argument(
        "-min_p",
        "--min_percentile",
        type=int,
        default=5,
        help="min percentile for percentile normalization of the image (default=5)",
    )
    parser.add_argument(
        "-max_p",
        "--max_percentile",
        type=int,
        default=95,
        help="max percentile for percentile normalization of the image (default=95)",
    )
    parser.add_argument(
        "-p_thresh", "--prob_thresh", type=float, default=0.01, help=" (default=0.01)"
    )
    parser.add_argument(
        "-nms_thresh",
        "--nms_thresh",
        type=float,
        default=0.001,
        help=" (default=0.001)",
    )
    parser.add_argument(
        "-x_min",
        "--x_min",
        type=int,
        default=0,
        help="Min x coordinate bounding box for ROI",
    )
    parser.add_argument(
        "-x_max",
        "--x_max",
        type=int,
        default=0,
        help="Max x coordinate bounding box for ROI",
    )
    parser.add_argument(
        "-y_min",
        "--y_min",
        type=int,
        default=0,
        help="Min y coordinate bounding box for ROI",
    )
    parser.add_argument(
        "-y_max",
        "--y_max",
        type=int,
        default=0,
        help="Max y coordinate bounding box for ROI",
    )
    parser.add_argument(
        "--use_bounding_box",
        help="Use bounding box (1) or all (0) for plotting, if ROI the specify bounding box (x_min, y_min, x_max, y_max) ",
        default=0,
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

    if (
        (arg.output_dir != "")
        & (not os.path.isdir(arg.output_dir))
        & (not Path(arg.output_dir).exists())
    ):
        os.mkdir(arg.output_dir)

    return arg


def nuclei_segmentation(
    image_path,
    model_name,
    model_dir,
    min_percentile,
    max_percentile,
    nms_thresh,
    prob_thresh,
):
    # Load the image file
    logging.info("Loading high-res image...")
    img = imread(image_path)

    # Load the pretrained model
    logging.info("Loading pretrained model...")
    # model = StarDist2D.from_pretrained('2D_versatile_he')
    model = StarDist2D(None, name=model_name, basedir=model_dir)

    # Percentile normalization of the image
    # Adjust min_percentile and max_percentile as needed
    logging.info("Normalize image using percentile normalization...")
    img = normalize(img, min_percentile, max_percentile)

    logging.info("Creating nuclei segmentation mask...")
    # Predict cell nuclei using the normalized image
    # Adjust nms_thresh and prob_thresh as needed
    # This step's execution time is long
    # Divides input image into blocks, processing htem individually using 'predict_instances' and then combines the results
    labels, polys = model.predict_instances_big(
        img,
        axes="YXC",
        block_size=4096,
        prob_thresh=prob_thresh,
        nms_thresh=nms_thresh,
        min_overlap=128,
        context=128,
        normalizer=None,
        n_tiles=(4, 4, 1),
    )
    # ONLY FOR TESTING
    # labels, polys = model.predict_instances(img)

    logging.info("Converting StarDist results into a Geodataframe...")
    # Creating a list to store Polygon geometries
    geometries = []

    # Iterating through each nuclei in the 'polys' DataFrame
    for nuclei in range(len(polys["coord"])):
        # Extracting coordinates for the current nuclei and converting them to (y, x) format
        coords = [
            (y, x) for x, y in zip(polys["coord"][nuclei][0], polys["coord"][nuclei][1])
        ]

        # Creating a Polygon geometry from the coordinates
        geometries.append(Polygon(coords))

    # Creating a GeoDataFrame using the Polygon geometries
    gdf = gpd.GeoDataFrame(geometry=geometries)
    gdf["id"] = [f"ID_{i+1}" for i, _ in enumerate(gdf.index)]

    return labels, gdf, img


def main(args):
    # Setup logging
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(message)s")

    labels, gdf, img = nuclei_segmentation(
        image_path=args.image_path,
        model_name=args.model_name,
        model_dir=args.model_dir,
        min_percentile=args.min_percentile,
        max_percentile=args.max_percentile,
        nms_thresh=args.nms_thresh,
        prob_thresh=args.prob_thresh,
    )

    save_tiff_imagej_compatible(
        Path(args.output_dir, f"{args.output_filename}_masks.tif"), labels, axes="YX"
    )

    logging.info(
        "Saving nuclei segmentation (geoDataFrame) as parquet file...")
    gdf.to_parquet(os.path.join(args.output_dir,
                                f"{args.output_filename}.parquet"))


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
