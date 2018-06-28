# National Intertidal Digital Elevation Model (NIDEM)

Date: June 2018

Author: Robbi Bishop-Taylor, Steven Sagar, Leo Lymburner

## Description

This script generates Geoscience Australia's (GA) National Intertidal Digital Elevation Model (NIDEM) datasets, which provide continuous elevation data for Australia's intertidal zone. It initially imports layers from the DEA Intertidal Extents Model (ITEM v2.0) and median tidal elevations for each tidal interval, computes elevations at interval boundaries, extracts contours around each tidal interval, and then interpolates between these contours using TIN/Delaunay triangulation linear interpolation. This interpolation method preserves the tidal interval boundaries of ITEM v2.0.

## To generate NIDEM

To generate NIDEM datasets:

 1. Set the locations to input datasets in the `NIDEM_configuration.ini` configuration .ini file
 2. On the NCI, run the `NIDEM_pbs_submit.sh` shell script which iterates through a set of ITEM polygon tile IDs in parallel. This script calls `NIDEM_generation.py`, which conducts the actual analysis.

## Output files

NIDEM consists of several output files:

1. Contour line shapefiles (`NIDEM_contours_XXX.shp`, with `XXX` representing the ITEM polygon tile ID) used for the interpolation. These datasets facilitate re-analysis by allowing DEMs to be generated using alternative interpolation methods.
2. Mask rasters (`NIDEM_mask_XXX.tif`) that flag cells with elevations greater than 25 m (value = 1), less than -25 m (value = 2), and ITEM confidence NDWI standard deviation greater than 0.25 (value = 3). These masks were used to filter the output NIDEM layers.
3. Unfiltered NIDEM rasters (`NIDEM_unfiltered_XXX.tif`) with elevations in metres relative to Mean Sea Level.
4. Filtered NIDEM rasters (`NIDEM_dem_XXX.tif`) with elevations in metre units relative to Mean Sea Level that have been cleaned by masking out cells included in the mask layers. This is the primary output product, and is expected to be the default product used for most applications.

The mask, unfiltered and filtered NIDEM products are also exported as a combined NetCDF dataset (`NIDEM_XXX.nc`).
