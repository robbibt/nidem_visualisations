# coding: utf-8

# National Intertidal Digital Elevation Model (NIDEM)
# 
# This script generates Geoscience Australia's (GA) National Intertidal Digital Elevation Model (NIDEM) datasets,
# which provide continuous elevation data for Australia's intertidal zone. It initially imports layers from the DEA
# Intertidal Extents Model (ITEM v2.0) and median tidal elevations for each tidal interval, computes elevations at
# interval boundaries, extracts contours around each tidal interval, and then interpolates between these contours
# using TIN/Delaunay triangulation linear interpolation. This interpolation method preserves the tidal interval
# boundaries of ITEM v2.0. The notebook produces several output files:
# 
# 1. Contour line shapefiles (`NIDEM_contours_XXX.shp`) used for the interpolation. These datasets facilitate
#    re-analysis by allowing DEMs to be generated using alternative interpolation methods.
# 2. Mask rasters (`NIDEM_mask_XXX.tif` and `NIDEM_mask_XXX.netcdf`) that flag cells with elevations greater than
#    25m (value=1), less than -25m (value=2), and ITEM confidence NDWI standard deviation greater than 0.25 (value=3).
#    These masks were used to filter the output NIDEM layers.
# 3. Unfiltered NIDEM rasters (`NIDEM_unfiltered_XXX.tif` and `NIDEM_unfiltered_XXX.netcdf`) with elevations in metre
#    relative to Mean Sea Level.
# 4. Filtered NIDEM rasters (`NIDEM_dem_XXX.tif` and `NIDEM_dem_XXX.netcdf`) with elevations in metre units relative to
#    Mean Sea Level that have been cleaned by masking out cells included in the mask layers. This is the primary
#    output product, and is expected to be the default product used for most applications.
# 
# Date: June 2018
# Author: Robbi Bishop-Taylor


#####################################
# Load modules and define functions #
#####################################

import os
import glob
import itertools
import fiona
import numpy as np
import pandas as pd
import collections
import scipy.interpolate 
import skimage.measure
from osgeo import gdal
from scipy import ndimage as nd
from fiona.crs import from_epsg
from shapely.geometry import LineString, MultiLineString, mapping


def rasterize_vector(input_data, cols, rows, geo_transform, projection,
                     field=None, raster_path=None, array_dtype=gdal.GDT_UInt16):

    """
    Rasterize a vector file and return an array with values for cells that occur within the
    shapefile. Can be used to obtain a binary array (shapefile vs no shapefile), or can create
    an array containing values from a shapefile field. If 'raster_path' is provided, the
    resulting array can also be output as a geotiff raster.

    This function requires dimensions, projection data (in 'WKT' format) and geotransform info
    ('(upleft_x, x_size, x_rotation, upleft_y, y_rotation, y_size)') for the output array.

    Last modified: April 2018
    Author: Robbi Bishop-Taylor
    :param input_data:
        Input shapefile path or preloaded GDAL/OGR layer. This must be in the same
        projection system as the desired output raster (i.e. same as the 'projection'
        parameter below)

    :param cols:
        Desired width of output array in columns. This can be obtained from
        an existing array using '.shape[0]'

    :param rows:
        Desired height of output array in rows. This can be obtained from
        an existing array using '.shape[1]'

    :param geo_transform:
        Geotransform for output raster; e.g. "(upleft_x, x_size, x_rotation,
        upleft_y, y_rotation, y_size)"

    :param projection:
        Projection for output raster (in "WKT" format). This must be the same as the
        input shapefile's projection system (i.e. same projection as used by 'input_data')

    :param field:
        Shapefile field to rasterize values from. If None (default), this assigns a
        value of 1 to all array cells within the shapefile, and 0 to areas outside
        the shapefile

    :param raster_path:
        If a path is supplied, the resulting array will also be output as a geotiff raster.
        (defaults to None, which returns only the output array and does not write a file)

    :param array_dtype:
        Optionally set the dtype of the output array. This defaults to integers
        (gdal.GDT_UInt16), and should only be changed if rasterising float values from a
        shapefile field
    :return:
        A 'row x col' array containing values from vector (if Field is supplied), or binary
        values (1=shapefile data, 0=no shapefile)

    """

    # If input data is a string, import as shapefile layer
    if isinstance(input_data, str):
        # Open vector with gdal
        data_source = gdal.OpenEx(input_data, gdal.OF_VECTOR)
        input_data = data_source.GetLayer(0)

    # If raster path supplied, save rasterized file as a geotiff
    if raster_path:

        # Set up output raster
        print('Exporting raster to {}'.format(raster_path))
        driver = gdal.GetDriverByName('GTiff')
        target_ds = driver.Create(raster_path, cols, rows, 1, array_dtype)

    else:

        # If no raster path, create raster as memory object
        driver = gdal.GetDriverByName('MEM')  # In memory dataset
        target_ds = driver.Create('', cols, rows, 1, array_dtype)

    # Set geotransform and projection
    target_ds.SetGeoTransform(geo_transform)
    target_ds.SetProjection(projection)

    # Rasterize shapefile and extract array using field if supplied; else produce binary array
    if field:

        # Rasterise by taking attributes from supplied
        gdal.RasterizeLayer(target_ds, [1], input_data, options=["ATTRIBUTE=" + field])

    else:

        # Rasterise into binary raster (1=shapefile data, 0=no shapefile)
        gdal.RasterizeLayer(target_ds, [1], input_data)

        # Return array from raster
    band = target_ds.GetRasterBand(1)
    out_array = band.ReadAsArray()
    target_ds = None

    return out_array


def array_to_geotiff(fname, data, geo_transform, projection,
                     nodata_val=0, dtype=gdal.GDT_Float32):

    """
    Create a single band GeoTIFF file with data from an array.

    Because this works with simple arrays rather than xarray datasets from DEA, it requires
    geotransform info ("(upleft_x, x_size, x_rotation, upleft_y, y_rotation, y_size)") and
    projection data (in "WKT" format) for the output raster.

    Last modified: March 2018
    Author: Robbi Bishop-Taylor

    :param fname:
        Output geotiff file path including extension

    :param data:
        Input array to export as a geotiff

    :param geo_transform:
        Geotransform for output raster; e.g. "(upleft_x, x_size, x_rotation,
        upleft_y, y_rotation, y_size)"

    :param projection:
        Projection for output raster (in "WKT" format)

    :param nodata_val:
        Value to convert to nodata in the output raster; default 0

    :param dtype:
        Optionally set the dtype of the output raster; can be useful when exporting
        an array of float or integer values. Defaults to gdal.GDT_Float32

    """

    # Set up driver
    driver = gdal.GetDriverByName('GTiff')

    # Create raster of given size and projection
    rows, cols = data.shape
    dataset = driver.Create(fname, cols, rows, 1, dtype)
    dataset.SetGeoTransform(geo_transform)
    dataset.SetProjection(projection)

    # Write data to array and set nodata values
    band = dataset.GetRasterBand(1)
    band.WriteArray(data)
    band.SetNoDataValue(nodata_val)

    # Close file
    dataset = None


def reproject_to_template(input_raster, template_raster, output_raster, resolution=None,
                          resampling=gdal.GRA_Bilinear, nodata_val=0):
    """
    Reprojects a raster to match the extent, cell size, projection and dimensions of a template
    raster using GDAL. Optionally, can set custom resolution for output reprojected raster using
    'resolution'; this will affect raster dimensions/width/columns.

    Last modified: April 2018
    Author: Robbi Bishop-Taylor

    :param input_raster:
        Path to input geotiff raster to be reprojected (.tif)

    :param template_raster:
        Path to template geotiff raster (.tif) used to copy extent, projection etc

    :param output_raster:
        Output reprojected raster path with geotiff extension (.tif)

    :param resolution:
        Optionally set custom cell size for output reprojected raster; defaults to
        'None', or the cell size of template raster

    :param resampling:
        GDAL resampling method to use for reprojection; defaults to gdal.GRA_Bilinear

    :param nodata_val:
        Values in the output reprojected raster to set to nodata; defaults to 0

    :return:
        GDAL dataset for further analysis, and raster written to output_raster (if this
        dataset appears empty when loaded into a GIS, close the dataset like 'output_ds = None')

    """

    # Import raster to reproject
    print("Importing raster datasets")
    input_ds = gdal.Open(input_raster)
    input_proj = input_ds.GetProjection()
    input_geotrans = input_ds.GetGeoTransform()
    data_type = input_ds.GetRasterBand(1).DataType
    n_bands = input_ds.RasterCount

    # Import raster to use as template
    template_ds = gdal.Open(template_raster)
    template_proj = template_ds.GetProjection()
    template_geotrans = template_ds.GetGeoTransform()
    template_w = template_ds.RasterXSize
    template_h = template_ds.RasterYSize

    # Use custom resolution if supplied
    if resolution:
        template_geotrans[1] = float(resolution)
        template_geotrans[-1] = -float(resolution)

    # Create new output dataset to reproject into
    output_ds = gdal.GetDriverByName('Gtiff').Create(output_raster, template_w,
                                                     template_h, n_bands, data_type)
    output_ds.SetGeoTransform(template_geotrans)
    output_ds.SetProjection(template_proj)
    output_ds.GetRasterBand(1).SetNoDataValue(nodata_val)

    # Reproject raster into output dataset
    print("Reprojecting raster")
    gdal.ReprojectImage(input_ds, output_ds, input_proj, template_proj, resampling)

    # Close datasets
    input_ds = None
    template_ds = None

    print("Reprojected raster exported to {}".format(output_raster))
    return output_ds


###################
# Set up analysis #
###################

# Working directory
os.chdir('/g/data/r78/rt1527/nidem')

# Define path to ITEM offset product
item_offset_path = '/g/data2/v10/ITEM/offset_products'
item_relative_path = '/g/data2/v10/ITEM/rel_products'
item_conf_path = '/g/data2/v10/ITEM/conf_products'

# Define paths to external datasets
gbr30_raster = '/g/data/r78/rt1527/datasets/GBR30/02_ESRI_Raster/gbr30_ALL/gbr30_all'
ausbath09_raster = '/g/data/r78/rt1527/datasets/ausbath_09/ausbath_09_v4'
srtm30_raster = '/g/data/rr1/Elevation/1secSRTM_DEMs_v1.0/DEM/Mosaic/dem1sv1_0'
manually_included_shp = 'scratch/manual_mask.shp'

# Set ITEM polygon for analysis
polygon_ID = 96
# Issues with # 18, 131, 162, 172, 204, 220, 226, 238, 240, 301

# Print run details
print('Processing polygon {} from {}'.format(polygon_ID, item_offset_path))


###########################################
# ITEM interval boundary value extraction #
###########################################

# ITEM offset values represent the median tidal height for each tidal interval (Sagar et al. 2015,
# https://doi.org/10.1016/j.rse.2017.04.009). Because ITEM tidal intervals are linearly spaced by design,
# this code uses a simple linear model to compute new offset values for each interval boundary (e.g. the
# boundary between ITEM interval 1 and 2). This allows us to assign a more appropriate tidal height to the
# contours that divide the ITEM tidal intervals than would be possible through simply assigning median tidal
# heights to the downhill or uphill contours.

# Import ITEM offset values for each ITEM tidal interval, dividing by 1000 to give metre units
item_offsets = np.loadtxt('{}/elevation.txt'.format(item_offset_path), delimiter=',', dtype='str')
item_offsets = {int(key): [float(val) / 1000.0 for val in value.split(' ')] for (key, value) in item_offsets}
interval_offsets = item_offsets[polygon_ID]

# Create dataframe of offset values by ITEM interval
interval_offsets_df = pd.DataFrame({'item_interval': range(1, 10), 'offset': interval_offsets})

# Compute linear model and calculate ITEM offsets at the boundary of each ITEM interval to ensure that extracted
# contours are placed precisely on ITEM interval boundaries.
m, b = np.polyfit(interval_offsets_df['item_interval'], interval_offsets_df['offset'], 1)
interval_boundaries = np.arange(0.5, 10.5, 1.0)
contour_offsets = (m * interval_boundaries + b)

# Compute ITEM offset interval used to fill lowest class of ITEM relative layer
# (not used for interpolation, but ensures lowest contour is placed exactly on interval boundary)
interval_zero = (m * 0 + b)


#########################################
# Import and prepare ITEM offset raster #
#########################################

# Imports ITEM REL raster for given polygon, and use a lookup index array of offset values to classify into a new
# array of evenly-spaced ITEM offset values (in metre units relative to sea level) suitable for contour extraction.

# Import raster and extract shape, projection info and geotransform data
item_filename = glob.glob('{}/ITEM_REL_{}_*.tif'.format(item_relative_path, polygon_ID))[0]
item_ds = gdal.Open(item_filename)
item_array = item_ds.GetRasterBand(1).ReadAsArray() 
yrows, xcols = item_array.shape
prj = item_ds.GetProjection()
geotrans = item_ds.GetGeoTransform()
upleft_x, x_size, x_rotation, upleft_y, y_rotation, y_size = geotrans

# Assign nodata -6666 values to new class 10 prior to lookup classification
item_array[item_array == -6666] = 10
 
# Create lookup index array, and index by ITEM relative layer to classify ITEM classes into offset values
# (this method is resilient to ITEM layers with fewer than 9 classes)
lookup_index = np.array([interval_zero] + interval_offsets + [np.nan])
offset_array = lookup_index[item_array]

# Contours generated by `skimage.measure.find_contours` stop before the edge of nodata pixels. To prevent
# gaps from occurring between adjacent NIDEM tiles, the following steps 'fill' pixels directly on the boundary
# of two NIDEM tiles with the value of the nearest pixel with data. First, identify areas to be filled by dilating
# non-NaN pixels by two pixels (i.e. ensuring vertical, horizontal and diagonally adjacent pixels are filled):
dilated_mask = nd.morphology.binary_dilation(~np.isnan(offset_array), iterations=2)

# For every pixel, identify the indices of the nearest pixel with data (i.e. data pixels will return their own
# indices; nodata pixels will return the indices of the nearest data pixel). This output can be used to index back
# into the original array, returning a new array where data pixels remain the same, but every nodata pixel isn
# filled with the value of the nearest data pixel:
nearest_inds = nd.distance_transform_edt(input=np.isnan(offset_array), return_distances=False, return_indices=True)
offset_array = offset_array[tuple(nearest_inds)]

# Since we only want to fill pixels on the boundary of NIDEM tiles, set pixels outside the dilated area back to NaN:
offset_array[~dilated_mask] = np.nan


##########################################################
# Compute ITEM confidence and elevation/bathymetry masks #
##########################################################

# The following code applies a range of masks to remove pixels where elevation values are likely to be invalid:
# 
# 1. Non-coastal pixels with elevations greater than 25 m above MSL. This mask is computed using SRTM-derived
#    1 Second Digital Elevation Model data (http://pid.geoscience.gov.au/dataset/ga/69769).
# 2. Deep water pixels with bathymetry values less than -25 m below MSL. This mask is computed by identifying
#    any pixels that are < -25 m in both the national Australian Bathymetry and Topography Grid
#    (http://pid.geoscience.gov.au/dataset/ga/67703) and the  GBR30 bathymetry
#    (http://pid.geoscience.gov.au/dataset/ga/115066) datasets.
# 3. Pixels with high ITEM confidence NDWI standard deviation (i.e. areas where inundation patterns are not driven
#    by tidal influences). This mask is computed using ITEM v2.0 confidence layer data from DEA.
# 
# In a small number of locations where these masks remove valid data, a manual shapefile mask is used to preserve
# these pixels in the final datasets.

# Import ITEM confidence NDWI standard deviation array for polygon
conf_filename = glob.glob('{}/ITEM_STD_{}_*.tif'.format(item_conf_path, polygon_ID))[0]
conf_ds = gdal.Open(conf_filename)

# Reproject SRTM-derived 1 Second DEM to cell size and projection of NIDEM
srtm30_reproj = reproject_to_template(input_raster=srtm30_raster,
                                      template_raster=item_filename,
                                      output_raster='scratch/temp.tif',
                                      nodata_val=-9999)

# Reproject Australian Bathymetry and Topography Grid to cell size and projection of NIDEM
ausbath09_reproj = reproject_to_template(input_raster=ausbath09_raster,
                                         template_raster=item_filename,
                                         output_raster='scratch/temp.tif',
                                         nodata_val=-9999)

# Reproject GBR30 bathymetry to cell size and projection of NIDEM
gbr30_reproj = reproject_to_template(input_raster=gbr30_raster,
                                     template_raster=item_filename,
                                     output_raster='scratch/temp.tif',
                                     nodata_val=-9999)

# Import shapefile of areas to manually include in the output datasets and convert to 
# raster with same cell size and projection of NIDEM
manually_included = rasterize_vector(input_data=manually_included_shp,
                                     cols=xcols, rows=yrows,
                                     geo_transform=geotrans,
                                     projection=prj).astype(np.bool)

# Convert raster datasets to arrays
conf_array = conf_ds.GetRasterBand(1).ReadAsArray() 
srtm30_array = srtm30_reproj.GetRasterBand(1).ReadAsArray()
ausbath09_array = ausbath09_reproj.GetRasterBand(1).ReadAsArray()
gbr30_array = gbr30_reproj.GetRasterBand(1).ReadAsArray()

# Convert arrays to boolean masks:
#  For elevation: any elevations > 25 m in SRTM 30m DEM
#  For bathymetry: any depths < -25 m in both GBR30 and Ausbath09 bathymetry
#  For ITEM confidence: any cells with NDWI STD > 0.25
elev_mask = srtm30_array > 25
bathy_mask = (ausbath09_array < -25) & (gbr30_array < -25)
conf_mask = conf_array > 0.25

# Create a combined mask with -9999 nodata in unmasked areas and where:
#  1 = elevation mask
#  2 = bathymetry mask
#  3 = ITEM confidence mask
combined_mask = np.full(item_array.shape, -9999)
combined_mask[elev_mask] = 1
combined_mask[bathy_mask] = 2
combined_mask[conf_mask] = 3

# Set manually included pixels to -9999 to prevent masking
combined_mask[manually_included] = -9999


####################
# Extract contours #
####################

# Uses `skimage.measure.find_contours` to rapidly extract contour boundaries between ITEM tidal intervals, and
# assigns these contours with previously calculated elevation values. Contours are exported as line shapefiles
# to assist in subsequent assessment of output DEMs.

# Output dict to hold contours for each offset
contour_dict = collections.OrderedDict()

try:
    
    for contour_offset in contour_offsets:

        # Extract contours from array
        contours = skimage.measure.find_contours(offset_array, contour_offset)      
        print('Extracting contour {}'.format(contour_offset))
        
        # Iterate through each contour feature, remove NAs and fix coordinates
        contour_list = list()
        for contour in contours:

            # Convert index coordinates to spatial coordinates in-place
            contour[:, 0] = contour[:, 0] * float(y_size) + upleft_y + (float(y_size) / 2)
            contour[:, 1] = contour[:, 1] * float(x_size) + upleft_x + (float(x_size) / 2)        
            contour = np.insert(contour, 2, contour_offset, axis=1)

            # Remove contour points with NAs
            contour = contour[~np.isnan(contour).any(axis=1)] 
            contour_list.append(contour)

        # Add list of contour arrays to dict 
        contour_dict[contour_offset] = contour_list
        
except:
    
    print('Contour creation failed')

    
# Export contours to line shapefile to assist in evaluating DEMs
schema = {'geometry':  'MultiLineString', 
          'properties': {'elevation': 'float:9.2'}}

with fiona.open('output_data/contour/NIDEM_contours_{}.shp'.format(polygon_ID), 'w', 
                crs=from_epsg(3577),
                driver='ESRI Shapefile',
                schema=schema) as output:
    
    for elevation_value, contour_list in contour_dict.items():
        
        # Filter out contours with less than two points (i.e. non-lines)
        contour_list = [x for x in contour_list if len(x) > 1]
        
        # Create multiline string by first flipping coordinates then creating list of linestrings
        contour_linestrings = [LineString([(x, y) for (y, x, z) in contour_array]) 
                               for contour_array in contour_list]
        contour_multilinestring = MultiLineString(contour_linestrings)

        # Write output shapefile to file with elevation field        
        output.write({'properties': {'elevation': elevation_value},
                      'geometry': mapping(contour_multilinestring)})
        
# Chain and concatenate all arrays nested within array lists (i.e. individual collections of same 
# elevation contours) and dictionary entries (i.e. collections of all same-elevation contours)
all_contours = np.concatenate(list(itertools.chain.from_iterable(contour_dict.values())))
points = all_contours[:, 0:2]
values = all_contours[:, 2]


#######################################################################
# Interpolate contours using TIN/Delaunay triangulation interpolation #
#######################################################################

# Generates continuous elevation surfaces by interpolating between the extracted contours. This uses the linear
# method from `scipy.interpolate.griddata`, which computes a TIN/Delaunay triangulation of the input data using
# Qhull before performing linear barycentric interpolation on each triangle.
# 
# Because the lowest and highest ITEM intervals cannot be correctly interpolated as they have no lower or upper
# bounds, the filtered NIDEM is constrained to valid intertidal terrain (ITEM intervals 1-8).

# Calculate bottom right bounds of ITEM layer
bottomright_x = upleft_x + (x_size * xcols)
bottomright_y = upleft_y + (y_size * yrows)

# Create interpolation grid (from-to-by values in metre units)
grid_y, grid_x = np.mgrid[upleft_y:bottomright_y:1j * yrows, upleft_x:bottomright_x:1j * xcols]

# Interpolate between points onto grid. This uses the 'linear' method from 
# scipy.interpolate.griddata, which computes a TIN/Delaunay triangulation of the input 
# data with Qhull and performs linear barycentric interpolation on each triangle
print('Interpolating data for polygon {}'.format(polygon_ID))
interpolated_array = scipy.interpolate.griddata(points, values, (grid_y, grid_x), method='linear')

# Identify valid intertidal area by selecting pixels between the lowest and highest ITEM intervals
valid_intertidal_extent = np.where((item_array > 0) & (item_array < 9), 1, 0)

# Create filtered and unfiltered versions of NIDEM
nidem_unfiltered = np.where(valid_intertidal_extent, interpolated_array, -9999).astype(np.float32)
nidem_filtered = np.where(combined_mask > 0, -9999, nidem_unfiltered).astype(np.float32)


#######################
# Export geoTIFF data #
#######################

# NIDEM is exported as two DEMs: an unfiltered version, and a version filtered to remove > 25 m or < -25 m terrain
# and pixels with high ITEM confidence NDWI standard deviation.

# Export NIDEM mask as a GeoTIFF
print('Exporting NIDEM mask for polygon {}'.format(polygon_ID))
array_to_geotiff(fname='output_data/mask/geotiff/NIDEM_mask_{}.tif'.format(polygon_ID),
                 data=combined_mask.astype(int),
                 geo_transform=geotrans,
                 projection=prj,
                 nodata_val=-9999)

# Export unfiltered NIDEM as a GeoTIFF
print('Exporting unfiltered NIDEM for polygon {}'.format(polygon_ID))
array_to_geotiff(fname='output_data/geotiff/dem_unfiltered/NIDEM_unfiltered_{}.tif'.format(polygon_ID),
                 data=nidem_unfiltered,
                 geo_transform=geotrans,
                 projection=prj,
                 nodata_val=-9999)

# Export filtered NIDEM as a GeoTIFF
print('Exporting filtered NIDEM for polygon {}'.format(polygon_ID))
array_to_geotiff(fname='output_data/geotiff/dem/NIDEM_dem_{}.tif'.format(polygon_ID),
                 data=nidem_filtered,
                 geo_transform=geotrans,
                 projection=prj,
                 nodata_val=-9999)


######################
# Export NetCDF data #
######################

# os.remove('/g/data/r78/rt1527/nidem/output_data/netcdf/NIDEM_new.nc')

from datacube.model import Variable
from datacube.utils.geometry import Coordinate
from datacube.utils.geometry import CRS
from datacube.storage.storage import create_netcdf_storage_unit
from datacube.storage import netcdf_writer

# Compute coords
x_coords = netcdf_writer.netcdfy_coord(np.linspace(upleft_x + 12.5, bottomright_x - 12.5, num=xcols))
y_coords = netcdf_writer.netcdfy_coord(np.linspace(upleft_y - 12.5, bottomright_y + 12.5, num=yrows))

# Create new dataset
output_netcdf = create_netcdf_storage_unit(filename='/g/data/r78/rt1527/nidem/output_data/netcdf/NIDEM_new.nc',
                                           crs=CRS('EPSG:3577'),
                                           coordinates={'x': Coordinate(x_coords, 'metres'),
                                                        'y': Coordinate(y_coords, 'metres')},
                                           variables={'elevation': Variable(dtype=np.dtype('float32'),
                                                       nodata=-9999, dims=('y', 'x'), units='metres')},
                                           variable_params={'elevation': {}})

# Assign data and set variable attributes
output_netcdf['elevation'][:] = netcdf_writer.netcdfy_data(nidem_unfiltered)
output_netcdf['elevation'].valid_range = [-25.0, 25.0]
output_netcdf['elevation'].standard_name = 'height_above_mean_sea_level'
output_netcdf['elevation'].coverage_content_type = 'modelResult'
output_netcdf['elevation'].long_name = 'NIDEM filtered by ITEM confidence (< 0.25), ' \ 
                                       'bathymetry (> -25 m) and elevation (< 25 m)'

# Add global attributes
output_netcdf.title = 'National Intertidal Digital Elevation Model (NIDEM) 25m v 0.1.0'
output_netcdf.institution = 'Commonwealth of Australia (Geoscience Australia)'
output_netcdf.product_version = '0.1.0'
output_netcdf.license = 'CC BY Attribution 4.0 International License'
output_netcdf.time_coverage_start = '1986-01-01'
output_netcdf.time_coverage_end = '2016-10-31'
output_netcdf.cdm_data_type = 'Grid'
output_netcdf.contact = 'clientservices@ga.gov.au'
output_netcdf.publisher_email = 'earth.observation@ga.gov.au'
output_netcdf.keywords = 'Tidal, Topography, Landsat, Elevation, Intertidal, MSL, ITEM'
output_netcdf.summary = "The National Intertidal Digital Elevation Model (NIDEM) is a continental-scale dataset " \
                "providing a three-dimensional representation of Australia's exposed intertidal zone (the land " \
                "between the observed highest and lowest tide) at 25 metre resolution. The model is based on the " \
                "full 30 year archive of Landsat satellite data managed within the Digital Earth Australia (DEA) " \ 
                "platform that provides spatially and spectrally calibrated earth observation data to enable " \ 
                "time-series analysis on a per-pixel basis across the entire Australian continent. NIDEM builds " \ 
                "upon the improved tidal modelling framework of the Intertidal Extents Model v2.0 (ITEM), " \ 
                "allowing each satellite observation in the 30 year time series to be more accurately associated " \ 
                "with modelled tide heights from a multi-resolution global tidal model (OTPS TPX08). Using these " \ 
                "modelled tide heights and a spatially consistent and automated triangulated irregular network " \ 
                "(TIN) interpolation procedure, each pixel of exposed intertidal extent in ITEM was assigned " \
                "an absolute elevation in metre units relative to mean sea level."

# Close dataset
output_netcdf.close()


###############
# Close files #
###############

conf_ds = None
item_ds = None
gbr30_reproj = None
ausbath09_reproj = None
srtm30_reproj = None













from netCDF4 import Dataset
f = Dataset("/g/data/r78/rt1527/ITEM_REL_1_118.96_-19.95.nc")
f













def compliant_netcdf(filename, input_array, x_coords, y_coords, long_name, valid_range=None):



# import numpy as np
# from datetime import datetime
#
#
# # Create dataset
# os.remove('/g/data/r78/rt1527/NIDEM.nc')
# f_out = Dataset('/g/data/r78/rt1527/NIDEM.nc', 'w', format='NETCDF4')
#
# # Add global attributes
# f_out.date_created = datetime.now().isoformat()
# f_out.Conventions = 'CF-1.6, ACDD-1.3'
# f_out.title = 'National Intertidal Digital Elevation Model (NIDEM) 25m v 0.1.0'
# f_out.institution = 'Commonwealth of Australia (Geoscience Australia)'
# f_out.product_version = '0.1.0'
# f_out.license = 'CC BY Attribution 4.0 International License'
# f_out.acknowledgment = 'Insert list of acknowledgements'
# f_out.time_coverage_start = '1986-01-01'
# f_out.time_coverage_end = '2016-10-31'
# # f_out.geospatial_bounds =
# f_out.geospatial_bounds_crs = 'EPSG:4326'
# # f_out.geospatial_lat_min = -20.7385955147
# # f_out.geospatial_lat_max = -19.4964947991
# f_out.geospatial_lat_units = 'degrees_north'
# # f_out.geospatial_lon_min = 118.634108413
# # f_out.geospatial_lon_max = 119.437004501
# f_out.geospatial_lon_units = 'degrees_east'
# f_out.cdm_data_type = 'Grid'
# f_out.contact = 'clientservices@ga.gov.au'
# f_out.publisher_email = 'earth.observation@ga.gov.au'
# f_out.keywords = 'Tidal, Topography, Landsat, Elevation, Intertidal, MSL, ITEM'
# f_out.summary = "The National Intertidal Digital Elevation Model is a continental-scale dataset providing a " \
#             "three-dimensional representation of Australia's vast exposed intertidal zone at 25 metre " \
#             "resolution. This model was based on the full 30 year archive of Landsat satellite data managed " \
#             "within the Digital Earth Australia (DEA) platform that provides spatially and spectrally " \
#             "calibrated earth observation data to enable time-series analysis on a per-pixel basis across the " \
#             "entire Australian continent. The Intertidal Digital Elevation Model builds upon the improved " \
#             "tidal modelling framework of the Intertidal Extents Model v2.0 (ITEM 2.0), enabling us to more " \
#             "accurately associate each satellite observation in the 30 year time series with modelled tide " \
#             "heights from a global ocean tidal model. Using these modelled tide heights and a spatially " \
#             "consistent and automated triangulated irregular network (TIN) interpolation procedure, we " \
#             "assigned relative ITEM v2.0 tidal intervals with absolute elevations relative to mean sea level."
#
# # Create dimensions
# f_out.createDimension('y', yrows)
# f_out.createDimension('x', xcols)
# f_out.createDimension('time', 1)
#
# # Create variables
# x = f_out.createVariable(varname='x', datatype='f4', dimensions='x')
# y = f_out.createVariable(varname='y', datatype='f4', dimensions='y')
# elevation = f_out.createVariable(varname='elevation', datatype='f4', dimensions=('time', 'y', 'x'), fill_value=-9999.)
# time = f_out.createVariable(varname='time', datatype='i4', dimensions='time')
# crs = f_out.createVariable(varname='crs', datatype='i4', dimensions=())
#
# # Assign values to created variables
# x[:] = np.linspace(upleft_x + 12.5, bottomright_x - 12.5, num=xcols)
# y[:] = np.linspace(upleft_y - 12.5, bottomright_y + 12.5, num=yrows)
# time[:] = 1.47795840e+09
# elevation[0, :, :] = nidem_filtered
#
# # Set attributes
# crs.grid_mapping_name = 'albers_conical_equal_area'
# crs.standard_parallel = [-18., -36.]
# crs.longitude_of_central_meridian = 132.0
# crs.latitude_of_projection_origin = 0.0
# crs.false_easting = 0.0
# crs.false_northing = 0.0
# crs.long_name = 'GDA94 / Australian Albers'
# crs.semi_major_axis = 6378137.0
# crs.semi_minor_axis = 6356752.31414
# crs.inverse_flattening = 298.257222101
# crs.spatial_ref='PROJCS["GDA94 / Australian Albers",GEOGCS["GDA94",DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6283"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4283"]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["standard_parallel_1",-18],PARAMETER["standard_parallel_2",-36],PARAMETER["latitude_of_center",0],PARAMETER["longitude_of_center",132],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","3577"]]'
# crs.crs_wkt='PROJCS["GDA94 / Australian Albers",GEOGCS["GDA94",DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6283"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4283"]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["standard_parallel_1",-18],PARAMETER["standard_parallel_2",-36],PARAMETER["latitude_of_center",0],PARAMETER["longitude_of_center",132],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","3577"]]'
# crs.GeoTransform = geotrans
#
# # Set attributes
# x.standard_name = 'projection_x_coordinate'
# x.long_name = 'x coordinate of projection'
# x.units = 'meters'
#
# y.standard_name = 'projection_y_coordinate'
# y.long_name = 'y coordinate of projection'
# y.units = 'meters'
#
# time.axis = 'T'
# time.calendar = 'standard'
# time.long_name = 'Time, unix time-stamp'
# time.standard_name = 'time'
# time.units = 'seconds since 1970-01-01 00:00:00'
#
# elevation.standard_name = 'height_above_mean_sea_level'
# elevation.long_name = 'National Intertidal Digital Elevation Model (NIDEM)'
# elevation.grid_mapping = 'crs'
# elevation.coverage_content_type = 'modelResult'
# elevation.valid_range = [-25.0, 25.0]
# elevation.units = 'meters'
#
# # Close file
# f_out.close()






