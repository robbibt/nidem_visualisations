# coding: utf-8

# NIDEM LiDAR extraction
#
# This script imports LiDAR .las files, extracts a 1% sample of LiDAR points from the point cloud,
# and exports these as .csv files of xyz points with a datestamp based on the time of each point's acquisition.
# This is then used as an input to `lidar_tidaltagging.py`, where non-inundated points are identified by
# selecting points located above the water's surface at the time each point was acquired.
#
# Date: October 2018
# Author: Robbi Bishop-Taylor, Steven Sagar, Leo Lymburner

import os
import math
import itertools
import pandas as pd
from pyproj import Proj, transform


#########
# Setup #
#########

# Read in setup parameters from file
os.chdir('C:/Users/u69654/Projects/nidem-GA/')
study_areas_df = pd.read_csv('lidar_study_areas.csv', index_col=0)
study_areas = study_areas_df.to_dict('index')

for name in ['Fraser', 'Whitsunday', 'Gladstone', 'Rockhampton', 'Isaac', 'Mackay', 'Fraser', 'Kaurumba']:

    # Read in study area details
    input_location = study_areas[name]['input_loc']

    # Create list of all lidar tiles that occur within given bounding box
    ul_lon, ul_lat = [float(coord) for coord in study_areas[name]['bbox_ul'].split(",")]
    br_lon, br_lat = [float(coord) for coord in study_areas[name]['bbox_br'].split(",")]

    # Convert lat-lon coordinates to local MGA zone
    mga_zone = study_areas[name]['mga_zone']
    proj_crs = {54: 'EPSG:28354', 55: 'EPSG:28355', 56: 'EPSG:28356'}[mga_zone]
    (ul_x, br_x), (ul_y, br_y) = transform(p1=Proj(init='EPSG:4326'), p2=Proj(init=proj_crs),
                                           x=[ul_lon, br_lon], y=[ul_lat, br_lat])

    # For each unique combination of 1x1km coordinates, produce file string
    all_combs = list(itertools.product(range(int(ul_x), int(br_x), 1000),
                                       range(int(ul_y), int(br_y), -1000)))
    loc_strings = set([str(math.floor(x / 1000)) + str(math.floor(y / 1000)) for x, y in all_combs])
    file_keys = [study_areas[name]['input_name'].format(loc) for loc in loc_strings]
    print(len(file_keys))

    ########################
    # Extract LAS into csv #
    ########################

    for file_key in file_keys:

        input_filename = '{}{}.las'.format(input_location, file_key)
        output_dir = os.path.normpath('{}/raw_data/validation'.format(os.getcwd()))
        output_filename = "{}_{}.csv".format(mga_zone, file_key)
        print('Downloading and extracting {}, MGA zone {}'.format(file_key, mga_zone))

        # If input file exists and not already processed, extract from LAS
        if os.path.isfile(input_filename) and not os.path.isfile('raw_data/validation/{}'.format(output_filename)):

            try:

                # Use lastools to convert data to text
                las2text_string = 'C:/Users/u69654/Desktop/lastools/LAStools/bin/las2txt.exe ' \
                                  '-i "{0}" ' \
                                  '-keep_random_fraction 0.01 ' \
                                  '-drop_classification 7 ' \
                                  '-odir "{1}" -o "temp.txt" ' \
                                  '-parse xyzcpt -sep comma'.format(input_filename, output_dir)
                os.system(las2text_string)

                # Read temporary file in and convert coordinates to lat/long
                points_df = pd.read_csv('{}/temp.txt'.format(output_dir), sep=',', header=None,
                                        names=['point_x', 'point_y', 'point_z', 'point_cat', 'point_path', 'point_time'])

                # Assign tide point
                tidepoint_lon, tidepoint_lat = [float(coord) for coord in study_areas[name]['tide_point'].split(",")]

                # Compute lon-lat coordinates for each point
                point_lon, point_lat = transform(p1=Proj(init=proj_crs), p2=Proj(init='EPSG:4326'),
                                                 x=points_df['point_x'].values, y=points_df['point_y'].values)

                # Assign tidepoint and point lon/lat to columns
                points_df['tidepoint_lon'] = tidepoint_lon
                points_df['tidepoint_lat'] = tidepoint_lat
                points_df['point_lon'] = point_lon
                points_df['point_lat'] = point_lat

                # Export to file
                points_df.to_csv('raw_data/validation/{}'.format(output_filename), index=False)

                # Clean up files
                points_df = None
                os.remove('raw_data/validation/temp.txt')

            except:

                print('Failed tile {}'.format(file_key))


#####################
# Export shapefiles #
#####################

# # Export Lidar sample sites to file
# lidar_val_sites = pd.DataFrame.from_dict({key:value["tide_point"].split(",") for key, value in study_areas.items()}, orient='index')
# lidar_val_sites.columns = ['longitude', 'latitude']
# lidar_val_sites.to_csv("lidar_validation_sites.csv")

# from shapely.geometry import mapping, Polygon
# import fiona
#
# # Import and plit columns
# lidar_val_sites = pd.DataFrame.from_dict(study_areas, orient='index')
# lidar_val_sites['bbox_ul'] = lidar_val_sites['bbox_ul'].str.split(',')
# lidar_val_sites['bbox_br'] = lidar_val_sites['bbox_br'].str.split(',')
#
# # Define a polygon feature geometry with one attribute
# schema = {'geometry': 'Polygon',
#           'properties': {'id': 'str'}}
#
# # Write a new Shapefile
# with fiona.open('lidar_validation_sites.shp', 'w', 'ESRI Shapefile', schema,
#                 crs={'init': 'epsg:4326', 'no_defs': True}) as shapefile:
#
#     for row in range(0, len(lidar_val_sites), 1):
#
#         # Here's an example Shapely geometry
#         (a, b), (c, d) = lidar_val_sites.iloc[row, 2:4]
#         point_list = list(itertools.product([float(a), float(c)], [float(b), float(d)]))
#
#         poly = Polygon([point_list[i] for i in [0, 2, 3, 1]])
#
#         # If there are multiple geometries, put the "for" loop here
#         shapefile.write({'geometry': mapping(poly),
#                          'properties': {'id': lidar_val_sites.index[row]}})
