import os
import math
import itertools
import pandas as pd
from pyproj import Proj, transform


#########
# Setup #
#########

# User input
name = "Whitsunday"
output_location = "C:/Users/u69654/Projects/lidar_processing/"

# Set up output location to read in setup parameters from file
study_areas_df = pd.read_csv('study_areas.csv', index_col=0)
study_areas = study_areas_df.to_dict('index')
input_location = study_areas[name]['input_loc']

# Create list of all lidar tiles that occur within given bounding box
ul_x, ul_y = study_areas[name]['bbox_ul'].split(", ")
br_x, br_y = study_areas[name]['bbox_br'].split(", ")
all_combs = list(itertools.product(range(int(float(ul_x)), int(float(br_x)), 1000),
                                   range(int(float(ul_y)), int(float(br_y)), -1000)))
loc_strings = set([str(math.floor(x / 1000)) + str(math.floor(y / 1000)) for x, y in all_combs])
file_keys = [study_areas[name]['input_name'].format(loc) for loc in loc_strings]


########################
# Extract LAS into csv #
########################

for file_key in file_keys:

    mga_zone = file_key[-12:-10]
    proj_crs = {'54': 'EPSG:28354', '55': 'EPSG:28355', '56': 'EPSG:28356'}[mga_zone]
    input_filename = "{}{}.las".format(input_location, file_key)
    output_dir = "{}output_data".format(output_location)
    output_filename = "{}_{}.csv".format(mga_zone, file_key)
    print("Downloading and extracting {}, MGA zone {}".format(file_key, mga_zone))

    try:

        # Use lastools to convert data to text
        las2text_string = 'C:/Users/u69654/Desktop/lastools/LAStools/bin/las2txt.exe ' \
                          '-i "{0}" ' \
                          '-keep_random_fraction 0.005 ' \
                          '-drop_classification 7 ' \
                          '-odir "{1}" -o "temp.txt" ' \
                          '-parse xyzcpt -sep comma'.format(input_filename, output_dir)
        os.system(las2text_string)

        # Read temporary file in and convert coordinates to lat/long
        points_df = pd.read_csv("{}/temp.txt".format(output_dir), sep=",", header=None,
                                names=["point_x", "point_y", "point_z", "point_cat", "point_path", "point_time"])

        # Assign tide point
        tidepoint_lat, tidepoint_lon = [float(coord) for coord in study_areas[name]['tide_point'].split(", ")]

        # Compute lon-lat coordinates for each point
        point_lon, point_lat = transform(p1=Proj(init=proj_crs), p2=Proj(init='EPSG:4326'),
                                         x=points_df['point_x'].values, y=points_df['point_y'].values)

        # Assign tidepoint and point lon/lat to columns
        points_df['tidepoint_lon'] = tile_lon
        points_df['tidepoint_lat'] = tile_lat
        points_df['point_lon'] = point_lon
        points_df['point_lat'] = point_lat

        # Export to file
        points_df.to_csv("{}/test_{}".format(output_dir, output_filename), index=False)

    except:

        print("Failed tile {}".format(raw_basename))
