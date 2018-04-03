import glob
import os
import pandas as pd
import datetime
from pyproj import Proj, transform


def gps_to_epoch(gps_time):
    """
    Computes an approximate conversion between GPS time and UNIX epoch.
    Does not account for leap seconds at this stage, as this level of accuracy
    is not required for tidal modelling.
    """

    unix_epoch = int(315964782 + (gps_time + 1000000000))

    return unix_epoch


class TimePoint:
    """
    Stand-in for real TimePoint class
    """

    def __init__(self, lon, lat, dt):
        self.lon = lon
        self.lat = lat
        self.dt = dt


def predict_tide(tp):

    return float(tp.lon) / float(tp.lat)


# Home directory
hdir = "/g/data/r78/rt1527/item_dem/validation_data/point_clouds/"

# Dict to convert MGA zones to EPSG
proj_dict = {'54': 'EPSG:28354', '55': 'EPSG:28355', '56': 'EPSG:28356'}

# List of input file
point_files = [os.path.basename(file) for file in glob.glob("{}output_data/*.txt".format(hdir))]

# List of dataframes
df_list = list()

# Iterate through each file
for input_file in point_files:

    # Projection
    mga_zone = input_file[0:2]
    proj_crs = proj_dict[mga_zone]

    # Read in with pandas
    points_df = pd.read_csv("{}output_data/{}".format(hdir, input_file), sep=",", header=None,
                            names=["point_x", "point_y", "point_z", "point_cat", "point_path", "point_time"])

    # Compute coordinates for entire tile
    tile_x, tile_y = input_file[-25:-11].split("_")
    tile_lon, tile_lat = transform(p1=Proj(init=proj_crs), p2=Proj(init='EPSG:4326'), x=tile_x, y=tile_y)

    # Compute coordinates for each point
    point_lon, point_lat = transform(p1=Proj(init=proj_crs), p2=Proj(init='EPSG:4326'),
                                     x=points_df['point_x'].values, y=points_df['point_y'].values)

    # Assign dataset lon/lat to columns
    points_df['tile_lon'] = tile_lon
    points_df['tile_lat'] = tile_lat
    points_df['point_lon'] = point_lon
    points_df['point_lat'] = point_lat

    # Create dataframe
    df_list.append(points_df)

# Merge lists into single dataframe
points_df = pd.concat(df_list)


################
# Convert time #
################

# Convert GPS time to datetime, and round to nearest hour
points_df['point_time'] = points_df['point_time'].apply(lambda ts: datetime.datetime.utcfromtimestamp(gps_to_epoch(ts)))
points_df['point_timeagg'] = points_df['point_time'].dt.round('15min')  # 30min


#################
# Compute tides #
#################

# Group into unique times and locations, create TimePoints and model tides
grouped_series = points_df.groupby(['tile_lat', 'tile_lon', 'point_timeagg'])
grouped_series = grouped_series.apply(lambda row: TimePoint(lon=row.sample(1)['tile_lon'],
                                                            lat=row.sample(1)['tile_lat'],
                                                            dt=row.sample(1)['point_timeagg']))
grouped_series = grouped_series.apply(predict_tide)

# Convert grouped data to dataframe and join back into main dataframe
grouped_df = grouped_series.to_frame(name="point_tidal")
points_df = points_df.join(grouped_df, on=['tile_lat', 'tile_lon', 'point_timeagg'], rsuffix="_test")
print(list(points_df.columns.values))

# Select output columns and export to file
points_df = points_df[['point_lon', 'point_lat', 'point_z', 'point_tidal',
                       'point_cat', 'point_path', 'point_time', 'point_timeagg']]
points_df.to_csv("{}output_data/output_points.csv".format(hdir), index=False)
