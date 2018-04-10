import glob
import os
import pandas as pd
import numpy as np
import datetime as dt
from otps.predict_wrapper import predict_tide
from otps import TimePoint
from pytz import timezone


def gps_week(input_datetime):
    """
    Computes GPS week number since start of GPS time epoch (6 Jan 1980).
    Currently does not account for leap seconds (only affects 10-15 second
    window around midnight Saturday night/Sunday morning each week)

    :param input_datetime: Datetime object used to identify GPS week
    :return: GPS weeks since start of GPS epoch (6 Jan 1980)
    """

    # Identify GPS week from GPS epoch using floor division
    gps_epoch = dt.datetime(1980, 1, 6)
    delta = input_datetime - gps_epoch
    gps_week_num = int(np.floor(delta.total_seconds() / 86400 / 7))

    return gps_week_num


def gps_adj_utc(gps_adj, leap_seconds=10):
    """
    Converts between adjusted GPS time and UTC, returning a datetime object.
    This assumes adjusted GPS time has already had - 1 billion subtracted from it;
    if you have unadjusted GPS time instead, subtract 1 billion before inputting
    it into this function.

    :param gps_adj: Adjusted GPS time
    :param leap_seconds: Leap seconds since start of GPS epoch; default 10
    :return: Datetime object with converted time in UTC
    """

    # Identify UTC and GPS epochs and compute offset between them
    utc_epoch = dt.datetime(1970, 1, 1)
    gps_epoch = dt.datetime(1980, 1, 6)
    utc_offset = (gps_epoch - utc_epoch).total_seconds() - leap_seconds

    # Convert to unix time then UTC by adding 1 billion + UTC offset to GPS time
    unix_timestamp = utc_offset + (int(gps_adj) + 1000000000)
    utc_time = dt.datetime.utcfromtimestamp(unix_timestamp)

    # Set UTC timezone info
    utc_time = utc_time.replace(tzinfo=timezone('UTC'))

    return utc_time


def gps_sotw_utc(gps_sotw, reference_date, leap_seconds=10):
    """
    Computes UTC time from GPS Seconds of Week format time

    :param gps_sotw: GPS seconds-of-the-week value
    :param reference_date: Date used to compute current GPS week number
    :param leap_seconds: Leap seconds since start of GPS epoch; default 10
    :return: Datetime object with converted time in UTC
    """

    # First test if GPS seconds-of-week fall within 0 and 604800 seconds
    if 0 <= int(gps_sotw) <= dt.timedelta(days=7).total_seconds():

        # Identify UTC and GPS epochs and compute offset between them
        utc_epoch = dt.datetime(1970, 1, 1)
        gps_epoch = dt.datetime(1980, 1, 6)
        utc_offset = (gps_epoch - utc_epoch).total_seconds() - leap_seconds

        # Identify GPS week
        gps_week_num = gps_week(reference_date)

        # Compute difference between UTC epoch and GPS time, then add GPS week days
        unix_timestamp = utc_offset + int(gps_sotw)
        utc_basetime = dt.datetime.utcfromtimestamp(unix_timestamp)
        utc_time = utc_basetime + dt.timedelta(days=gps_week_num * 7)

        # Set UTC timezone info
        utc_time = utc_time.replace(tzinfo=timezone('UTC'))

        return utc_time

    else:
        print("GPS seconds-of-week must be between 0 and 604800 seconds")
        return None


# User input to read in setup parameters from file
name = "Whitsunday"
hdir = "/g/data/r78/rt1527/item_dem/validation_data/point_clouds/"

# Dict to convert MGA zones to EPSG
proj_dict = {'54': 'EPSG:28354', '55': 'EPSG:28355', '56': 'EPSG:28356'}

# Dict of study areas and files to process
study_areas_df = pd.read_csv('study_areas.csv', index_col=0)
study_areas = study_areas_df.to_dict('index')
point_files = [os.path.basename(file) for file in glob.glob("{}output_data/*{}*.txt".format(hdir, name))]

# List of dataframes
df_list = list()

# Iterate through each file
for input_file in point_files:

    # Read in with pandas
    points_df = pd.read_csv("{}output_data/{}".format(hdir, input_file), sep=",", header=None,
                            names=["point_x", "point_y", "point_z", "point_cat", "point_path", "point_time"])

    # Create dataframe
    df_list.append(points_df)

# Merge lists into single dataframe
points_df = pd.concat(df_list)


################
# Convert time #
################

# Convert GPS time to datetime, and round to nearest hour
points_df['point_time'] = points_df['point_time'].apply(lambda ts: gps_sotw_utc(ts, study_areas[name]['ref_date']))
points_df['point_timeagg'] = points_df['point_time'].dt.round('5min')  # 30min


#################
# Compute tides #
#################

# Group into unique times and locations, create TimePoints and model tides
grouped_series = points_df.groupby(['tidepoint_lat', 'tidepoint_lon', 'point_timeagg'])
grouped_series = grouped_series.apply(lambda row: TimePoint(lon=row.iloc[0]['tidepoint_lon'],
                                                            lat=row.iloc[0]['tidepoint_lat'],
                                                            timestamp=row.iloc[0]['point_timeagg']))

# Convert grouped data to dataframe and compute tides
grouped_df = grouped_series.to_frame(name="point_tidal")
grouped_df['point_tidal'] = [float(tp.tide_m) for tp in predict_tide(list(grouped_series))]  # 8 seconds

# Join back into main dataframe
points_df = points_df.join(grouped_df, on=['tidepoint_lat', 'tidepoint_lon', 'point_timeagg'], rsuffix="_test")
print(list(points_df.columns.values))

# Select output columns and export to file
points_df = points_df[['point_lon', 'point_lat', 'point_z', 'point_tidal',
                       'point_cat', 'point_path', 'point_time', 'point_timeagg']]
points_df.to_csv("{}output_data/output_points_whitsunday2.csv".format(hdir), index=False)
