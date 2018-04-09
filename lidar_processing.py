import boto3
import botocore
import os
import zipfile


# Home directory
hdir = "/g/data/r78/rt1527/item_dem/validation_data/point_clouds/"

# Dict to convert MGA zones to EPSG
proj_dict = {'54': 'EPSG:28354', '55': 'EPSG:28355', '56': 'EPSG:28356'}


####################
# Set up Amazon S3 #
####################

# Set up resource
s3 = boto3.resource('s3')

# Print out bucket names
for bucket in s3.buckets.all():
    print(bucket.name)

# Select QLD bucket
bucket = s3.Bucket('qld.elvis')


##########################
# Create spatial indexes #
##########################

# for mga_zone in proj_dict.keys():
#
#     # Get files information from bucket
#     files = bucket.objects.filter(Prefix="z{}/".format(mga_zone))
#
#     # Return only paths matching list
#     files_information = [file.key for file in files if file.key[-7:] == 'Las.zip']
#
#     # Extract coordinates and key value
#     coords = [(file_key[-25:-19], file_key[-18:-11], file_key) for file_key in files_information]
#
#     # Write to file
#     with open('{}output_data/lidar_index_{}.txt'.format(hdir, mga_zone), 'w') as fp:
#
#         # Write tile coordinates and file name/key to file
#         fp.write('\n'.join('%s %s %s' % x for x in coords))


#################################
# Download and extract from LAZ #
#################################

# Test download file from list (points = bottom left)
for i, file_key in enumerate(["z55/Cairns_2010_Prj_SW_330000_8181000_1K_Las.zip",
                              "z55/Cairns_2010_Prj_SW_330000_8180000_1K_Las.zip",
                              "z55/Cairns_2010_Prj_SW_330000_8179000_1K_Las.zip",
                              "z55/Cairns_2010_Prj_SW_330000_8178000_1K_Las.zip",
                              "z55/Cairns_2010_Prj_SW_330000_8177000_1K_Las.zip"]):

    mga_zone = file_key[1:3]
    raw_basename = file_key[4:-4]
    raw_filename = "{}raw_data/{}.zip".format(hdir, raw_basename)
    unzipped_filename = "{}raw_data/{}.laz".format(hdir, raw_basename)
    output_dir = "{}output_data".format(hdir)
    output_filename = "{}_{}.txt".format(mga_zone, raw_basename)
    print("Downloading and extracting {}, MGA zone {}".format(raw_basename, mga_zone))

    try:

        # Download file from S#
        s3.Bucket('qld.elvis').download_file(file_key, raw_filename)

        # Unzip zip file containing .LAZ point cloud
        with zipfile.ZipFile(raw_filename, "r") as zip_ref:

            zip_ref.extractall("{}raw_data/".format(hdir))

        # Use lastools to convert data to text
        las2text_string = 'C:/Users/u69654/Desktop/lastools/LAStools/bin/las2txt.exe ' \
                          '-i "{0}" ' \
                          '-drop_classification 7 ' \
                          '-keep_random_fraction 0.02 ' \
                          '-odir "{1}" -o "{2}" ' \
                          '-parse xyzcpt -sep comma'.format(unzipped_filename, output_dir, output_filename)

        os.system(las2text_string)

    except botocore.exceptions.ClientError as e:

        if e.response['Error']['Code'] == "404":
            print("The object does not exist.")

        else:
            raise
