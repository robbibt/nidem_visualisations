#!/bin/bash

# Low memory polygons, mem=8GB, walltime=0:45:00, jobfs=1GB
for polygon in 1 2

# # Higher memory polygons, mem=16GB, walltime=0:45:00, jobfs=1GB
# for polygon in 25 64 91 102 153 236 282

# # Very high memory and long wall time polygons, mem=64GB, walltime=4:00:00, jobfs=2GB
# for polygon in 8 178 220 301

# # Originally failed due to no data, mem=8GB, walltime=0:45:00, jobfs=1GB
# for polygon in 18 131 172 204 226 238 240

do

    PBS="#!/bin/bash\n\
    #PBS -N NIDEM_${polygon}\n\
	#PBS -P r78\n\
	#PBS -q express\n\
	#PBS -l walltime=0:30:00\n\
	#PBS -l mem=8GB\n\
	#PBS -l jobfs=1GB\n\
	#PBS -l ncpus=1\n\
	#PBS -l wd\n\
	module use /g/data/v10/public/modules/modulefiles\n\
    module load dea\n\
    python /g/data/r78/rt1527/nidem/NIDEM_generation.py $polygon"

    echo -e ${PBS} | qsub
    sleep 0.1
    echo "Submitting NIDEM polygon $polygon"

done
