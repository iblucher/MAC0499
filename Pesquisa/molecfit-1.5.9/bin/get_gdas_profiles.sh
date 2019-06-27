#!/bin/bash
#
# Simple script to extract GDAS profiles from an archive 'gdas_profiles.tar.gz'.
# Data are currently available only from 2004 12 01 - 2011 01 31.
# (see http://ready.arl.noaa.gov/data/archives/gdas1/missing.txt for a list of missing data)
# (based on extract_grib.sh by ESO)
#
# NOTE: Currently, the tarball contains >24000 single files! Hence, the archive
#       must not be extracted. We extract the required profiles directly from
#       it.
#
# Input format has to be yyyymmdd hh. Ex: May 01 2009 at 06h15 you have to write:
#     get_gdas_profiles.sh <basedir> 20090501 06.25
#
# basedir = installation directory containing bin/ data/ config/ et al. as subdirs.
#
#
# Version 1.0
# Author: Wolfgang Kausch
# In final form: 28 January 2011
# Last modification: 05th April 2011
#
# Version 1.1
#   - Updating GDAS archive --> ENDDATE = 2011 01 31
#   - minor bugfixes
#
# Version 1.2
#   - Updating GDAS archive --> ENDDATE = 2011 02 28
#
# Version 1.3
#   - Updating GDAS archive --> ENDDATE = 2011 03 31
#
# Version 1.4
#   - Updating GDAS archive --> ENDDATE = 2012 01 31
#
# Version 1.5
#   - Updating GDAS archive --> ENDDATE = 2012 08 31
#
# Version 1.6
#   - Updating GDAS archive --> ENDDATE = 2012 12 31
#
# Version 1.7
#   - Updating GDAS archive --> ENDDATE = 2013 02 28
#
# Version 1.8
#   - Updating GDAS archive --> ENDDATE = 2013 04 30
#
# Version 2.0
#   - Change of tarballfilename including coordinates: gdas_profiles_C<lon><lat>.tar.gz
#   - Updating GDAS archive --> ENDDATE = 2013 09 30
#
# Input parameters:
#     1st parameter: basedir of SM-02; MUST contain directories /bin, /data/profiles/grib
#     2nd parameter: date in the form YYYYMMDD
#     3rd parameter: hour in the form HH.hh        (i.e., hour and decimal hour)
#                                                  MUST be 5 characters in total!!
#     4th parameter: longitude in format <sign><degree> with one(!) decimal place
#     5th parameter: latitude in format <sign><degree> with one(!) decimal place
#
# Output files of this script will be in format yyyy-mm-ddThh.gdas and will be put in the directory dir
#


#
# DO NOT MODIFY ANYTHING ELSE BELOW
#
if [ $# -le 2 ]
then
    echo "Missing input. Stop. "
    echo "    Usage: get_gdas_profiles.sh <dir-of-data> <date YYYYMMDD> <hour HH.hh>"
    exit
fi

basedir=$1
bindir=$basedir/bin
datadir=$basedir/data/profiles/gdas
outputdir=$basedir/data/profiles/grib

date=$2
hour=$3
lon=$4
lat=$5

archive="$datadir/gdas_profiles_C${lon}${lat}.tar.gz"

# Extracting date
year=${date:0:4}
month=${date:4:2}
day=${date:6:2}

year_begin=${date:0:4}
month_begin=${date:4:2}
day_begin=${date:6:2}

year_end=${date:0:4}
month_end=${date:4:2}
day_end=${date:6:2}

date_begin=$year_begin$month_begin$day_begin
date_end=$year_end$month_end$day_end

hh=${hour:0:2}
hh=$(echo $hour |cut -d'.' -f1) # skipping decimal places from hour variable


if [ $hh -lt 06 ]
    then
    hour_begin="00"
    hour_end="06"
fi

if [ $hh -ge 06 ] && [ $hh -lt 12 ]
    then
    hour_begin="06"
    hour_end="12"
fi

if [ $hh -ge 12 ] && [ $hh -lt 18 ]
    then
    hour_begin="12"
    hour_end="18"
fi

# Checking changes of day/month/year
if [ $hh -ge 18 ]
    then
    hour_begin="18"
    hour_end="00"
    juldate=`echo  ${year_begin} ${month_begin} ${day_begin} | $bindir/Gregorian2JD`
    let juldate=juldate+1
    new_date=`echo $juldate | $bindir/JD2Gregorian`
    year_end=${new_date:9:4}
    month_end=${new_date:14:2}
    day_end=${new_date:27:2}

    if [ $day_end -le 9 ]
    then
        day_end=0${new_date:28:1}
    fi
fi

date_begin=$year_begin$month_begin$day_begin
date_end=$year_end$month_end$day_end
echo " "
echo "############################################################"
echo "#                                                          #"
echo "# GETTING GDAS FILES FROM ${year_begin} ${month_begin} ${day_begin} ${hour_begin}h to ${year_end} ${month_end} ${day_end} ${hour_end}h #"
echo "#                                                          #"
echo "############################################################"
echo " "
file1="C${lon}${lat}D${year_begin}-${month_begin}-${day_begin}T${hour_begin}.gdas"
file2="C${lon}${lat}D${year_end}-${month_end}-${day_end}T${hour_end}.gdas"
echo "file1=$file1"
echo "file2=$file2"

echo "tar zxf $archive $file1  2>/dev/null"
echo "tar zxf $archive $file2  2>/dev/null"
tar zxf $archive $file1  2>/dev/null
tar zxf $archive $file2  2>/dev/null

mkdir -p "$outputdir"
mv "$file1" "$outputdir"/ 2>/dev/null
mv "$file2" "$outputdir"/ 2>/dev/null

if [ -e $outputdir/$file1 -a -e $outputdir/$file2 ]; then
    echo -n "Profiles '$file1' and '$file2' found, "
    echo "copied to $outputdir/."
else
    echo "WARNING: No GDAS profiles found in archive!"
fi
echo " "



