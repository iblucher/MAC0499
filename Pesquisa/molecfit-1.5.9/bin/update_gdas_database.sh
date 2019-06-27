#!/bin/bash
#
#   Author: W. Kausch
#
#   Version 1.0:  18 July 2013
#   First version
#
#   Version 1.1: 28 October 2013
#   Usage of new GDAS file format including site selection (adding
#       input parameters <site>)
#
# Script to update GDAS database
#
#   Usage:    ./update_gdas_database.sh <months> <site>
#
#       format: <months> abbreviated 3 characters describing the month,
#               followed by two digits describing the year,
#               e.g. jan11 (January 2011),
#                    nov12 (November 2012),....
#               allowed values: jan, feb, mar, apr, may, jun, jul, aug, sep,
#                               oct, nov, dec
#
#               <lon>: longitude in decimal degrees in the range [-180,180]
#               <lat>: latitude  in decimal degrees in the range [-90,90]

# Examples:
#   ./update_gdas_database nov12 P
#       - adds November 2012 profiles for Cerro Paranal to gdas database file
#       gdas_profiles_C-70.4-24.6.tar.gz
#
#
#
txtbld=$(tput bold)       # Bold
txtrst=$(tput sgr0)       # Text reset

# bulk update:
# for y in $(seq 04 16); do \
#    for m in jan feb mar apr may jun jul aug sep oct nov dec; do \
#      echo ./bin/update_gdas_database.sh $(printf "%s%02d" $m $y) L; \
# done; done | NO_DOWNLOAD=1 BATCH_UPDATE=1 parallel --gnu
# cd data/downloads/<site>
# ls | GZIP='-9 --rsyncable -n' tar cfz ../"gdas_profiles_C<coordinates>.tar.gz" --files-from -

if [ ! $# -eq 2 ] && [ ! $# -eq 4 ]
then
    echo "Wrong number of input parameters: Two requested "
    echo "Usage: update_gdas_database.sh <month> <site>"
    echo "    "
    echo "being <month> = \"<mon><yr>\", e.g. apr08, jan12,...."
    echo "      <site> = P  Cerro Paranal "
    echo "             = L  La Silla"
    echo "             = A  Cerro Armazones"
    echo "             = C  Cajnantor"
    exit
fi

month=$1
site=$2


if [ "$site" != "P" -a "$site" != "L" -a "$site" != "A" -a "$site" != "C" -a "$site" != "custom" ];then
    echo "Observing site unknown."
    echo "Usage: update_gdas_database.sh <month> <site>"
    echo "    "
    echo "being <month> = \"<mon><yr>\", e.g. apr08, jan12,...."
    echo "      <site> = P  Cerro Paranal "
    echo "             = L  La Silla"
    echo "             = A  Cerro Armazones"
    echo "             = C  Cajnantor"
    exit
    exit
fi

if [ "$site" == "P" ];then
    lon=-70.4
    lat=-24.6
    echo "GDAS data for Cerro Paranal requested: Lon${lon}, Lat=${lat}"
    datdir="downloads/paranal"
fi
if [ "$site" == "L" ];then
    lon=-70.7
    lat=-29.3
    echo "GDAS data for La Silla requested: Lon${lon}, Lat=${lat}"
    datdir="downloads/lasilla"
fi
if [ "$site" == "C" ];then
    lon=-67.8
    lat=-23.0
    echo "GDAS data for Cajnantor requested: Lon${lon}, Lat=${lat}"
    datdir="downloads/cajnantor"
fi
if [ "$site" == "A" ];then
    lon=-70.2
    lat=-24.6
    echo "GDAS data for Cerro Armazones requested: Lon${lon}, Lat=${lat}"
    datdir="downloads/armazones"
fi
if [ "$site" == "custom" ];then
    lon=$3
    lat=$4
    echo "GDAS data for custom site requested: Lon${lon}, Lat=${lat}"
    datdir="downloads/${lon}${lat}"
fi


currdate=`date`

if [ ! "${month:0:3}" == "jan" -a ! "${month:0:3}" == "feb" -a \
     ! "${month:0:3}" == "mar" -a ! "${month:0:3}" == "apr" -a \
     ! "${month:0:3}" == "may" -a ! "${month:0:3}" == "jun" -a \
     ! "${month:0:3}" == "jul" -a ! "${month:0:3}" == "aug" -a \
     ! "${month:0:3}" == "sep" -a ! "${month:0:3}" == "oct" -a \
     ! "${month:0:3}" == "nov" -a ! "${month:0:3}" == "dec" ];then
    echo "Invalid month: ${month:0:3}"
    echo "    Usage: update_gdas_database.sh <month>"
    echo "    "
    echo "being <month> = \"<mon><yr>\", e.g. apr08, jan12,...."
    exit

fi

currdir=`pwd`

if [ ! -d data/ ]; then
    mkdir data/
fi
if [ ! -d data/downloads/ ]; then
    mkdir data/downloads/
fi
if [ ! -d data/downloads/gdas ]; then
    mkdir data/downloads/gdas
fi
if [ ! -d data/${datdir} ]; then
    mkdir data/${datdir}
fi
cd data/downloads/gdas

# for batch processing the wget to get the listing takes significant time
# which is not necessary when the files are already present locally
if [ -z "$NO_DOWNLOAD" ]; then
# work in tmpdir to allow for simple parallel processing
tmpdir="$(mktemp -d gdas_XXXXXX || exit 1)"
cd "$tmpdir"
# lock month to avoid duplicate download
(
    flock 9
    # to skip downloading already present files
    ls -1 ../gdas1.$month.* | xargs -n 1 ln -s 2>/dev/null
    echo "------------------------------------------------------------------------------------------------------"
    echo "${txtbld}Downloading: $month${txtrst}"
    wget -nc ftp://ftp.arl.noaa.gov/archives/gdas1/gdas1.$month.*
    echo "------------------------------------------------------------------------------------------------------"
    for f in gdas1.$month.*; do
        [ ! -L "$f" ] && mv "$f" ../
    done
    rm -f ./*
) 9> ${XDG_RUNTIME_DIR:-/tmp}/gdas.$month.lock
cd ..
rmdir "$tmpdir"
fi

cd $currdir

filelist="$(mktemp || exit 1)"

ls -1 data/downloads/gdas/gdas1.$month.w? >> "$filelist"

if [ ! -d data/${datdir} ]; then
    mkdir data/${datdir}
fi

if [ "$site" = "custom" ]; then
./bin/extract_gdas_profiles "$filelist" ${site} ${lon} ${lat} data/${datdir}
else
./bin/extract_gdas_profiles "$filelist" ${site} data/${datdir}
fi
rm "$filelist"

cd data/${datdir}
# update the tarball, better done in a final step when adding many month
if [ -z "$BATCH_UPDATE" ]; then
  if [ -e ./gdas_profiles_C${lon}${lat}.tar.gz ]; then
      tar zxf gdas_profiles_C${lon}${lat}.tar.gz
      rm gdas_profiles_C${lon}${lat}.tar.gz
  fi

  GZIP='-9 --rsyncable -n' tar cfz gdas_profiles_C${lon}${lat}.tar.gz *.gdas
fi

cd ..

cd $currdir
if [ -z "$BATCH_UPDATE" ]; then
  echo "GDAS database updated: $currdate Data from ${month} added" >> data/${datdir}/last_month.txt

  echo " "
  echo "GDAS database updated, month: ${month}"
  echo " "
  echo "Please copy tarball data/${datdir}/gdas_profiles_C${lon}${lat}.tar.gz  "
  echo "to your molecfit installation path <INST_DIR>/molecfit/bin/get_gdas_profile.sh"
  echo "(see Appendix 'Maintenance' in the Molecfit User Manual)."
fi
