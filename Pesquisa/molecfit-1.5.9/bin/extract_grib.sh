#!/bin/bash
set -e
#
# Simple script to download parameters at all height levels from the nomad3.ncep.noaa server.
# Data are available only from 2006 02 01.
#
# Input format has to be yyyymmdd hh. Ex: May 01 2009 at 06h15 you have to write:
#     extract_grib.sh 20090501 06.25
#
# Version 1.1
# Author: Alexandre Gallenne <agallenn@eso.org>
# In final form: 04 december 2009
# Last modification: 10 december 2009
#
# Version 1.2
# Modified by H. Sana to save file names  in 3rd argument
#
# Version 1.3
# Modified by A. Smette:
#     -  longitude [between 0 and 360]  in 3rd argument
#     -  latitude  [between -90 and 90] in 4th argument
#     -  directory where gdas files are saved in 5th argument
#     -  filename of the file that contains the names and hour of gdas files in 6th argument
#
# Version 2.0
# Last modification: 28th February 2011
# Modified by Innsbruck ESO in-kind Team:
#     - Data rows are reversed due to changes of format on nomad3 server after
#       15-12-2009 --> Script restricted to dates after 15-12-2009
#     - Adding argument #6 / #7 for output options
#     - Bugfix: hour 18:00 < 00:00 caused problem calculating second gdas name
#     - Bugfix: Change 2009 --> 2010 caused problem calculating second gdas name
#
#
# Version 2.1
# Last modifications: 28 October 2013
# Modified by Innsbruck ESO in-kind Team:
#     - Change to use new GDAS filenames (adding param #6)
#
# This script use shell, fortran and perl codes.
#
# Output files of this script will be in format yyyy-mm-ddThh.gdas and will be put in the directory dir of you choice.
#
##################################################################################################

# Path for cnvgrib, wgrib2, wgrib, Gregorian2JD, JD2Gregorian and the output directory
#cnv="/Users/alexandregallenne/Desktop/Projet/grid/cnvgrib-1.1.9/cnvgrib"

#  Input parameters:
#     1st parameter: installation directory of grib
#     2nd parameter: date in the form YYYYMMDD
#     3rd parameter: hour in the form HH.hh        (i.e., hour and decimal hour)
#                                                  MUST be 5 characters in total!!
#     4th parameter: longitude in decimal degrees in the range [0,360]
#     5th parameter: latitude  in decimal degrees in the range [-90,90]
#     6th parameter: longitude in decimal degrees in the range [-180,180] (for filename use only)
#     7th parameter: optional, directory where to save the output files
#     8th parameter: optinal, name of the file that contains the output files and corresponding hour
#

if [ $2 -le 20091216 ]
then
    echo "Script 'bin/extract_grib.sh' not applicable for dates before 17-12-2009, exiting."
    exit
fi

tmpdir=$(mktemp -d molec_XXXXXX || exit 1)
trap "rm -rf $tmpdir" EXIT

extractgribdir="$1"
file="$2"
hour="$3"
longitude="$4"
latitude="$5"
long_fn="$6"
dir="./"
if [ $# -ge 7 ]
then
    dir="$7""/"
fi
mkdir -p $dir

if [ $# -ge 8 ]
then
    outputname="$8"
fi

tac="tac"
# macos has tail -r instead of tac
if ! command -v tac >/dev/null 2>/dev/null; then
  tac="tail -r"
fi


grib2="$extractgribdir""/wgrib2"
grib1="$extractgribdir""/wgrib"
GtoJD="$extractgribdir""/Gregorian2JD"
JDtoG="$extractgribdir""/JD2Gregorian"

hh=${hour:0:2}
hh=$(echo $hour |cut -d'.' -f1)

##### Files selection #####

if [ $file -lt 20060201 ]
then
    echo "no archives for this date."
    exit
fi

if [ $hh -lt 06 ]
    then
    file1=gdas"$file"00f00
    file2=gdas"$file"06f00
fi

if [ $hh -ge 06 ] && [ $hh -lt 12 ]
    then
    file1=gdas"$file"06f00
    file2=gdas"$file"12f00
fi

if [ $hh -ge 12 ] && [ $hh -lt 18 ]
    then
    file1=gdas"$file"12f00
    file2=gdas"$file"18f00
fi

var=0
if [ $hh -ge 18 ]
    then
    var=1

    file1=gdas"$file"18f00
    juldate=`echo  ${file:0:4} ${file:4:2} ${file:6:2} | $GtoJD`

    new_juldate=$(($juldate + 1))
    new_date=`echo $new_juldate | $JDtoG`

    file4=20${new_date:11:2}${new_date:14:2}
    if [ ${new_date:27:2} -lt 10 ]    # adding zero if day < 10
    then
        file2=gdas20"${new_date:11:2}""${new_date:14:2}"0"${new_date:28:1}"00f00
    else
        file2=gdas20"${new_date:11:2}""${new_date:14:2}""${new_date:27:2}"00f00
    fi
fi

if [ $file -ge 20081031 ]
then
    file2=$file2.grib2
    if [ $file -ge 20081101 ]
    then
        file1=$file1.grib2
    fi
fi

if [ $var == 0 ]
then
    file3=${file:0:6}
    file4=${file:0:6}
else
    file3=${file:0:6}
fi

echo "###################################################"
echo "#                                                 #"
echo "# GET FILES FROM ${file1:4:4} ${file1:8:2} ${file1:10:2} ${file1:12:2}h to ${file2:4:4} ${file2:8:2} ${file2:10:2} ${file2:12:2}h #"
echo "#                                                 #"
echo "###################################################"

##### Temporary output files #####
height="HGT"
temperature="TMP"
humidity="RH"

outfile_temperature_1="$tmpdir/temperature_1"
outfile_humidity_1="$tmpdir/humidity_1"
outfile_height_1="$tmpdir/height_1"
outfile_temperature_2="$tmpdir/temperature_2"
outfile_humidity_2="$tmpdir/humidity_2"
outfile_height_2="$tmpdir/height_2"

##url="http://nomad3.ncep.noaa.gov/pub/gdas/rotating/200904/$file"
url1="http://nomad3.ncep.noaa.gov/pub/gdas/rotating/$file3/$file1"
url2="http://nomad3.ncep.noaa.gov/pub/gdas/rotating/$file4/$file2"

if [ -n "$outputname" ]; then
    echo $file1 ${file1:12:2} >  $outputname
    echo $file2 ${file2:12:2} >> $outputname
fi

set +e

echo $url1
$extractgribdir/get_inv.pl "$url1.inv" > $tmpdir/my_inv1
res1=$?
echo $url2
$extractgribdir/get_inv.pl "$url2.inv" > $tmpdir/my_inv2
res2=$?

# url is different > 201307
if [ $res1 -ne 0 ] || [ ! -s $tmpdir/my_inv1 ]; then
    url1="http://nomad3.ncep.noaa.gov/pub/gdas/rotating/$file1"
#    url1="http://140.90.198.158/pub/gdas/rotating/$file1"

    echo $url1
    $extractgribdir/get_inv.pl "$url1.inv" > $tmpdir/my_inv1
    res1=$?
    if [ $res1 -ne 0 ] || [ ! -s $tmpdir/my_inv1 ]; then
        echo "no $file1 available yet\n"
        exit
    fi
fi
if [ $res2 -ne 0 ] || [ ! -s $tmpdir/my_inv2 ]; then
    url2="http://nomad3.ncep.noaa.gov/pub/gdas/rotating/$file2"

    echo $url2
    $extractgribdir/get_inv.pl "$url2.inv" > $tmpdir/my_inv2
    res2=$?

    if [ $res2 -ne 0 ] || [ ! -s $tmpdir/my_inv2 ]; then
        echo "no $file2 available yet\n"
        exit
    fi
fi

set -e

#### Extraction from the server and creation of new grib files#####

grep ":$temperature:" < $tmpdir/my_inv1 | $extractgribdir/get_grib.pl "${url1}" $outfile_temperature_1
grep ":$humidity:"    < $tmpdir/my_inv1 | $extractgribdir/get_grib.pl "${url1}" $outfile_humidity_1
grep ":$height:"      < $tmpdir/my_inv1 | $extractgribdir/get_grib.pl "${url1}" $outfile_height_1

grep ":$temperature:" < $tmpdir/my_inv2 | $extractgribdir/get_grib.pl "${url2}" $outfile_temperature_2
grep ":$humidity:"    < $tmpdir/my_inv2 | $extractgribdir/get_grib.pl "${url2}" $outfile_humidity_2
grep ":$height:"      < $tmpdir/my_inv2 | $extractgribdir/get_grib.pl "${url2}" $outfile_height_2


mv $outfile_temperature_1 "$outfile_temperature_1".grib2
mv $outfile_humidity_1    "$outfile_humidity_1".grib2
mv $outfile_height_1      "$outfile_height_1".grib2
mv $outfile_temperature_2 "$outfile_temperature_2".grib2
mv $outfile_humidity_2    "$outfile_humidity_2".grib2
mv $outfile_height_2      "$outfile_height_2".grib2

echo "################################"
echo "#                              #"
echo "# LATITUDE/LONGITUDE SELECTION #"
echo "#                              #"
echo "################################"

index=1
while [ $index -le 2 ]
do
    echo
    echo "PROCESSING OUTPUT FILE" $index

    if [ $index == 1 ]
    then
        outfile_temperature="$tmpdir/temperature_1.grib2"
        outfile_humidity="$tmpdir/humidity_1.grib2"
        outfile_height="$tmpdir/height_1.grib2"
        file_out="C${long_fn}${latitude}D${file1:4:4}"-"${file1:8:2}"-"${file1:10:2}"T"${file1:12:2}".gdas
    fi
    if [ $index == 2 ]
    then
        outfile_temperature="$tmpdir/temperature_2.grib2"
        outfile_humidity="$tmpdir/humidity_2.grib2"
        outfile_height="$tmpdir/height_2.grib2"
        file_out="C${long_fn}${latitude}D${file2:4:4}"-"${file2:8:2}"-"${file2:10:2}"T"${file2:12:2}".gdas
    fi

    ind=1
    P=( 10 20 30 50 70 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 925 950 975 1000 )

    echo "#######################"
    echo "#                     #"
    echo "# WRITING OUTPUT FILE #"
    echo "#                     #"
    echo "#######################"

    while  [ $ind -le 26 ]
    do
        $grib2 $outfile_temperature -d $ind -lola $longitude:1:0.1 $latitude:1:0.1 $tmpdir/tmp2.txt text
        $grib2 $outfile_height      -d $ind -lola $longitude:1:0.1 $latitude:1:0.1 $tmpdir/tmp1.txt text

        if [ $ind -ge 5 ]
        then
            $grib2 $outfile_humidity -d $ind -lola 289.6:1:0.1 -24.6:1:0.1 $tmpdir/tmp3.txt text

            sed '1d' $tmpdir/tmp3.txt > $tmpdir/tmp33.txt
        else
            echo "0" > $tmpdir/tmp33.txt
        fi

        sed '1d' $tmpdir/tmp1.txt > $tmpdir/tmp11.txt
        sed '1d' $tmpdir/tmp2.txt > $tmpdir/tmp22.txt
        ind1=$(($ind-1))

        echo ${P[ind1]} > $tmpdir/tmp00.txt
        paste $tmpdir/tmp00.txt $tmpdir/tmp11.txt $tmpdir/tmp22.txt $tmpdir/tmp33.txt >> $tmpdir/tmp44.txt
        ind=$(($ind+1))

    done

    index=$(($index+1))
#     echo " " >> $tmpdir/tmp44.txt
    echo "# P[hPa] HGT[m] T[K]    RELHUM[%]" >> $tmpdir/tmp44.txt
    $tac $tmpdir/tmp44.txt > $dir$file_out
    rm $tmpdir/tmp44.txt

done
