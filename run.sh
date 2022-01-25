#!/bin/sh
echo ""
echo "################################################################"
echo "##                                                            ##"
echo "##      This is to run an assimilation system that performs   ##"
echo "##      3dvar cycle using the Lorenz 40-var model.            ##"
echo "##                                                            ##"
echo "##      Author: Chanh Q. Kieu, Research Associate             ##"
echo "##      Dept. of Atmospheric and Oceanic Science              ##"
echo "##      Univ. of Maryland, College Park                       ##"
echo "##      email: kieucq@atmos.umd.edu                           ##"
echo "##                                                            ##"
echo "################################################################"
echo ""
#echo "Enter the number of cycles that you want to perform the exp"
#read ncycle
ncycle=21
npair=0
nensemble=$((2*$npair+1))
ninterval=50
cold_start="F"
echo "SUMMARY OF ETKF CONFIGURATION"
echo "==========================================================="
echo " Number of cycles that the script will run is    : $ncycle"
echo " Number of the ensemble members for this test is : $nensemble"
echo " Interval for doing the analysis is              : $ninterval"
echo " Cold start option is                            : $cold_start"
echo "==========================================================="
#
# start looping over the entire cycles now
#
rm -f *.dat
icycle=1
filetime=0
echo ""
echo "Perfoming a cycle run now ..."
while [ "$icycle" -le "$ncycle" ]
do
 echo ""
 echo " Start doing the cycle $icycle"th ...
#
# create an indexing for the file name that indicate
# the cycle we are running. This step is needed only
# for $icycle -eq 1. After that the ofile will be
# updated at the end of this loop
#
 echo " indexing filetime ..."
 if [ "$filetime" -lt 10 ]; then
  ofile="000$filetime.dat"
 elif [ "$filetime" -lt 100 ]; then
  ofile="00$filetime.dat"
 elif [ "$filetime" -lt 1000 ]; then
  ofile="0$filetime.dat"
 else
  ofile="$filetime.dat"
 fi
 echo " file extention is $ofile"
 sleep 1 
#
# checking the background data first
#
 echo " checking background data ..."
 bgdfile="bgd$ofile"
 if [ -f "./ana/$bgdfile" ]; then
  echo " Background file $bgdfile exists...linking now"
  ln -sf ./ana/$bgdfile ./bgd.dat
 else
  echo " background file $bgdfile does not exist...exit 2"
  exit 2
 fi
#
# checking the observation data now
#
 echo " linking observational file now ..."
 obsfile="obs$ofile"
 if [ -f "./obs/$obsfile" ]; then
  echo " obs $obsfile file exists ...linking now"
  ln -sf ./obs/$obsfile ./obs.dat
 else
  echo " obs file $obsfile does not exist. exit 3"
  exit 3
 fi
#
# perform a 3dvar assimilation now and remove the links
#
 echo " 3dvar is called now ..."
 ./3dvar.exe
 rm obs.dat bgd.dat
#
# running the model now
#
 echo " running the model now..."
 ./L40.exe
#
# back up the forecast from the model runs
#
 echo " backing up the forecast and analysis output from L40 model"
 fscfile="fsc$ofile"
 anafile="ana$ofile"
 echo " forecast file that will be stored is $fscfile ..."
 mv fsc.dat ./fsc/$fscfile
 echo " analysis file that will be stored is $anafile ..."
 mv ana.dat ./ana/$anafile
#
# update the next cycle now
#
 icycle=$(($icycle+1))
 filetime=$(($filetime+$ninterval))
#
# update the background at the next cycle now
#
 if [ "$filetime" -lt 10 ]; then
  ofile="000$filetime.dat"
 elif [ "$filetime" -lt 100 ]; then
  ofile="00$filetime.dat"
 elif [ "$filetime" -lt 1000 ]; then
  ofile="0$filetime.dat"
 else
  ofile="$filetime.dat"
 fi
 echo " updated file extention is $ofile"
 bgdfile="bgd$ofile"
 echo " background for the next analysis is $bgdfile"
 mv bgd.dat ./ana/$bgdfile
done
#
# cleaning some mess now
#
echo "Cleaning some unused data"
echo ""
rm -f 3dvar.dat debug.txt log
echo "DONE THE ENSEMBLE RUN SUCCESSFULLY"
exit 0
