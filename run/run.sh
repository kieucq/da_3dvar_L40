#!/bin/sh
#
# NOTE: This is to run an assimilation system that performs 3dvar cycle 
#       using the Lorenz 40-var model.
#
# HIST: - Aug 01, 2009: created by CKi during postdoc period at UMD
#       - Apr 09, 2022: modified for G594 DA class by CK
#
# AUTH:  Chanh Kieu (ckieu@iu.edu)            
#
#=======================================================================
ncycle=21
npair=0
nensemble=$((2*$npair+1))
ninterval=50
cold_start="T"
homedir="/N/u/ckieu/Karst/model/da_3dvar_L40"
echo ""
echo "              SUMMARY OF THE 3DVAR CONFIGURATION           "
echo "==========================================================="
echo " Number of cycles that the script will run is    : $ncycle"
echo " Number of the ensemble members for this test is : $nensemble"
echo " Interval for doing the analysis is              : $ninterval"
echo " Cold start option is                            : $cold_start"
echo " Home dir for running this 3DVAR is: ${homedir}"
echo "==========================================================="
#
# start looping over the entire cycles now
#
rm -f *.dat
icycle=1
filetime=0
if [ "$cold_start" = "T"  ]; then
    echo "Initialize a cold start with a OSSE mode now"
    #
    # running the truth integration for error calculation
    sed 's/RESTART/'$ninterval'/' namelist.template > namelist.L40 
    cd ${homedir}/truth
    time ./truth.exe >& ${homedir}/run/log.txt
    #
    # creating OSSE obs from the truth output
    cd ${homedir}/obs
    ln -sf ${homedir}/truth/truth*.dat ./
    time ./obs.exe >> ${homedir}/run/log.txt
    #
    # running a ctl from a pertubed truth for reference
    cd ${homedir}/ctl
    ln -sf ${homedir}/run/namelist.L40 ./
    time ./ctl.exe >> ${homedir}/run/log.txt
fi
echo ""
echo "Perfoming a 3DVAR DA cycle now..."
cd ${homedir}/run/
while [ "$icycle" -le "$ncycle" ]
do
    echo "Start doing the cycle $icycle"th ...
    #
    # create an indexing for the file name that indicate
    # the cycle we are running. This step is needed at  
    # begining of each cycle. The ofile will be
    # updated again at the end of this loop for next bgd
    #
    echo " First, indexing the file timestamp..."
    if [ "$filetime" -lt 10 ]; then
        ofile="000$filetime.dat"
    elif [ "$filetime" -lt 100 ]; then
        ofile="00$filetime.dat"
    elif [ "$filetime" -lt 1000 ]; then
        ofile="0$filetime.dat"
    else
        ofile="$filetime.dat"
    fi
    echo " File extention is $ofile"
    #
    # checking the background data first
    #
    echo " checking background data ..."
    bgdfile="bgd$ofile"
    if [ -f "${homedir}/ana/$bgdfile" ]; then
        echo " Background file $bgdfile exists...linking now"
        ln -sf ${homedir}/ana/$bgdfile ./bgd.dat
    else
        echo " background file $bgdfile does not exist...exit 2"
        exit 2
    fi
    #
    # checking the observation data now
    #
    echo " linking observational file now ..."
    obsfile="obs$ofile"
    if [ -f "${homedir}/obs/$obsfile" ]; then
        echo " obs $obsfile file exists ...linking now"
        ln -sf ${homedir}/obs/$obsfile ./obs.dat
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
     mv fsc.dat ${homedir}/fsc/$fscfile
     echo " analysis file that will be stored is $anafile ..."
     mv ana.dat ${homedir}/ana/$anafile
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
     mv bgd.dat ${homedir}/ana/$bgdfile
done
#
# cleaning some mess now
#
echo "Cleaning some unused data"
echo ""
rm -f 3dvar.dat debug.txt log
echo "DONE THE ENSEMBLE RUN SUCCESSFULLY"
exit 0
