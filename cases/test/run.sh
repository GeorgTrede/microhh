# Run the input python file
python3 ${1}_input.py

# find number of processes needed
# -> stored in .ini file under [master] -> npx, npy
#    but make sure that commented lines are not included, i.e. lines starting with #
npx=$(awk -F "=" '/npx=/ && !/#/ {print $2}' ${1}.ini)
npy=$(awk -F "=" '/npy=/ && !/#/ {print $2}' ${1}.ini)
proc=$(($npx * $npy))

# check which microhh version to use:
# if the hostname is zilxap03, then use microhh_hera, else use microhh
if [ "$HOSTNAME" == "zilxap03" ]; then
    # run microhh_hera
    mpiexec -n $proc microhh_hera init ${1}
    mpiexec -n $proc microhh_hera run ${1}
else
    # run microhh
    mpiexec -n $proc microhh init ${1}
    mpiexec -n $proc microhh run ${1}
fi

# find resolution of the simulation
# -> stored in .ini file under [grid] -> itot, jtot, ktot
#    but make sure that commented lines are not included, i.e. lines starting with #
itot=$(awk -F "=" '/itot=/ && !/#/ {print $2}' ${1}.ini)
jtot=$(awk -F "=" '/jtot=/ && !/#/ {print $2}' ${1}.ini)
ktot=$(awk -F "=" '/ktot=/ && !/#/ {print $2}' ${1}.ini)
res="${itot}_${jtot}_${ktot}"

# find uflux of the simulation
# -> stored in .ini file under [force] -> uflux
#    but make sure that commented lines are not included, i.e. lines starting with #
#    also, this only exists if [force] -> swlspres=uflux,
#    otherwise set to "X"

# check whether string of swlspres is identical to "uflux"
if [ $(awk -F "=" '/swlspres=/ && !/#/ {print $2}' ${1}.ini) == "uflux" ]; then
    uflux=$(awk -F "=" '/uflux=/ && !/#/ {print $2}' ${1}.ini)
    # round uflux to nearest integer
    uflux=$(echo $uflux | awk '{print int($1+0.5)}')
else
    uflux="X"
fi
uflux="uflux${uflux}"

# create folders snaps_res_uflux/nc_files and snaps_res_uflux/npy_files
mkdir -p snaps_${res}_${uflux}/nc_files
mkdir -p snaps_${res}_${uflux}/npy_files

# run cross_to_nc.py
python3 cross_to_nc.py -n ${proc}

# move the created .nc files to snaps_res_uflux/nc_files
mv *.nc snaps_${res}_${uflux}/nc_files

# run eval.py and suppress the output violently
python3 eval.py $@ > /dev/null 2>&1

# move all files *.00* and *.01* to snaps_res_uflux/
mv *.00* snaps_${res}_${uflux}
mv *.01* snaps_${res}_${uflux}

# remove .out files
rm *.out