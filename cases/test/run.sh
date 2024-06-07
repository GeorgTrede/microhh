# Run the input python file
python3 ${1}_input.py

# find number of processes needed
# -> stored in .ini file under [master] -> npx, npy
#    but make sure that commented lines are not included, i.e. lines starting with #
npx=$(awk -F "=" '/npx=/ && !/#/ {print $2}' ${1}.ini)
npy=$(awk -F "=" '/npy=/ && !/#/ {print $2}' ${1}.ini)
proc=$(($npx * $npy))

# Determine which version of microhh to use based on the hostname
if [ "$HOSTNAME" == "zilxap03" ]; then
    command_prefix="microhh_hera"
else
    command_prefix="microhh"
fi

# Create a temporary file for storing command output
temp_file=$(mktemp)

# Define a function to run a command and display its most recent output
run_and_display() {
    local command=$1
    # Execute the command and redirect output to the temp file
    $command > "$temp_file" &
    local cmd_pid=$!

    # Loop to update the display with the most recent output
    while kill -0 $cmd_pid 2>/dev/null; do
        # Capture the most recent line of output without introducing a new line
        last_line=$(tail -n 1 "$temp_file")
        
        # Clear the line and print the last line without a newline at the end
        echo -ne "\r\033[K${last_line}"
        
        # Wait for 0.25 seconds before checking the output again
        sleep 0.25
    done

    # Print a newline to separate the output from the next command
    echo
}

# Run the commands with the determined prefix and display their most recent output
run_and_display "mpiexec -n $proc ${command_prefix} init ${1}"
run_and_display "mpiexec -n $proc ${command_prefix} run ${1}"

# Clean up
rm "$temp_file"

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

# find sbot[thl] of the simulation
# -> stored in .ini file under [boundary] -> sbot[thl]
#    but make sure that commented lines are not included, i.e. lines starting with #
#    also, this only exists optionally
#    if it does not exist, set to "X"
if [ $(awk -F "=" '/sbot\[thl\]=/ && !/#/ {print $2}' ${1}.ini) ]; then
    sbot=$(awk -F "=" '/sbot\[thl\]=/ && !/#/ {print $2}' ${1}.ini)
else
    sbot="X"
fi
sbot="thl${sbot}"

# create folders snaps_res_uflux_sbot/nc_files and snaps_res_uflux_sbot/npy_files
mkdir -p snaps_${res}_${uflux}_${sbot}/nc_files
mkdir -p snaps_${res}_${uflux}_${sbot}/npy_files

# run cross_to_nc.py and suppress the output violently
python3 cross_to_nc.py -n ${proc} > /dev/null 2>&1

# move the created .nc files to snaps_res_uflux/nc_files
mv *.nc snaps_${res}_${uflux}_${sbot}/nc_files

# clear frames folder
rm frames/*

# run eval.py
# python3 eval.py $@

# move all files *.00* and *.01* to snaps_res_uflux/
mv *.00* snaps_${res}_${uflux}_${sbot}
mv *.01* snaps_${res}_${uflux}_${sbot}
mv *.02* snaps_${res}_${uflux}_${sbot}

# remove .out files
rm *.out

# echo that the run is finished and the files are stored in the respective folders
echo "Run finished. Files are stored in snaps_${res}_${uflux}_${sbot}/"