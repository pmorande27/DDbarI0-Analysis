#!/bin/bash
echo "Getting all amplitudes"
output=$(python3 /home/pm843/bin/list_all_amplitude_plots.py)
# Convert the output into an array
arr=($output)

# Iterate over the array and echo each element
for item in "${arr[@]}"; do

    plot_amplitude "./plots/$item"
    new_filename="${$item%.xml}.dat"
    new_filename2="${new_filename/Amplitude/}"
    mv "plot.dat" "./plots/$new_filename2"
done
