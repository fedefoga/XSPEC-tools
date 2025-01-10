# Improved version of the script
proc fetch_data {type {index 1}} {
    # Helper function to fetch and process XSPEC output data
    set raw [tcloutr plot ldata $type $index]
    return [lsearch -all -inline -not -exact [split $raw " "] {}]
}

proc validate_model {model} {
    # Validate the model name to avoid issues with file handling
    if {[regexp {[^a-zA-Z0-9_]} $model]} {
        puts stderr "Error: Invalid model name '$model'. \
        Only alphanumeric and underscore characters are allowed."
        exit 1
    }
}

proc extract_residuals {model {energy_min 0.3} {energy_max 12.0}} {
    # Validate input parameters
    validate_model $model

    if {![file exists "${model}.xcm"]} {
        puts stderr "Error: Base model file ${model}.xcm not found!"
        exit 1
    }

    # Configure XSPEC environment
    chatter 0
    query yes

    # Load base data and configure energy ranges
    @data.xcm
    setp en
    notice all
    ignore bad
    ignore *:**-$energy_min $energy_max-**
    energies 0.01 1000 1000 log

    # Perform the fit for the base model
    @${model}
    fit

    set xd [fetch_data "x"]
    set xe [fetch_data "xerr"]
    set yd0 [fetch_data "y" 1]
    set ye0 [fetch_data "yerr" 1]
    set mo0 [fetch_data "model" 1]

    set nchan [llength $xd]

    # Write the residuals to the output file
    set output_file "residuals_pn_${model}.txt"
    set fp [open $output_file w]

    for {set k 0} {$k<$nchan} {incr k} {
    	puts $fp [format "%15.8e\t%15.8e\t%15.8e\t%15.8e\t%15.8e" \
        		[lindex $xd $k] [lindex $xe $k] \
        		[lindex $yd0 $k] [lindex $ye0 $k] [lindex $mo0 $k] \
       	 		]
    }
    close $fp

    puts "Residuals saved to $output_file"
}

# Example usage
# extract_residuals "my_model" "./output" 0.3 12.0

