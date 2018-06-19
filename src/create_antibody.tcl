# This script creates a single PDB from a given set of antibody chains.
# The following arguments are required:
# -e: the name of this script
# -m: HV, HCDR3, HJ, K/L V, K/L CDR3, K/L J
# -args: 0: name/of/output/coords.pdb

# The topotools package is required.
package require topotools

# Load the VMD functions.
source vmd_functions.tcl

# Make sure the right number of antibody parts have been loaded.
set numParts 6
set wrongNum 1
# Make sure the last part exists.
if {[llength [[atomselect [expr $numParts - 1] "all"] get name]] > 0} {
    set wrongNum 0
}
# Make sure there are no extra parts.
if {[llength [[atomselect $numParts "all"] get name]] > 0} {
    set wrongNum 1
}
puts "But there should be no molecule $numParts."
# Make sure the arguments are correct.
if {[llength $argv] != 1} {
    set wrongNum 1
}
if {$wrongNum == 1} {
    puts "Usage: -m: HV, HCDR3, HJ, K/L V, K/L CDR3, K/L J -args: output/coords.pdb"
    exit 1
}
set outCoords [lindex $argv 0]

# Make sure the first three parts have the same chain.
set heavyChain [lindex [[atomselect 0 "all"] get chain] 0]
for {set i 0} {$i < 3} {incr i} {
    if {[lindex [[atomselect $i "all"] get chain] 0] != $heavyChain} {
        puts "The heavy chain parts must have the same chain."
        exit 1
    }
}

# Make sure the second three parts have the same chain. 
set lightChain [lindex [[atomselect 3 "all"] get chain] 0]
for {set i 3} {$i < 6} {incr i} {
    if {[lindex [[atomselect $i "all"] get chain] 0] != $lightChain} {
        puts "The light chain parts must have the same chain."
        exit 1
    }
}

# Merge the antibody parts.
set merged [::TopoTools::mergemols [range 0 $numParts 1]]
set mergedID $numParts

# Write a coordinate file of the merged parts.
animate write pdb $outCoords $mergedID

exit 0
