# This VMD script moves an antigen to a specific position specified by a z rotation angle and an x, y, z epitope center.
# The following arguments are required:
# -e: the name of this script.
# -m: 0: antigen/coords.pdb
# -args: 0: epitope/file.txt
#        1: antigenFirstChain
#        2: output/coords.pdb
#        3: zAngle
#        4: x
#        5: y
#        6: z

# Load the VMD functions.
source vmd_functions.tcl

# Load the arguments.
if {[llength $argv] != 7} {
    puts "Usage: -args epitope/file.txt antigenFirstChain output/file.dat zAngle x y z"
    exit 1
}

set epitopeFile [lindex $argv 0]
set antigenFirstChain [lindex $argv 1]
set outputFile [lindex $argv 2]
set zAngle [lindex $argv 3]
set x [lindex $argv 4]
set y [lindex $argv 5]
set z [lindex $argv 6]

# Select every atom in each molecule.
set Ag [atomselect 0 "all"]

# Make sure the antigen exists.
if {[llength [$Ag get name]] == 0} {
    puts "One or more molecules failed to load."
    exit 1
}

# Select the epitope.
set epitope [atomselect 0 [readFile $epitopeFile]]

# Move the antigen to the correct position.
positionAntigen $Ag $epitope $antigenFirstChain $zAngle $x $y $z

# Write a file for the antigen.
$Ag writepdb $outputFile

exit 0
