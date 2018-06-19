# This VMD script is the NEWEST one (2.1.0) that  calculates the interaction energy between the antigen and the MAPs part for each antigen position.
# The following arguments are required:
# -e: the name of this script
# -args: 0: the antigen+part PSF
#        1: the antigen+part PDB
#        2: the positions file
#        3: the epitope file
#        4: antigenFirstChain
#        5: the name of the output file
#        6: the first parameter file
#        7: the path to the NAMD executable
#        8: (optional) the second parameter file
#        etc.

# Load the VMD functions.
source vmd_functions.tcl
package require namdenergy

# Load the arguments.
if {[llength $argv] < 8} {
	puts "Usage: structure.psf coordinates.pdb positions/file.dat epitope/file.txt antigenFirstChain output/energies.dat experiment/details.txt path/to/namd parameter/file1.prm parameter/file2.prm (optional) ..."
	exit 1
}
set inStruct [lindex $argv 0]
set inCoords [lindex $argv 1]
set positionsFile [lindex $argv 2]
set epitopeFile [lindex $argv 3]
set antigenFirstChain [lindex $argv 4]
set outEnergy [lindex $argv 5]
set NamdCommand [lindex $argv 6]
set parameterFiles []
for {set i 7} {$i < [llength $argv]} {incr i} {
	lappend parameterFiles [lindex $argv $i]
}

# Load the antigen and the MAPs part.
mol new $inCoords
mol addfile $inStruct

set agSeg [antigenSegment]
set mapsSeg [mapsSegment]

# Select every atom in each segment. These segments were named by merge_antigen_part.tcl.
set Ag [atomselect 0 "segname $agSeg"]
set MAPsPart [atomselect 0 "segname $mapsSeg"]

# Select the epitope.
set epitope [atomselect 0 [readFile $epitopeFile]]

# Make a file of interaction energies.
set eFile [open $outEnergy "w"]
close $eFile

# Calculate the energy at all positiions.
set positions [open $positionsFile]
while {[gets $positions p] >= 0} {
    # Move the antigen to the correct position.
    set p [split $p]
    set zAngle [lindex $p 0]
    set x [lindex $p 1]
    set y [lindex $p 2]
    set z [lindex $p 3]
    puts "Position: $zAngle $x $y $z"
    repositionAntigen $Ag $epitope $antigenFirstChain $zAngle $x $y $z
    # Calculate the interaction energy (electrostatic and van der Waals) between the antigen and the MAPs part.
    set energyList [lindex [namdenergy -sel $Ag $MAPsPart -elec -vdw -par $parameterFiles -exe $NamdCommand] 0]
    # Read the energy from the fourth position in the list.
    set energy [lindex [split $energyList] 4]
    puts "Energy: $energy"
    # Write the energy to the file of all energies.
    set eFile [open $outEnergy "a"]
    puts $eFile "$zAngle $x $y $z $energy"
    close $eFile
}
close $positions

exit 0
