# Merge an antigen and a MAPs part into a single molecule and generate a PDB and PSF.
# The following command-line arguments are required:
# -e: the name of this script
# -args: 0: the coordinates of the antigen
#        1: the coordinates of the MAPs part
#        2: the prefix of the combined structure
#        3: the first topology file
#        4: (optional) the second topology file
#        etc.

# Load the VMD functions.
source vmd_functions.tcl

# The psfgen package is required.
package require psfgen

# Load the arguments.
if {[llength $argv] < 4} {
	puts "Usage: -args antigen/coordinates.pdb MAPs/part/coordinates.pdb combined/structure/prefix topology/file1.rtf (optional) topology/file2.rtf ..."
	exit 1
}
set AgPDB [lindex $argv 0]
set MAPsPDB [lindex $argv 1]
set combined [lindex $argv 2]

# Read the topology file and alias the atoms.
for {set i 3} {$i < [llength $argv]} {incr i} {
	topology [lindex $argv $i]
}
pdbalias residue HIS HSE
pdbalias atom ILE CD1 CD

# Build a segment for the antigen.
segment [antigenSegment] {
	pdb $AgPDB
}

# Read antigen coordinates into the antigen segment.
coordpdb $AgPDB [antigenSegment]

# Build a segment for the MAPs part. Do not patch the terminal residues.
segment [mapsSegment] {
    first none
    last none
	pdb $MAPsPDB
}

# Read the immunoglobulin coordinates into the MAPs segment.
coordpdb $MAPsPDB [mapsSegment]

# Create the PDB and PSF.
writepdb "$combined.pdb"
writepsf "$combined.psf"

exit 0
