# This module defines all of the common functions used by VMD processes.


proc floatZero {} {
    return 0.001
}


proc isZero {x} {
    # Is x close enough to zero to be called zero?
    if {[expr {abs($x)}] < [floatZero]} {
        return 1
    } else {
        return 0
    }
}


proc pi {} {
    return 3.14159265359
}


# Convert radians to degrees.
proc rad2deg {rad} { 
    return [expr $rad * 180 / [pi]]
}


# Convert degrees to radians.
proc deg2rad {deg} { 
    return [expr $deg * [pi] / 180]
}


# Calculate the angle (in radians, from 0 to pi) between two vectors (as of now, there is no single function to do this in VMD).
proc vectorAngle {v1 v2} {
    # Use the formula angle = acos((v1 dot v2) / (|v1||v2|)).
    return [expr acos([vecdot $v1 $v2] / ([veclength $v1] * [veclength $v2]))]
}


# The names of the MAPs and Antigen segments.
proc antigenSegment {} {
    return "AGEN"
}


proc mapsSegment {} {
    return "MAPS"
}


# Read an epitope from a file. Should be formatted as:
# Chain1 Residue1
# Chain2 Residue2
# ...
proc readFile {fileName} {
	set f [open $fileName]
	set data [read $f]
	close $f
	return $data
}


# Compute the angle (in radians, from 0 to 2pi) between the x axis and the projection of a vector onto the x-y plane.
proc angleToX {vec} {
    # Obtain the x and y coordinates of the vector.
    set xCoor [lindex $vec 0]
    set yCoor [lindex $vec 1]
    # Calculate the angle between the vector's projection onto the x-y plane and the x axis.
    set rotAngle [expr acos($xCoor / sqrt($xCoor * $xCoor + $yCoor * $yCoor))]
    # Convert to degrees.
    set rotAngle [rad2deg $rotAngle]
    if {$yCoor < 0} {
        # If the y coordinate is negative, subtract the angle from a full rotation.
        set rotAngle [expr {360 - $rotAngle}]
    }
    return $rotAngle
}


# Calculate the transformation matrix needed to minimize the z coordinates of the epitope.
proc transformationZMin {antigen epitope} {
    # Calculate the "arm" vector from the center of geometry of the antigen to that of the epitope.
    set AgCenter [measure center $antigen]
    set arm [vecsub [measure center $epitope] $AgCenter]
    # Create the matrix needed to rotate around the axis formed by the cross product of the arm vector and the desired direction vector (0, 0, -1) by the angle between those two vectors.
    if {[veclength $arm] > [floatZero]} {
        return [trans center $AgCenter offset $AgCenter axis [veccross $arm "0 0 -1"] [vectorAngle $arm "0 0 -1"] rad]
    } else {
        return { {1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1} }
    }
}


# Rotate the antigen's epitope to minimize its z coordinates.
proc rotateZMin {antigen epitope} {
    $antigen move [transformationZMin $antigen $epitope]
}


# Calculate the z angle (in degrees) at which an antigen is facing. This angle is defined arbitrarily as the angle between the positive x axis and the projection onto the x-y plane of the vector from the antigen's center of mass to the C-alpha atom of its first residue, when the antigen has been rotated such that the z coordinates of its epitope are minimized. This definition is in place to provide a standard way to calculate the z rotation of any antigen structure.
proc getZAngle {antigen epitope antigenFirstChain} {
    # The antigen must be molecule 0.
    set antigenMolID 0
    # Transform the coordinates of the C-alpha atom of the first residue with the transformation matrix needed to minimize the z coordinates of the epitope and center the antigen at the origin. Then find the angle its x-y projection makes with the x axis.
    set transformation [transformationZMin $antigen $epitope]
    set transformedAntigen [coordtrans $transformation [measure center $antigen]]
    set transformedCAlpha [coordtrans $transformation [measure center [atomselect $antigenMolID "chain $antigenFirstChain and residue 0 and name CA"]]]
    return [angleToX [vecsub $transformedCAlpha $transformedAntigen]]
}


# Rotate the antigen around the z axis by zRot degrees, ASSUMING that the epitope has had its z coordinates minimized.
proc rotateZAxis {antigen zRot} {
    set AgCenter [measure center $antigen]
    $antigen move [trans center $AgCenter offset $AgCenter axis z $zRot deg]
}


# Calculate the z rotation angle and (x, y, z) position of an antigen's epitope, ASSUMING that the epitope has had its z coordinates minimized.
proc getAntigenPosition {antigen epitope antigenFirstChain} {
    return "[getZAngle $antigen $epitope $antigenFirstChain] [measure center $epitope]"
}


# Perform complete antigen positioning, including minimizing epitope z coordinates.
proc positionAntigen {antigen epitope antigenFirstChain zRot x y z} {
    # Rotate the epitope to minimize its z coordinates.
    rotateZMin $antigen $epitope
    # Rotate the antigen around the z axis so that it points in the direction given by zRot.
    rotateZAxis $antigen [expr $zRot - [getZAngle $antigen $epitope $antigenFirstChain]]
    # Translate the antigen so that its epitope is centered at (x, y, z).
    $antigen moveby [vecsub "$x $y $z" [measure center $epitope]]
}


# Perform partial antigen positioning, assuming that epitope z coordinates have been minimized.
proc repositionAntigen {antigen epitope antigenFirstChain zRot x y z} {
    # Rotate the antigen around the z axis so that it points in the direction given by zRot.
    set dzRot [expr $zRot - [getZAngle $antigen $epitope $antigenFirstChain]]
    if {![isZero $dzRot]} {
        rotateZAxis $antigen $dzRot
    }
    # Translate the antigen so that its epitope is centered at (x, y, z).
    $antigen moveby [vecsub "$x $y $z" [measure center $epitope]]
}


# Move the antigen to the position with the epitope centered at the origin and pointing towards the negative z axis and with the antigen facing zero degrees.
proc mountAntigen {antigen epitope antigenFirstChain} {
    positionAntigen $antigen $epitope $antigenFirstChain 0 0 0 0
}


# Generate a range of numbers (just like range() in Python).
proc range {start end step} {
    set out {}
    set iMax [expr {1 + ((abs($end - $start) - 1) / abs($step))}]
    for {set i 0} {$i < $iMax} {incr i} {
        lappend out [expr {$start + ($i * $step)}]
    }
    return $out
}


# Read an attribute from an experiment.
proc readAttribute {fileName attribute number} {
    set f [open $fileName]
    set getting 0
    set count 0
    # Read each line of the attribute file.
    while {$count < $number && $getting != -1} {
        set getting [gets $f line]
        # If the line has a colon, then it contains an attribute.
        set div [split $line ":"]
        if {[llength $div] > 1} {
            set attr [lindex $div 0]
            # Count the number of times the attribute has been found.
            if {$attr == $attribute} {
                incr count
            }
        }
    }
    # Close the file.
    close $f
    # Make sure the attribute was found.
    if {$count < $number} {
        puts "$attribute not found."
        exit 1
    }
    # If so, return the associated value.
    return [string trim [lindex $div 1]]
}
