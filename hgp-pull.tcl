#================================================================================
# HARMONIC GUIDING POTENTIAL (HGP) along an arbitrary RC trajectory 
#  - trajectory defined by a given number of points assumed to be 
#    connected by straight segments
#  - angle restrain along z-axis
#  INPUT: - nopulls     -> number of pullings (j = 1..nopulls)
#         - N0          -> # of subsegments along a trajectory segment
#         - selFile$j   -> file with atomIDs of SMD atoms for pull $j
#         - HGPpath$j   -> file with trajectory points for pull $j
#         - HGPk        -> elastic const. of HGP for pulling
#         - restrFile$j -> file with atomdIDs of angle restriction
#                          (around z axis) for pull $j
#         - HGPkq       -> elastic const. of HGP for angle restriction
#================================================================================

# Define some useful functions
#-----------------------------
# Function:  veclength {v}
#  Returns:    the vector length
proc veclength {v} {
    set retval 0
    foreach term $v {
	set retval [expr $retval + $term * $term]
    }
    return [expr sqrt($retval)]
}

# Function:  vecnorm {v}
#  Returns:    the normal vector pointing along v
proc vecnorm {v} {
    set sum 0
    foreach term $v {
	set sum [expr $sum + $term * $term]
    }
    set sum [expr sqrt($sum)]
    set retval {}
    foreach term $v {
	lappend retval [expr $term / $sum]
    }
    return $retval
}

# Function:  vecdot {v1 v2}
#  Returns:    the dot product of vectors v1 and v2
proc vecdot {v1 v2} {
    set sum 0
    foreach term1 $v1 term2 $v2 {
	set sum [expr $sum + $term1 * $term2]
    }
    return $sum
}


#---------------------------------------------------------#
# NOTE: VMD (NAMD-tclForces) uses base 0 (1) atom IDs !!! #
#---------------------------------------------------------#

# Initialize counters
set step   $firststep                  ;# initial timestep
set Nseg   [expr int($step / $N0)]     ;# initial segment;  $Nseg <= $NS !
set Tseg   [expr $step - $Nseg * $N0]  ;# initial position on segment $Nseg; $Tseg <= $N0

#---------------------------------------------------------#

for {set i 1} {$i <= $nopulls} {incr i} {

    # Read in atom IDs of the selection
    set selfile [subst $[subst selFile$i]]
    set id [open $selfile r]
    gets $id sel
    close $id
    print "HGP: \#No. of SMD atoms in $selfile = [llength $sel]"
    lappend sels $sel

    # Read in trajectory points
    set hgppathfile [subst $[subst HGPpath$i]]
    set id [open $hgppathfile r]
    set path {}
    while {[gets $id line] >= 0} { lappend path $line }
    close $id
    print "HGP: \#No. of trajectory points in $hgppathfile = [llength $path]"
    lappend paths $path

    # Calculate stepping vector dR0 along segments
    # dR0_i = (R0_{i+1} - R0_i) / N0
    set NS [expr [llength $path] - 1]
    if {$Nseg > $NS} {
	print "ERROR: Pull \#$i: Initial timestep ($step) larger than the total timesteps required for simulation ([expr $NS * $N0])!"
	print "ERROR: Pull \#$i: Initial timestep must be 0 unless continuing an interrupted pulling."
	exit
    }
    set j 0; set dR0 {}
    while {$j < $NS} {
	set j1 [expr $j + 1]
	lappend dR0 [vecscale [expr 1. / $N0] [vecsub [lindex $path $j1] [lindex $path $j]]]
	set j $j1
    }
    lappend NSs $NS; lappend dR0s $dR0

    # Create namd2 groups
    lappend gs [addgroup $sel]

    # Initialize number of segments for each path
    lappend Maxsegs [expr [llength $path] - 1] ;# no. of segments is (no. of points) - 1

    # Initialize target RC:  R0(t0) = [lindex $path 0]
    set R0 [vecadd [lindex $path $Nseg] [vecscale $Tseg [lindex $dR0 $Nseg]]]
    print "HGP: N0 = $N0 ; step = $step ; Nseg = $Nseg ; Tseg = $Tseg ; R0$i = $R0"
    lappend R0s $R0

    # Initialize external work
    lappend works 0.

    # Print headers
    print "HGP${i}1: \# TS   R$i           R0$i          F$i"
    print "HGP${i}2: \# TS   Work   Nseg   Tseg   (for path $i)"
}

print "HGP: outFreq = $outFreq"

#======================================================================

proc calcforces {} {
    global nopulls gs HGPk N0 step Nseg Tseg Maxsegs dR0s R0s outFreq works

    # load group COM coordinates R_i(t) , i = 0,...,N0-1
    loadcoords coor

    # calculate HGP force on selected atoms
    set Fs {}
    for {set i 0; set i1 1} {$i < $nopulls} {incr i; incr i1} {

	# pulling force
	set R      $coor([lindex $gs $i])
	set dR     [vecsub [lindex $R0s $i] $R]
	set n0     [vecnorm [lindex [lindex $dR0s $i] $Nseg]]
	set dRn    [vecdot $dR $n0]
	set F      [vecscale [expr $HGPk * $dRn] $n0]
	lappend Fs $F
	addforce   [lindex $gs $i] $F

	# output data 1
        if {$step % $outFreq == 0} {
	    print  "HGP${i1}1: $step $R [lindex $R0s $i] $F"
	}
    }

    # calculate new target RC: R0
    if {$Tseg == $N0} {
	set Tseg 0
	incr Nseg
    }
    for {set i 0; set i1 1} {$i < $nopulls} {incr i; incr i1} {
	if {$Nseg < [lindex $Maxsegs $i]} {
	    set dR0T [lindex [lindex $dR0s $i] $Nseg]
	    set R0s  [lreplace $R0s $i $i [vecadd [lindex $R0s $i] $dR0T]]

	    # update the external work
	    set work [expr [lindex $works $i] + [vecdot [lindex $Fs $i] $dR0T]]
	    set works [lreplace $works $i $i $work]
	    # output data 2
	    if {$step % $outFreq == 0} {
		print "HGP${i1}2: $step $work $Nseg $Tseg"
	    }
	} else {
	    print "HGP: pull ${i1}: Nseg: $Nseg; Maxseg: [lindex $Maxsegs $i]"
	    print "HGP: pull ${i1}: End of trajectory! Please ignore following HGP data (if any)."
	}
    }

    # increment time step
    incr Tseg
    incr step
}
