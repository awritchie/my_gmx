mol new rapral.pdb
mol new 3ns_2-E30D_K31E+N29C_180.gro
mol addfile 3ns_2-E30D_K31E+N29C_180.nojump.xtc waitfor all 1 
set outfile [open "output.coords.180" w]
set sel1 [atomselect 0 "backbone and resid 14 to 100 and index > 1334"]
set sel2 [atomselect 1 "backbone and resid 11 to 97"]
set all [atomselect 1 "all"]
set v1 [atomselect 1 "resname CNC and name NE"]
set v2 [atomselect 1 "resname CNC and name CD"
]
set nf [molinfo 1 get numframes]
for {set i 0} {$i < $nf} {incr i} {
     $sel2 frame $i 
     $all frame $i 
     set M [measure fit $sel2 $sel1 weight mass] 
     $all move $M 
     $v1 frame $i
     $v2 frame $i
     set V1 [lindex [ $v1 get {x y z}] 0]
     set V2 [lindex [ $v2 get {x y z}] 0]
     set simdata($i.r) [vecsub $V1 $V2] 
     puts $outfile "$simdata($i.r)"
     puts $simdata($i.r)
     }
puts "Done writing coordinates, please exit"
close $outfile
quit

