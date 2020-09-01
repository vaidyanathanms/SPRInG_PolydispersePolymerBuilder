package require psfgen 
topology lignin.top 	 
set name   poplar 

resetpsf
segment  li {
 residue 1 SYR
 residue 2 GUAI 
 residue 3 GUAI
 residue 4 SYR
 residue 5 SYR 
 residue 6 SYR
 residue 7 SYR
 residue 8 GUAI 
 residue 9 SYR
 residue 10 SYR
 residue 11 SYR
 residue 12 GUAI
 residue 13 SYR
 residue 14 GUAI
 residue 15 SYR
 residue 16 SYR 
 residue 17 GUAI
 residue 18 SYR
 residue 19 GUAI
 residue 20 GUAI
 residue 21 SYR
 residue 22 SYR
 residue 23 SYR
 residue 24 GUAI
 residue 25 GUAI
 residue 26 SYR
 residue 27 SYR
 residue 28 SYR
 residue 29 SYR
 residue 30 SYR
 residue 31 GUAI 
}

patch 4O5  li:1 li:2
patch BO4R li:2 li:3
patch  1BL li:3 li:4
patch O4BR li:4 li:5
patch O4BL li:5 li:6
patch O4BR li:6 li:7
patch 4O5  li:7 li:8 
patch BO4L li:8 li:9
patch BO4R li:9 li:10
patch BO4L li:10 li:11
patch BO4R li:11 li:12
patch  BB  li:12 li:13
patch O4BL li:13 li:14
patch O4BR li:14 li:15
patch O4BL li:15 li:16
patch 4O5  li:16 li:17 
patch BO4R li:17 li:18
patch BO4L li:18 li:19
patch BO4R li:19 li:20
patch  BB  li:20 li:21
patch O4BL li:21 li:22
patch O4BR li:22 li:23
patch 4O5  li:23 li:24 
patch BO4L li:24 li:25
patch BO4R li:25 li:26
patch BO4L li:26 li:27
patch BO4R li:27 li:28
patch  1BL li:28 li:29
patch O4BR li:29 li:30
patch 4O5  li:30 li:31

regenerate angles dihedrals
coordpdb poplar.pdb 
guesscoord 	 
writepdb $name.pdb	 
writepsf $name.psf


exit
