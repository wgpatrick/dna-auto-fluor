# I'm using this code to just test functions as I make them

import functions as f
import json

# Parametric variables of the helix
bps_turn = 10.5 # Number of base pairs per turn of the double helix
len_bp = .332 # nm, length of the helix per bp
dia = 2.0 # nm, diameter of the duplex
min_distance = 25 # nm, minimum distance separating fluorophore docking sites
max_distance = 6
min_num_docking_strands = 5


# Get the file from the user
json_name = "hexagon_prism_final_barcode_7.json"

# Open and read the file
caDNAno_file = open(json_name, 'r')
file_contents = eval(caDNAno_file.read())
caDNAno_file.close()

strands = file_contents["vstrands"]

location = f.location_matrix(strands)

staples = f.staple_positions(strands,location)

areas = [f.area_plane(f.xy_points(location,dia)), f.area_plane(f.xz_points(location,dia,len_bp)), f.area_plane(f.yz_points(location,dia,len_bp))]

if (areas[0] > areas[1]) & (areas[0] > areas[2]):
	suggested_plane = "xy"
elif (areas[1] > areas[0]) & (areas[1] > areas[2]):
	suggested_plane = "xz"
else:
	suggested_plane = "yz"

loc = [0,23]

face_pts = f.face_points(staples,location,loc,suggested_plane,dia,len_bp)

[first_point,hull_face_pts] = f.first_point(face_pts,staples, max_distance, location, len_bp, dia,min_num_docking_strands)


[suggested_sites,enough_sites,possible_combos]=f.remaining_docking_sites(hull_face_pts,first_point,staples,max_distance,location,len_bp,dia,min_num_docking_strands,min_distance)

[a,b,c]=f.find_docking_strands(suggested_sites,strands,staples,max_distance,location,len_bp,dia)

print a

#print face_centroid_hull

#print f.first_point(face_pts,suggested_plane,location,dia,len_bp)
