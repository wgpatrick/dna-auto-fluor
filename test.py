# Some code to use while testing

import functions as f
import json

# Parametric variables of the helix
bps_turn = 10.5 # Number of base pairs per turn of the double helix
len_bp = .332 # nm, length of the helix per bp
dia = 2.0 # nm, diameter of the duplex
min_distance = 25 # nm, minimum distance separating fluorophore docking sites
max_distance = 6
min_num_docking_strands = 5
docking_site_color_1 = 16528437 # red
docking_site_color_2 = 1870874 # green
docking_site_color_3 = 1250225 # blue
three_p_length = 7 # length of three prime docking strand extension

# Get the file from the user
json_name = "hexagon_prism_pre_docking.json"

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


[a,b,c]=f.find_docking_strands(suggested_sites,strands,staples,max_distance,location,len_bp,dia,docking_site_color_1,docking_site_color_2,docking_site_color_3)



new_strands = f.ghost_strands(a,b,c,location,strands,three_p_length)
test_strand = new_strands[-3]





file_contents_new = file_contents
file_contents_new["vstrands"] = strands
new_name = json_name[:-5] + "_new" + ".json"
with open(new_name, 'w') as outfile:
  json.dump(file_contents_new, outfile)

print "Complete. The new file is " + new_name + "."

#print face_centroid_hull

#print f.first_point(face_pts,suggested_plane,location,dia,len_bp)
