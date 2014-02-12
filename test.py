# I'm using this code to just test functions as I make them

import functions as f
import json

# Parametric variables of the helix
bps_turn = 10.5 # Number of base pairs per turn of the double helix
len_bp = .332 # nm, length of the helix per bp
dia = 2.0 # nm, diameter of the duplex
min_distance = 25 # nm, minimum distance separating fluorophore docking sites
max_distance = 6

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

face_centroid_hull = f.centroid_hull(face_pts)
from operator import itemgetter, attrgetter
sorted_face = sorted(face_centroid_hull[0], key=itemgetter(4),reverse=True)
print sorted_face

maximum = max([i[4] for i in face_centroid_hull[0]])
index = [i[4] for i in face_centroid_hull[0]].index(maximum)

print index
print maximum

first_point = f.first_point(face_pts,staples, max_distance, location, len_bp, dia)
print first_point

nearby_strands = f.nearby_threeprime(46,40,staples,max_distance,location,len_bp,dia)

print nearby_strands

#print face_centroid_hull

#print f.first_point(face_pts,suggested_plane,location,dia,len_bp)
