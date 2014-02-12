# AutoFluor 0.1
# This script prepares a DNA nanostructure for super resolution fluorescent microscopy.
# Originally written by Will Patrick - moonshot@mit.edu

# List of things this script cannot do:
# Cannot handle helices that have have a break

import functions as f
import json

# Parametric variables of the helix
len_bp = .332 # nm, length of the helix per bp
dia = 2.0 # nm, diameter of the duplex
min_distance = 25 # nm, minimum distance separating fluorophore docking sites
max_distance = 6
min_num_docking_strands = 5

## USER INTERFACE 

print
print "AutoFluor 0.1"
print
print "Welcome. The purpose of this script is to help prepare a caDNAno file for super resolution fluoresence microscopy."
print
print "The script analyzes a caDNAno JSON file and identifies potential docking sites for fluorescence tagging and nearby 3' staple strands that could be utilized as docking staple strands."
print 
print "In it's current form, this script cannot process all types of JSON files. Before using the script, please read over the brief documentation."
#print "In it's current form the script has several limitations: (1) the script is not currently written to handle stands that have multiple ending and starting points, (2) the viewing planes for the  microscp"

# Get the file from the user
json_name = raw_input("Enter the name of the caDNAno file: ")

# Open and read the file
caDNAno_file = open(json_name, 'r')
file_contents = eval(caDNAno_file.read())
caDNAno_file.close()

strands = file_contents["vstrands"]

print 
print "Got it."
print 
print "Step 1: Select the orientation of the nanostructure."
print 
# Creates the location matrix which has the number, location (row and column) and beginning and ending staple for each strand
location = f.location_matrix(strands)

# Calculates the surface area of the nanostructure projected onto the xy, xz, and yz planes 
areas = [f.area_plane(f.xy_points(location,dia)), f.area_plane(f.xz_points(location,dia,len_bp)), f.area_plane(f.yz_points(location,dia,len_bp))]

if (areas[0] > areas[1]) & (areas[0] > areas[2]):
	suggested_plane = "xy"
elif (areas[1] > areas[0]) & (areas[1] > areas[2]):
	suggested_plane = "xz"
else:
	suggested_plane = "yz"
	
print "Calculations suggest that the " + suggested_plane + " plane is the best orientation to view the nanostructure." 
print 
print "Please select a face of the nanostructure on the " + suggested_plane + " plane to bind to the surface."
print 

point_face = eval(raw_input("Enter a point on the nanostructure face that should be attached to surface. Use the following format [strand number, staple position]: "))

f.streptavidin_staples(suggested_plane,point_face,location,strands,dia,len_bp)

print "Great. That seemed to work. All of the staple strands that should be biotinylated are now colored magenta."
print
print "Now it is time to select the locations for the docking sites"

staples = f.staple_positions(strands,location)

face_pts = f.face_points(staples,location,point_face,suggested_plane,dia,len_bp)

[first_point,hull_face_pts] = f.first_point(face_pts,staples, max_distance, location, len_bp, dia,min_num_docking_strands)

[suggested_sites,enough_sites,possible_combos]=f.remaining_docking_sites(hull_face_pts,first_point,staples,max_distance,location,len_bp,dia,min_num_docking_strands,min_distance)

[a,b,c]=f.find_docking_strands(suggested_sites,strands,staples,max_distance,location,len_bp,dia)

print
print "Docking sites have been found... creating a new JSON file."
print

## Writing JSON file
file_contents_new = file_contents
file_contents_new["vstrands"] = strands
new_name = json_name[:-5] + "_new" + ".json"
with open(new_name, 'w') as outfile:
  json.dump(file_contents_new, outfile)

print "Complete. The new file is " + new_name + "."




