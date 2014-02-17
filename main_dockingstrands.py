# AutoFluor 0.1
# This script prepares a caDNAno design of a DNA nanostructure for super resolution fluorescent microscopy.
# Originally written by Will Patrick - moonshot@mit.edu - MIT Media Lab
# Feb 2014

####################
### CAPABILITIES ###
####################

# 1. Can take any given JSON file
# 2. The user can vary several parameters: (1) The minimum distance between docking sites, (2) 
# the minimum number of docking strands per docking sites, and  (3) the number of strands per 
# docking site
# 3. Optimal docking sites will automatically be selected on the structure.
# 4. Docking strands nearby docking sites are selected and their color is changed in the JSON file.
# 5. Strands requiring biotintlyation are found. Their color is change to magenta.

###################
### LIMITATIONS ###
###################

# 1. Cannot handle non-continous strands.
# 2. Can only site fluorophores at 3 docking sites 
# 3. Can only attach the nanostructure using the three cartesian planes: xy, xz, yz
# 4. The code cannot select the optimal attachment face onto the streptavidin surface on its own. It asks the user
# 5. Can only prepare files produced using caDNAno 1.0

###################
### LIBRARIES #####
###################

import functions as f # There are additional libraries in functions.py
import json

#####################
##GLOBAL VARIABLES ##
#####################

# Parametric variables of the helix
len_bp = .332 # nm, length of the helix per bp
dia = 2.0 # nm, diameter of the duplex
all_staples_color = 4013375 #dark grey
docking_site_color_1 = 16528437 # red
docking_site_color_2 = 1870874 # green
docking_site_color_3 = 1250225 # blue
three_p_length = 7 # length of three prime docking strand extension

###########################
## USER INTERFACE BEGINS ## 
###########################

print
print "AutoFluor 0.1"
print
print "Welcome. The purpose of this script is to help prepare a caDNAno file for super resolution fluoresence microscopy."
print
print "The script analyzes a caDNAno JSON file and identifies potential docking sites for fluorescence tagging and nearby 3' staple strands that could be utilized as docking staple strands."
print 
print "In it's current form, this script cannot process all types of JSON files. Before using the script, please read over the brief documentation."

#################
### INPUTS ######
#################

json_name = raw_input("Name of the caDNAno file: ")
min_distance = eval(raw_input("Minimum distance between fluorophore docking sites [suggestion: 25nm]: "))# suggestion: 25 nm, minimum distance separating fluorophore docking sites
max_distance = eval(raw_input("Maximum distance a 3' docking strand can be located away from a docking site [suggestion: 5-6nm]: "))# suggestion: 6 nm, max distance of a docking strand from the docking site
min_num_docking_strands = eval(raw_input("Minimum number of docking strands at each docking site [suggestion: 4-6]:")) #suggestion:5, the minimum number of docking strands at each docking site

#######################
### PREPARE FILE ######
#######################


caDNAno_file = open(json_name, 'r') # This opens the file specified by the user
file_contents = eval(caDNAno_file.read()) 
caDNAno_file.close()
strands = file_contents["vstrands"] # This selects all the strand data in the JSON file
f.change_color_of_all_strands(strands,all_staples_color) # Changes color of all staples to black. This will make it easier to spot the docking strands later on

#############################
### BEGIN ANALYZING FILE ####
#############################

print 
print "Got it."
print 
print "Step 1: Select the orientation of the nanostructure."
print 

# Creates the location matrix which has the number, location (row and column) and beginning and ending staple for each strand
location = f.location_matrix(strands)

##################################################################
#### SELECTING WHICH FACE (XY, XZ, YZ) TO VIEW THE STRUCTURE #####
##################################################################

# Calculates the surface area of the nanostructure projected onto the xy, xz, and yz planes 
areas = [f.area_plane(f.xy_points(location,dia)), f.area_plane(f.xz_points(location,dia,len_bp)), f.area_plane(f.yz_points(location,dia,len_bp))]

# Finds the area with the largest 
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

##################################################################
#### SELECTING THE FACE TO ATTACH TO THE STEPTAVIDIN SURFACE #####
##################################################################

staples = f.staple_positions(strands,location) # Finds all staples in the entire nanostructure

# The user is asked to select an appropriate attachment face.
point_face = eval(raw_input("Enter a point on the nanostructure face that should be attached to surface. Use the following format [strand number, staple position]: "))

f.streptavidin_staples(suggested_plane,point_face,location,strands,dia,len_bp)

print "Great. That seemed to work. All of the staple strands that should be biotinylated are now colored magenta."
print
print "Now it is time to select the locations for the docking sites"

####################################################################
#### ANALYZING THE STRUCTURE TO FIND DOCKING SITES AND STRANDS #####
####################################################################

face_pts = f.face_points(staples,location,point_face,suggested_plane,dia,len_bp,min_distance)
print "Found face points"
[first_point,hull_face_pts] = f.first_point(face_pts,staples, max_distance, location, len_bp, dia,min_num_docking_strands)
print "Found first site"

[suggested_sites,enough_sites,possible_combos]=f.remaining_docking_sites(hull_face_pts,first_point,staples,max_distance,location,len_bp,dia,min_num_docking_strands,min_distance)

print "Found suggested sites"

[a,b,c]=f.find_docking_strands(suggested_sites,strands,staples,max_distance,location,len_bp,dia,docking_site_color_1,docking_site_color_2,docking_site_color_3)

print "Found docking strands"

strands=f.ghost_strands(a,b,c,location,strands,three_p_length)



#############
## RESULTS ##
#############

print
print "Results:"
print
print "Three docking sites have been found. They are:"
print "Docking site 1: Strand " + str(f.strandindex_to_strandnum(location,suggested_sites[0][0])) + " at position " + str(suggested_sites[0][1])
print "Docking site 2: Strand " + str(f.strandindex_to_strandnum(location,suggested_sites[1][0])) + " at position " + str(suggested_sites[1][1])
print "Docking site 3: Strand " + str(f.strandindex_to_strandnum(location,suggested_sites[2][0])) + " at position " + str(suggested_sites[2][1])
print
print "Distance between site 1 and 2: " + str(suggested_sites[3][0]) + "nm"
print "Distance between site 2 and 3: " + str(suggested_sites[3][1]) + "nm"
print "Distance between site 3 and 1: " + str(suggested_sites[3][2]) + "nm"
print
print "Potential 3' strands were found at each site. These strands are all within " + str(max_distance) + "nm of the docking site."
print
print "Site 1: " + str(len(a)) + " docking strands colored red."
print "Site 1 docking strands [strand number, 3' position]:"
for i in a:
	print [f.strandindex_to_strandnum(location,i[0]), i[1]]
print
print "Site 2: " + str(len(b)) + " docking strands colored green." 
print "Site 2 docking strands [strand number, 3' position]:"
for i in b:
	print [f.strandindex_to_strandnum(location,i[0]), i[1]]
print
print "Site 3: " + str(len(c)) + " docking strands colored blue."
print "Site 3 docking strands [strand number, 3' position]:"
for i in c:
	print [f.strandindex_to_strandnum(location,i[0]), i[1]]
print

print "Creating new JSON file..."

###########################
## WRITING THE JSON FILE ##
###########################

## Writing JSON file
file_contents_new = file_contents
file_contents_new["vstrands"] = strands
new_name = json_name[:-5] + "_new" + ".json"
with open(new_name, 'w') as outfile:
  json.dump(file_contents_new, outfile)

print "Complete. The new file is " + new_name + "."


###############
## ALL DONE! ##
###############


