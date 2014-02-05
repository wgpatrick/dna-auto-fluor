# This script prepares a DNA nanostructure for DNA-PAINT.
# Originally written by Will Patrick - moonshot@mit.edu

# List of things this script cannot do:
# Cannot handle helices that have have a break

import math
import numpy as np
from pyhull.convex_hull import ConvexHull

# Parametric variables of the helix
bps_turn = 10.5 # Number of base pairs per turn of the double helix
len_bp = .332 # nm, length of the helix per bp
dia = 2.0 # nm, diameter of the duplex

# 
# FUNCTIONS
#

#Make a location matrix which contains column num, row, column, beginning strand and staple strand information

def location_matrix(strands):
	location=[]
	for strand in strands:
		row = strand['row']
		col = strand['col']
		num = strand['num']
		i = 0
		# Finds first strand of the 
		for scaf in strand['scaf']:
			if scaf[0] != -1:
				beginning_staple = i
				break
			else:
				i += 1
		j = 0
		for scaf in reversed(strand['scaf']):		
			if scaf[0] != -1:
				end_staple = len(strand['scaf']) - j - 1
				break
			else:
				j += 1
		location.append([num,row,col,beginning_staple,end_staple])
	return location

# Calculate the distance between two points in the nanostructure.
# This function outputs the distance in the x, y, z, xy, xz, yz and xyz components.

def p2pdist(strand_1,position_1,strand_2,position_2):
	row_1 = location[strand_1][1]
	col_1 = location[strand_1][2]
	row_2 = location[strand_2][1]
	col_2 = location[strand_2][2]
	x = col2x(col_2)-col2x(col_1)
	y = row2y(row_2,col_2)-row2y(row_1,col_1)
	z = (position_2 - position_1) * len_bp
	xy = math.sqrt(x**2 + y**2)
	xz = math.sqrt(x**2 + z**2)
	yz = math.sqrt(y**2 + z**2)
	xyz = math.sqrt(xy**2 + z**2)
	return [x,y,z,xy,xz,yz,xyz]

#### This calculates the x position of a strand given it's row and column in the hexagonal lattice
def row2y(row, col):
	if row % 2 ==0:
		return (row * dia * 1.5 + (col%2 * dia/2))
	else:
		return (row * dia * 1.5) + (dia/2.0 - col%2 * dia/2.0)

# This calculates the y positon of a strand given it's column number in the hexagonal lattice
def col2x(col):
	return (col * dia * math.sqrt(3) / 2.0) + dia/2

def staple2z(staple_pos):
	z = staple_pos*len_bp
	return z

# calculate points on xy_plane
def xy_points(location):
	xy_points=[]
	for strand in location:
		x_point = col2x(strand[2])
		y_point = row2y(strand[1],strand[2])
		xy_points.append([x_point,y_point])
	return xy_points

def xz_points(location):
	xz_points=[]
	for strand in location:
		x_point = col2x(strand[2])
		z_point_1 = staple2z(strand[3])
		z_point_2 = staple2z(strand[4])
		xz_points.append([x_point,z_point_1])
		xz_points.append([x_point,z_point_2])
	return xz_points

def yz_points(location):
	yz_points=[]
	for strand in location:
		y_point = row2y(strand[1],strand[2])
		z_point_1 = staple2z(strand[3])
		z_point_2 = staple2z(strand[4])
		yz_points.append([y_point,z_point_1])
		yz_points.append([y_point,z_point_2])
	return yz_points

def area_plane(points):
	convex_hull = ConvexHull(points)
	area = 0.0
	vertices = convex_hull.vertices
	hull_points = convex_hull.points
	for vertex in vertices:
		area += hull_points[vertex[0]][0] * hull_points[vertex[1]][1] - hull_points[vertex[0]][1] * points[vertex[1]][0]
	area = area / 2.0
	return area

## MAIN 

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
json = raw_input("Enter the name of the caDNAno file: ")

# Open and read the file
CDFilename = json
caDNAno_file = open(CDFilename, 'r')
file_contents = eval(caDNAno_file.read())
caDNAno_file.close()
strands = file_contents["vstrands"]

print 
print "Got it."
print 
print "Step 1: Select the orientation of the nanostructure."
print 
# Creates the location matrix which has the number, location (row and column) and beginning and ending staple for each strand
location = location_matrix(strands)

# Calculates the surface area of the nanostructure projected onto the xy, xz, and yz planes 
areas = [area_plane(xy_points(location)), area_plane(xz_points(location)), area_plane(yz_points(location))]

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



print point_face



