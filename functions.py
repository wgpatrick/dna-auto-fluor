#!/usr/bin/python
# -*- coding: ascii -*-

import math
import numpy as np
from pyhull.convex_hull import ConvexHull
from operator import itemgetter, attrgetter

#Make a location matrix which contains column num, row, column, beginning strand and staple strand information

def location_matrix(strands):
	location_matrix=[]
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
		location_matrix.append([num,row,col,beginning_staple,end_staple])
	return location_matrix

# Calculate the distance between two points in the nanostructure.
# This function outputs the distance in the x, y, z, xy, xz, yz and xyz components.

def p2pdist(strand_index_1,position_1,strand_index_2,position_2, location_matrix, len_bp, dia):
	row_1 = location_matrix[strand_index_1][1]
	col_1 = location_matrix[strand_index_1][2]
	row_2 = location_matrix[strand_index_2][1]
	col_2 = location_matrix[strand_index_2][2]
	x = col2x(col_2,dia)-col2x(col_1,dia)
	y = row2y(row_2,col_2,dia)-row2y(row_1,col_1,dia)
	z = (position_2 - position_1) * len_bp
	xy = math.sqrt(x**2 + y**2)
	xz = math.sqrt(x**2 + z**2)
	yz = math.sqrt(y**2 + z**2)
	xyz = math.sqrt(xy**2 + z**2)
	return [x,y,z,xy,xz,yz,xyz]

#### This calculates the x position of a strand given it's row and column in the hexagonal lattice
def row2y(row, col, dia):
	if row % 2 ==0:
		return (row * dia * 1.5 + (col%2 * dia/2))
	else:
		return (row * dia * 1.5) + (dia/2.0 - col%2 * dia/2.0)

# This calculates the y positon of a strand given it's column number in the hexagonal lattice
def col2x(col, dia):
	return (col * dia * math.sqrt(3) / 2.0) + dia/2

def staple2z(staple_pos,len_bp):
	z = staple_pos*len_bp
	return z

# calculate points on xy_plane
def xy_points(location_matrix,dia):
	xy_points=[]
	for strand in location_matrix:
		x_point = col2x(strand[2],dia)
		y_point = row2y(strand[1],strand[2],dia)
		xy_points.append([x_point,y_point])
	return xy_points

def xz_points(location_matrix,dia,len_bp):
	xz_points=[]
	for strand in location_matrix:
		x_point = col2x(strand[2],dia)
		z_point_1 = staple2z(strand[3],len_bp)
		z_point_2 = staple2z(strand[4],len_bp)
		xz_points.append([x_point,z_point_1])
		xz_points.append([x_point,z_point_2])
	return xz_points

def yz_points(location_matrix,dia,len_bp):
	yz_points=[]
	for strand in location_matrix:
		y_point = row2y(strand[1],strand[2],dia)
		z_point_1 = staple2z(strand[3],len_bp)
		z_point_2 = staple2z(strand[4],len_bp)
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

## Biotin tagging staple strands.
# Find staple strands on the attachment face and change their color.
# Will have to be different for different planes

def streptavidin_staples(suggested_plane, loc, location_matrix, strands,dia,len_bp):
	strand_num = loc[0]
	staple_loc = loc[1]

	if suggested_plane == "xy":
		i = 0
		for strand in strands:
			for five_prime_end in strand['stap_colors']:
				if abs(five_prime_end[0] - staple_loc) < 7:
					five_prime_end[1] = 9109606
					i =+ 1
		print i + " staple stands are now magenta within 7 bp of the attachment location"

	elif suggested_plane == "xz":
		y_pos = row2y(location_matrix[strand_num][1],location_matrix[strand_num][2],dia)
		for strand in strands:
			y_strand = row2y(strand['row'],strand['col'],dia)
			if y_strand == y_pos:
				for five_prime_end in strand['stap_colors']:
					five_prime_end[1] = 9109606
				print "Strand " + str(strand['num']) + " has " + str(len(strand['stap_colors'])) + " staple strands that are now magenta."
	else:
		x_pos = col2x(location_matrix[strand][2])
		for strand in strands:
			x_strand = col2x(strand['col'],dia)
			if x_strand == x_pos:
				for five_prime_end in strand['stap_colors']:
					five_prime_end[1] = 9109606
				print "Strand " + str(strand['num']) + " has " + str(len(strand['stap_colors'])) + " staple strands that are now magenta."

	print suggested_plane

# Finds the location of the 5' prime and 3' ends of every staple strand. 
#The output is a list of all staples with [3' strand index, 3' staple position, 5' strand index, 5' staple position]

def staple_positions(strands,location_matrix):
	# Calculates a list which has all of the three prime locations . 
	three_primes=[]
	for strand in strands:
		prime_build=[]
		for i in range(0,len(strand["stap"]),1):
			if (strand["stap"][i][2] == -1) & (strand["stap"][i][0] >= 0):
				prime_build.append(i)
		three_primes.append(prime_build)

	# Finds all of the corresponding 5' ends and builds the list
	staple_positions = []
	for i in range(0,len(three_primes),1):
		for j in range(0,len(three_primes[i]),1):
			strand_index = i
			current_location = three_primes[i][j]
			while strand_index != None:
				old_strand_index = strand_index
				old_location = current_location
				strand_index = strandnum_to_strandindex(location_matrix,strands[old_strand_index]["stap"][old_location][0])
				current_location = strands[old_strand_index]["stap"][old_location][1]
			staple_positions.append([i,three_primes[i][j],old_strand_index,old_location])
	return staple_positions

# finds all of the potential docking sites (places with a staple strand with a 3' end) on the suggested face

def face_points(staple_positions,location_matrix,loc,suggested_plane,dia,len_bp):
	strand_index = loc[0]
	staple_loc = loc[1]

	face_points=[]

	if suggested_plane == "xy":
		
		#Notice this is an imperfect approach for figuring out which "side" of the structure the user selected the point
		
		for i in range(0,len(location_matrix),1):
			x_point = col2x(location_matrix[i][2],dia)
			y_point = row2y(location_matrix[i][1],location_matrix[i][2],dia)
			if abs(staple_loc - staple_start) < abs(staple_loc - staple_end):
				face_points.append([i,location_matrix[i][4],x_point,y_point])
			else:
				face_points.append([i,location_matrix[i][5],x_point,y_point])

	elif suggested_plane == "xz":

		for i in range(0,len(location_matrix),1):
			x_point = col2x(location_matrix[i][2],dia)
			z_point_1 = staple2z(location_matrix[i][3],len_bp)
			z_point_2 = staple2z(location_matrix[i][4],len_bp)
			face_points.append([i,location_matrix[i][3],x_point,z_point_1])
			face_points.append([i,location_matrix[i][4],x_point,z_point_2])
	
	else:

		for i in range(0,len(location_matrix),1):
			y_point = row2y(location_matrix[i][1],location_matrix[i][2],dia)
			z_point_1 = staple2z(location_matrix[i][3],len_bp)
			z_point_2 = staple2z(location_matrix[i][4],len_bp)
			face_points.append([i,location_matrix[i][3],y_point,z_point_1])
			face_points.append([i,location_matrix[i][4],y_point,z_point_2])
	return face_points

## Finds the location of the centroid on the suggest plane

def centroid_hull(face_points):
	points_for_hull = []
	for point in face_points:
		points_for_hull.append([point[2],point[3]])
	convex_hull = ConvexHull(points_for_hull)
	hull_points = convex_hull.points
	hull_vertices = convex_hull.vertices
	for_centroid = [0,0]
	face_points_hull = []
	for i in hull_vertices:
		for_centroid[0] += hull_points[i[0]][0]
		for_centroid[1] += hull_points[i[0]][1]
		face_points_hull.append([face_points[i[0]][0],face_points[i[0]][1],hull_points[i[0]][0],hull_points[i[0]][1]])
	centroid = [for_centroid[0]/float(len(hull_vertices)),for_centroid[1]/float(len(hull_vertices))]

	for point in face_points_hull:
		d_centroid = math.sqrt((point[2] - centroid[0])**2 + (point[3] - centroid[1])**2)
		point.append(d_centroid)

	return [face_points_hull,centroid]
	
## This formula finds the strand index given a number of a strand

def strandnum_to_strandindex(location_matrix,num):
	for i in range(0,len(location_matrix),1):
		if location_matrix[i][0]==num:
			return i

# Returns a list of 3' strands within a certain distance of a specific point on the nanostructure. 
# Arguments are: the position of interest (strand_index & staple pos), the list of all staples computed by staple_positions(), 
# the max_distance allowed, location matrix and length per base pair and the helical diamater


def nearby_threeprime(strand_index,staple_pos,all_staples,max_distance,location_matrix,len_bp,dia):
	nearby_staples=[]
	for staple in all_staples:
		#print p2pdist(strand_index,staple_pos,staple[0],staple[1],location_matrix,len_bp,dia)[6]
		if p2pdist(strand_index,staple_pos,staple[0],staple[1],location_matrix,len_bp,dia)[6] < max_distance:
			nearby_staples.append(staple)
	return nearby_staples

### This function finds the first docking site on the nanostructure

def first_point(face_points,staples, max_distance, location, len_bp, dia):
	face_centroid_hull=centroid_hull(face_points)
	maximum = max([i[4] for i in face_centroid_hull[0]])
	index = [i[4] for i in face_centroid_hull[0]].index(maximum)
	sorted_face_hull = sorted(face_centroid_hull[0], key=itemgetter(4),reverse=True)
	for point in sorted_face_hull:
		num_nearby_docking_strands = len(nearby_threeprime(point[0],point[1],staples,max_distance,location,len_bp,dia))
		if num_nearby_docking_strands >= 5:
			first_point = point
			break
	return first_point		








## Find the locations of the docking sites




