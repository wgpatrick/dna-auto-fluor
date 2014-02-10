#!/usr/bin/python
# -*- coding: ascii -*-

import math
import numpy as np
from pyhull.convex_hull import ConvexHull

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

def p2pdist(strand_1,position_1,strand_2,position_2, location_matrix, len_bp):
	row_1 = location_matrix[strand_1][1]
	col_1 = location_matrix[strand_1][2]
	row_2 = location_matrix[strand_2][1]
	col_2 = location_matrix[strand_2][2]
	x = col2x(col_2)-col2x(col_1)
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
		print suggested_plane
		z_pos = staple2z(strand_num,len_bp)
		for strand in strands:
			
		print "z =" +  str(z_pos)
		for 

	elif suggested_plane == "xz":
		y_pos = row2y(location_matrix[strand_num][1],location_matrix[strand_num][2],dia)
		print "y =" + str(y_pos)
		for strand in strands:
			y_strand = row2y(strand['row'],strand['col'],dia)
			if y_strand == y_pos:
				for five_prime_end in strand['stap_colors']:
					five_prime_end[1] = 3333333
				print "Strand " + str(strand['num']) + " has " + str(len(strand['stap_colors'])) + " staple strands that are now magenta."
	else:
		x_pos = col2x(location_matrix[strand][2])
		print "x =" + str(x_pos)
		for strand in strands:
			x_strand = col2x(strand['col'],dia)
			if x_strand == x_pos:
				for five_prime_end in strand['stap_colors']:
					five_prime_end[1] = 3333333
				print "Strand " + str(strand['num']) + " has " + str(len(strand['stap_colors'])) + " staple strands that are now magenta."

	print suggested_plane
