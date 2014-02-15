#!/usr/bin/python
# -*- coding: ascii -*-

# List of functions used in AutoFluor 0.1
# This script prepares a DNA nanostructure for super resolution fluorescent microscopy.
# Originally written by Will Patrick - moonshot@mit.edu - MIT Media Lab
# Feb 2014

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

	elif suggested_plane == "xz":
		y_pos = row2y(location_matrix[strand_num][1],location_matrix[strand_num][2],dia)
		for strand in strands:
			y_strand = row2y(strand['row'],strand['col'],dia)
			if y_strand == y_pos:
				for five_prime_end in strand['stap_colors']:
					five_prime_end[1] = 9109606
	else:
		x_pos = col2x(location_matrix[strand][2])
		for strand in strands:
			x_strand = col2x(strand['col'],dia)
			if x_strand == x_pos:
				for five_prime_end in strand['stap_colors']:
					five_prime_end[1] = 9109606
				
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

## Finds the location of the centroid on the points of face

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

# Returns a list of 3' strands within a certain distance of a specific point on the nanostructure. 
# Arguments are: the position of interest (strand_index & staple pos), the list of all staples computed by staple_positions(), 
# the max_distance allowed, location matrix and length per base pair and the helical diamater


def nearby_threeprime(strand_index,staple_pos,all_staples,max_distance,location_matrix,len_bp,dia):
	nearby_staples=[]
	for staple in all_staples:
		if p2pdist(strand_index,staple_pos,staple[0],staple[1],location_matrix,len_bp,dia)[6] < max_distance:
			nearby_staples.append(staple)
	return nearby_staples

### This function finds the first docking site on the nanostructure

def first_point(face_points,staples, max_distance, location, len_bp, dia, min_num_docking_strands):
	face_centroid_hull=centroid_hull(face_points)

	maximum = max([i[4] for i in face_centroid_hull[0]])
	index = [i[4] for i in face_centroid_hull[0]].index(maximum)
	sorted_face_hull = sorted(face_centroid_hull[0], key=itemgetter(4),reverse=True)
	
	for point in sorted_face_hull:
		num_nearby_docking_strands = len(nearby_threeprime(point[0],point[1],staples,max_distance,location,len_bp,dia))
		if num_nearby_docking_strands >= min_num_docking_strands:
			first_point = point
			break
	if first_point ==[]:
		print "Could not find a first docking site with 5 or more nearby 3' staple strands. Try increasing the minimum distance between the fluorophore docking strand and the site"
	return [first_point,sorted_face_hull]		


## This finds the remaining docking sites

def remaining_docking_sites(hulled_face_points,first_point,all_staples,max_distance,location_matrix,len_bp,dia,min_num_docking_strands,min_distance):
	sites_with_enough_nearby_docking_strands =[]
	sites_with_enough=[]
	for site in hulled_face_points:
		threeprimes=nearby_threeprime(site[0],site[1],all_staples,max_distance,location_matrix,len_bp,dia)
		if len(threeprimes) >= min_num_docking_strands:
			sites_with_enough.append(site)

	possible_combos=[]
	for i in range(0,len(sites_with_enough),1):
		for j in range(0,len(sites_with_enough),1):
			d1 = math.sqrt( (sites_with_enough[i][2]-first_point[2])**2 + (sites_with_enough[i][3]-first_point[3])**2 )
			d2 = math.sqrt( (sites_with_enough[i][2]-sites_with_enough[j][2])**2 + (sites_with_enough[i][3]-sites_with_enough[j][3])**2 )
			d3 = math.sqrt( (sites_with_enough[j][2]-first_point[2])**2 + (sites_with_enough[j][3]-first_point[3])**2 )	
			if ((d3 > min_distance) & (d2 > min_distance) & (d1 > min_distance)):			
				site_1 = first_point[0:2]
				site_2 = sites_with_enough[i][0:2]
				site_3 = sites_with_enough[j][0:2]
				d_perim = d1+d2+d3
				if not ( [site_1,site_3,site_2,[d3,d2,d1],d_perim] in possible_combos):
					possible_combos.append([site_1,site_2,site_3,[d1,d2,d3],d_perim])
	sorted_combos = sorted(possible_combos,key=itemgetter(4),reverse=True)
	suggested_sites = sorted_combos[0]

	return [suggested_sites,sites_with_enough,sorted_combos]

# find_docking_strands() formula does two things. (1) It changes the color of all the strands at sites 1, 2 and 3. (2) It outputs all the docking strands.

def find_docking_strands(suggested_sites,strands,all_staples,max_distance,location_matrix,len_bp,dia,color1,color2,color3):
	A = suggested_sites[0]
	B = suggested_sites[1]
	C = suggested_sites[2]

	A_docking_strands = nearby_threeprime(A[0],A[1],all_staples,max_distance,location_matrix,len_bp,dia)
	B_docking_strands = nearby_threeprime(B[0],B[1],all_staples,max_distance,location_matrix,len_bp,dia)
	C_docking_strands = nearby_threeprime(C[0],C[1],all_staples,max_distance,location_matrix,len_bp,dia)


	for i in A_docking_strands:
		for staple in strands[i[2]]["stap_colors"]:
			if staple[0] == i[3]:
				staple[1] = color1 
	for i in B_docking_strands:
		for staple in strands[i[2]]["stap_colors"]:
			if staple[0] == i[3]:
				staple[1] = color2 
	for i in C_docking_strands:
		for staple in strands[i[2]]["stap_colors"]:
			if staple[0] == i[3]:
				staple[1] = color3 

	return [A_docking_strands,B_docking_strands,C_docking_strands]

#### Ghost strands

def ghost_strands(a_staples,b_staples,c_staples,location,strands,three_p_length):
	all_strands = strands
	total_num_strands = len(location)
	total_len_strand = len(all_strands[0]["stap"])
	staple_extension =[a_staples,b_staples,c_staples]
	
	# making ghost strands
	new_strands=[]
	for i in range(0,len(staple_extension),1):
		new_strand={}
		new_strand['row'] = 0
		total_num_strands =len(all_strands)+len(new_strands)
		if total_num_strands % 2 == 0:
			new_strand['num'] = total_num_strands + 1
		else:
			new_strand['num'] = total_num_strands - 1
		new_strand['col'] = i
		new_strand['scaf']=[]
		for a in range(0,total_len_strand,1):
			new_strand['scaf'].append([-1,-1,-1,-1])
		new_strand['stap']=[]
		for a in range(0,total_len_strand,1):
			new_strand['stap'].append([-1,-1,-1,-1])
		new_strand['stap_colors']=[]
		new_strand['skip']=[]
		for j in range(0,total_len_strand,1):
			new_strand['skip'].append(0)
		new_strand['loop']=[]
		for j in range(0,total_len_strand,1):
			new_strand['loop'].append(0)
		new_strand['stap_Loop']=[]
		new_strand['scafLoop']=[]
		new_strands.append(new_strand)
	all_strands.extend(new_strands)

	new_location = location_matrix(all_strands)

	# adding staple extensions to ghost strands

	num_a = len(a_staples)
	staple_locs =[]	
	for i in range(0,len(staple_extension),1):
		staples = staple_extension[i]
		position = 0
		staple_locs.append([])	
		if  (len(all_strands) - 2 + i)  % 2 == 1:
			for staple in staples:
				staple_locs[i].append([strandindex_to_strandnum(new_location,len(all_strands)-3+i),position])
				all_strands[-3+i]['stap'][position] = [staple[0],staple[1],all_strands[-3+i]['num'], position+1,]
				for bp in range(0,three_p_length-2,1):
					position += 1
					all_strands[-3+i]['stap'][position] = [all_strands[-3+i]['num'], position - 1, all_strands[-3+i]['num'], position + 1]
				position += 1
				all_strands[-3+i]['stap'][position] = [all_strands[-3+i]['num'], position - 1,-1,-1]
				position += 3
		else:
			position = len(staples) * (three_p_length + 3)
			for staple in staples:
				staple_locs[i].append([strandindex_to_strandnum(new_location,len(all_strands)-3+i),position])
				all_strands[-3+i]['stap'][position] = [staple[0],staple[1],all_strands[-3+i]['num'], position-1]
				for bp in range(0,three_p_length-2,1):
					position += -1
					all_strands[-3+i]['stap'][position] = [all_strands[-3+i]['num'], position + 1, all_strands[-3+i]['num'], position - 1]
				position += -1
				all_strands[-3+i]['stap'][position] = [all_strands[-3+i]['num'], position + 1,-1,-1]
				position += -3

	# # adding links to the ghost strands
	print staple_locs
	for site in staple_extension:
		for staple in site:
			strand_index = staple[0]
			position = staple[1]
			all_strands[strand_index]["stap"][position][2] = staple_locs[staple_extension.index(site)][site.index(staple)][0]
			all_strands[strand_index]["stap"][position][3] = staple_locs[staple_extension.index(site)][site.index(staple)][1]
	

	return all_strands


#### This formula changes the color of staple strands in the caDNAno file to the specific color

def change_color_of_all_strands(strands,color):
	for strand in strands:
		for staple in strand["stap_colors"]:
			staple[1] = color

##### ##### ##### ##### ##### ##### ##### ##### #####
#### Find the strand index given the strand number ##
##### ##### ##### ##### ##### ##### ##### ##### ##### 

def strandnum_to_strandindex(location_matrix,num):
	for i in range(0,len(location_matrix),1):
		if location_matrix[i][0]==num:
			return i

##### ##### ##### ##### ##### ##### ##### ##### ##### ###
##### Find the strand index given the strand number ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ###

def strandindex_to_strandnum(location_matrix,index):
	return location_matrix[index][0]



