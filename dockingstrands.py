# This script prepares a DNA nanostructure for DNA-PAINT.
# Originally written by Will Patrick - moonshot@mit.edu

# List of things this script cannot do:
# Cannot handle helices that have have a break

import math
from pyhull.convex_hull import ConvexHull


# Open and read the file
CDFilename = "caDNAnoTrackDesign.json"
caDNAno_file = open(CDFilename, 'r')
file_contents = eval(caDNAno_file.read())
caDNAno_file.close()
strands = file_contents["vstrands"]

# Parametric variables of the helix
bps_turn = 10.5 # Number of base pairs per turn of the double helix
len_bp = .332 # nm, length of the helix per bp
dia = 2.0 # nm, diameter of the duplex


# Make a location matrix which contains column num, row, column, beginning strand and staple strand information
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
print location

# Calculate the distance between two points in the nanostructure.
# This function outputs the distance in the x, y, z, xy, xz, yz and xyz components.

def p2pdist(strand_1,position_1,strand_2,position_2):
	row_1 = location[strand_1][1]
	col_1 = location[strand_1][2]
	row_2 = location[strand_2][1]
	col_2 = location[strand_2][2]
	x = (col_2-col_1) * dia * math.sqrt(3) / 2.0
	y = (row_2-row_1) * dia * 1.5 + (dia/4.0 - (col_2%2 * dia/2)) - (-dia/4.0 + col_1%2 *dia/2)
	z = (position_2 - position_1) * len_bp
	xy = math.sqrt(x**2 + y**2)
	xz = math.sqrt(x**2 + z**2)
	yz = math.sqrt(y**2 + z**2)
	xyz = math.sqrt(xy**2 + z**2)
	return [x,y,z,xy,xz,yz,xyz]



a = p2pdist(0,10,8,100)
print a
