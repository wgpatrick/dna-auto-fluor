# This script prepares a DNA nanostructure for DNA-PAINT imaging.
# The script inserts docker staple strands into the nanostructure design
# Originally written by Will Patrick - moonshot@mit.edu

# Open and read the file
CDFilename = "caDNAnoTrackDesign.json"
caDNAno_file = open(CDFilename, 'r')
file_contents = eval(caDNAno_file.read())
caDNAno_file.close()

strands = file_contents["vstrands"]


num_strands = len(strands)


# Make a location matrix which contains column num, row, column, beginning strand and staple strand information
location=[]
for strand in strands:
	row = strand['row']
	col = strand['col']
	num = strand['num']
	i = 0
	# Finds first strand 
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
	
   
