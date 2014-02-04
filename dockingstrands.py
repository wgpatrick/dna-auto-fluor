# This script prepares a DNA nanostructure for DNA-PAINT imaging.
# The script inserts docker staple strands into the nanostructure design
# Originally written by Will Patrick - moonshot@mit.edu

# Open and read the file
CDFilename = "/users/willpatrick/Downloads/caDNAnoTrackDesign.json"
caDNAno_file = open(CDFilename, 'r')
file_contents = eval(caDNAno_file.read())
caDNAno_file.close()

strands = file_contents["vstrands"]


num_strands = len(strands)
print strands[0]['scaf'][0]
for strand in strands:
	print start_staple
	end_staple= []
	print str(strand['row']) + ', ' + str(strand['col']) + ', ' + str(len(strand['scaf']))


