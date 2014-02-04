# This script prepares a DNA nanostructure for DNA-PAINT imaging.
# The script inserts docker staple strands into the nanostructure design
# Originally written by Will Patrick - moonshot@mit.edu

# Open and read the file
CDFilename = "/users/willpatrick/Downloads/caDNAnoTrackDesign.json"
caDNAno_file = open(CDFilename, 'r')
file_contents = eval(caDNAno_file.read())
caDNAno_file.close()

name_line = file_contents["vstrands"]

print name_line[1]['col']