import functions as f
import json

# Parametric variables of the helix
bps_turn = 10.5 # Number of base pairs per turn of the double helix
len_bp = .332 # nm, length of the helix per bp
dia = 2.0 # nm, diameter of the duplex

# Get the file from the user
json_name = "hexagon_prism_final_barcode_7.json"

# Open and read the file
caDNAno_file = open(json_name, 'r')
file_contents = eval(caDNAno_file.read())
caDNAno_file.close()

strands = file_contents["vstrands"]

location = f.location_matrix(strands)

print f.find_threeprime_staples(strands)
print len(f.find_threeprime_staples(strands))

strand_index = 0
current_location = 23
while strand_index != None:
	old_strand_index = strand_index
	old_location = current_location
	strand_index = f.strandnum_to_strandindex(location,strands[old_strand_index]["stap"][old_location][0])
	current_location = strands[old_strand_index]["stap"][old_location][1]
print old_strand_index
print old_location


print f.strandnum_to_strandindex(location,55)