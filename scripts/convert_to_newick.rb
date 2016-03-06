# Michael Matschiner, 2016-03-06

# This script converts a tree file in NEXUS format with translation
# table into plain newick format. The input file will be overwritten
# with the output.

# This script should be run e.g. with
# ruby convert_to_newick summary.tre
# where 'summary.tre' should be replaced with the actual path to the
# tree file that should be converted.

# Check if command line arguments are provided, and print help
# text if they aren't.
if ARGV == [] or ARGV.include?("-h") or ARGV.include?("--help")
	puts
	puts "convert_to_newick.rb"
	puts
	puts "This script converts a tree file in NEXUS format with translation"
	puts "table into plain newick format. The input file will be overwritten"
	puts "with the output."
	puts
	puts "This script should be run e.g. with"
	puts "ruby convert_to_newick summary.tre"
	puts "where 'summary.tre' should be replaced with the actual path to the"
	puts "tree file that should be converted."
	exit
end

# Get the command line arguments.
input_file_name = ARGV[0]
output_file_name = input_file_name

# Read the input file.
input_file = File.open(input_file_name)
input_lines = input_file.readlines
input_file.close
translation_number = []
translation_id = []
in_translation = false
tree_string = ""
input_lines.each do |l|
	if l.strip.downcase == "translate"
		in_translation = true
	elsif l.strip == ";"
		in_translation = false
	elsif in_translation
		translation_number << l.strip.split[0]
		translation_id << l.strip.split[1].chomp(",")
	elsif l.strip[0..3].downcase == "tree"
		tree_string = l.split(" = ")[1].strip.chomp(";").chomp(":0.0").gsub(/\[.+?\]/,"")
	end
end

# Reverse translate the taxon ids in the tree string.
translation_number.size.times do |x|
	tree_string.sub!("(#{translation_number[x]}:","(#{translation_id[x]}:")
	tree_string.sub!(",#{translation_number[x]}:",",#{translation_id[x]}:")
end

# Use the tree string to overwrite the input file.
output_file = File.open(output_file_name,"w")
output_file.write("#{tree_string};\n")
output_file.close
