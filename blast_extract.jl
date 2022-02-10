# v.1.0

# REQUIRED FILE STRUCTURE

# . 					
# ├── /reference       	Directory containing blast search subject sequences in fasta format
# ├── /blast_out       	Directory containing *_out.txt blast output in outfmt6 format
# └── blast_extract.jl 	

# INFO

# This script was designed to extract and return query blast hits from a reference genome in a reciprocal blast annotation pipeline. \
# Sequences are output in fasta format after optional filtering steps to handle negative sense hits (negative_sense()) \
# and merge those within a specified number of basepairs (sseqid_merge(n)). 

# referenceIO() and blastIO() default to the directory the blast_extract.jl script is located in ("."). \
# A argument can be added in the command line along with referenceIO(ARGS[1]) and blastIO(ARGS[1]) (ex. "julia blast_extract.jl YOUR_PATH/") if \
# the blast_extract.jl file is located elsewhere

# Currently functional in blast searches with nucleotide subject

using JuliaDB
function referenceIO(dir::AbstractString = ".") #DONE
	global reference = Dict{String, Any}()
	for file in readdir("$dir/reference")
		file_name = SubString(file, 1, Integer(findfirst('.', file))-1)
		println("Starting referenceIO for $file_name")
		temp_dict = Dict{String, String}()
		temp_key = ""
		temp_value = ""
		for line in eachline("$dir/reference/$file")
			if startswith(line, ">") === true
				if findfirst(' ', line) === nothing
					temp_key = SubString(line, 2:length(line))
				else
					temp_index = findfirst(' ', line)
					temp_key = SubString(line, 2:(temp_index-1))
				end
			else
				temp_value = temp_value * line
				merge!(temp_dict, Dict{String, String}(temp_key => temp_value))
			end
		end
		merge!(reference, Dict(file_name => temp_dict))
	end
	return reference
end

function blastIO(dir::AbstractString = ".") #DONE
	global blast_out = Dict{String, IndexedTable}()
	for file in readdir("$dir/blast_out")
		file_name = SubString(file, 1, Integer(findfirst('.', file))-1)
		println("Starting blastIO for $file_name")
		merge!(blast_out, Dict(file_name => loadtable("$dir/blast_out/$file", spacedelim = true, header_exists = false,
			colnames = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])))
		blast_out[file_name] = insertcols(blast_out[file_name], 1, :index => range(1, length(blast_out[file_name])))
	end
	return blast_out
end

function negative_sense() #DONE
	for key in keys(reference)
		key_bridge = SubString(key, 1, Integer(findfirst('_', key))-1) # Rewrite so it works with multiple out and one db
		try
			for i in range(1, length(blast_out["$key_bridge"*"_out"]))
				if blast_out["$key_bridge"*"_out"][i][:sstart] > blast_out["$key_bridge"*"_out"][i][:send]
					temp_table = table([i], [blast_out["$key_bridge"*"_out"][i][:qseqid]], [blast_out["$key_bridge"*"_out"][i][:sseqid]], [blast_out["$key_bridge"*"_out"][i][:pident]], 
					[blast_out["$key_bridge"*"_out"][i][:length]], [blast_out["$key_bridge"*"_out"][i][:mismatch]], [blast_out["$key_bridge"*"_out"][i][:gapopen]], [blast_out["$key_bridge"*"_out"][i][:qstart]], 
					[blast_out["$key_bridge"*"_out"][i][:qend]], [blast_out["$key_bridge"*"_out"][i][:send]], [blast_out["$key_bridge"*"_out"][i][:sstart]], [blast_out["$key_bridge"*"_out"][i][:evalue]], [blast_out["$key_bridge"*"_out"][i][:bitscore]],  
					names = [:index, :qseqid, :sseqid, :pident, :length, :mismatch, :gapopen, :qstart, :qend, :sstart, :send, :evalue, :bitscore])
					blast_out["$key_bridge"*"_out"] = merge(blast_out["$key_bridge"*"_out"], temp_table, pkey = :index)
				else
					continue
				end
				blast_out["$key_bridge"*"_out"] = filter(r -> r.sstart < r.send, blast_out["$key_bridge"*"_out"])
			end
		catch
			println("*** A corresponding output file was not found for $key.fasta ***")
		end
	end
end

function sseqid_merge(n::Integer) #DONE
	global sseqid_parse = Dict{String, Any}()
	temp_counter = 0
	for key in keys(reference)
		temp_dict = Dict{String, Vector}()
		key_bridge = SubString(key, 1, Integer(findfirst('_', key))-1) # Rewrite so it works with multiple out and one db
		for sid in Set(select(blast_out["$key_bridge"*"_out"], :sseqid))
			temp_table_s = filter(val -> (val.sseqid == sid), blast_out["$key_bridge"*"_out"])
			for qid in Set(select(blast_out["$key_bridge"*"_out"], :qseqid))
				temp_table_sq = filter(val -> (val.qseqid == qid), temp_table_s)
				temp_table_sq = sort(temp_table_sq, :send)
				if length(temp_table_sq) == 1
					if !haskey(temp_dict, sid)
						merge!(temp_dict, Dict{String, Array}(sid => [temp_table_sq[1][:sstart], temp_table_sq[1][:send]]))
					else
						append!(temp_dict[sid], [temp_table_sq[1][:sstart], temp_table_sq[1][:send]])
					end
				elseif length(temp_table_sq) == 2
					if temp_table_sq[1][:send]+n >= temp_table_sq[2][:sstart]
						if !haskey(temp_dict, sid)
							merge!(temp_dict, Dict{String, Array}(sid => [temp_table_sq[1][:sstart], temp_table_sq[2][:send]]))
						else
							append!(temp_dict[sid], [temp_table_sq[1][:sstart], temp_table_sq[2][:send]])
						end
					else
						if !haskey(temp_dict, sid)
							merge!(temp_dict, Dict{String, Array}(sid => [temp_table_sq[1][:sstart], temp_table_sq[1][:send]]))
						else
							append!(temp_dict[sid], [temp_table_sq[1][:sstart], temp_table_sq[1][:send]])
						end
					end
				else
					temp_send = nothing
					for i in range(1, length(temp_table_sq))
						if i+1 > length(temp_table_sq)
							append!(temp_dict[sid], [temp_table_sq[i][:sstart], temp_table_sq[i][:send]])
							break
						elseif temp_table_sq[i][:send]+n >= temp_table_sq[i+1][:sstart]
							temp_send = temp_table_sq[i+1][:send]
							temp_counter += 1
						elseif temp_table_sq[i][:send]+n < temp_table_sq[i+1][:sstart]
							if temp_send === nothing
								if !haskey(temp_dict, sid)
									merge!(temp_dict, Dict{String, Array}(sid => [temp_table_sq[i][:sstart], temp_table_sq[i][:send]]))
								else					
									append!(temp_dict[sid], [temp_table_sq[i][:sstart], temp_table_sq[i][:send]])
								end
							else
								try
									append!(temp_dict[sid], [temp_table_sq[i-temp_counter][:sstart], temp_send])
									temp_send = nothing
									temp_counter = 0
								catch
									append!(temp_dict[sid], [temp_table_sq[1][:sstart], temp_send])
									temp_send = nothing
									temp_counter = 0
								end
							end
						end
					end
				end
			end
			merge!(sseqid_parse, Dict(key => temp_dict))
		end
	end
	return sseqid_parse
end

using IterTools
function parseIO() #DONE
	for file in keys(sseqid_parse)
		println("Parsing sequences for hits found in $file")
		file_bridge = SubString(file, 1, Integer(findfirst('_', file))-1)
		io = open("$file_bridge"*"_parsed.fa", "w")
		for id in keys(sseqid_parse[file])
			for (i, (sstart, send)) in enumerate(partition(sseqid_parse[file][id], 2, 2))
				temp_substring = SubString(reference[file][id], sstart:send)
				write(io, ">$id"*"_$i\n") + write(io, "$temp_substring\n")
			end
		end
		close(io)
	end
end

referenceIO() #upload subject sequences
println("referenceIO completed!")
blastIO() #upload blast output
println("blastIO completed!")
negative_sense() #filter hits to ensure sstart < send
println("Successfully filtered for negative sense hits!")
sseqid_merge(2000) #merge hits within n basepairs and extract hits
println("Successfully merged hits within the specified number of basepairs!")
parseIO() #output extracted hits