# v.1.2

# CHANGE LOG
# negativesense() function depreciated and merged into parseIO

# REQUIRED FILE STRUCTURE
# . 					
# ├── /reference       	Directory containing blast search subject sequences in fasta format
# ├── /blast_out       	Directory containing *_out.txt blast output in outfmt6 format
# └── blast_extract.jl 	

# INFO
# Updated README to come..

using Dates
function logIO(str::AbstractString)
	if @isdefined log
		println(log, str, now())
		flush(log)
	else
		global log = open(ARGS[1]*"/blast_extract.log", "w")
		println(log, "log file for blast_extract.jl - ", now())
		flush(log)
	end
end

using JuliaDB
logIO("Initiating referenceIO function - ")
function referenceIO(dir::AbstractString = ".") #DONE
	global reference = Dict{String, Any}()
	for file in readdir("$dir/reference")
		temp_counter = 0
		file_name = SubString(file, 1, Integer(findfirst('.', file))-1)
		logIO("Starting referenceIO for $file_name - ")
		temp_dict = Dict{String, String}()
		temp_key = ""
		temp_value = ""
		for line in eachline("$dir/reference/$file")
			if temp_counter <= 6 #ONLY LOADS FIRST 7 SCAFFOLDS
				if startswith(line, ">") === true
					if findfirst(' ', line) === nothing
						temp_key = SubString(line, 2:length(line))
					else
						temp_index = findfirst(' ', line)
						temp_key = SubString(line, 2:(temp_index-1))
					end
					temp_counter = temp_counter + 1
				else
					temp_value = temp_value * line
				end
			else
				break
			end
			merge!(temp_dict, Dict{String, String}(temp_key => temp_value))
		end
		merge!(reference, Dict(file_name => temp_dict))
	end
	return reference
end
logIO("referenceIO function initiated! - ")

logIO("Initiating blastIO function - ")
function blastIO(dir::AbstractString = ".") #DONE
	global blast_out = Dict{String, IndexedTable}()
	for file in readdir("$dir/blast_out")
		file_name = SubString(file, 1, Integer(findfirst('.', file))-1)
		logIO("Starting blastIO for $file_name - ")
		merge!(blast_out, Dict(file_name => loadtable("$dir/blast_out/$file", spacedelim = true, header_exists = false,
			colnames = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])))
		blast_out[file_name] = insertcols(blast_out[file_name], 1, :index => range(1, length(blast_out[file_name])))
	end
	return blast_out
end
logIO("blastIO function initiated! - ")

logIO("Initiating sseqid_merge function - ")
function sseqid_merge(n::Integer) #DONE
	global sseqid_parse = Dict{String, Any}()
	temp_counter = 0
	if length(keys(reference)) == 1 #FOR ONE REFERENCE AND MULTIPLE BLAST_OUTS
		for key in keys(reference)
			for n in 1:length(keys(blast_out))
				temp_dict = Dict{String, Vector}()
				key_bridge = SubString(key, 1, Integer(findfirst('_', key))-1)
				for sid in Set(select(blast_out["$key_bridge"*"$n"*"_out"], :sseqid))
					temp_table_s = filter(val -> (val.sseqid == sid), blast_out["$key_bridge"*"$n"*"_out"])
					for qid in Set(select(blast_out["$key_bridge"*"$n"*"_out"], :qseqid))
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
				end
				merge!(sseqid_parse, Dict(key => temp_dict))
			end
		end
	else #MULTIPLE REFERENCE AND MULTIPLE BLAST_OUTS
		for key in keys(reference)
			temp_dict = Dict{String, Vector}()
			key_bridge = SubString(key, 1, Integer(findfirst('_', key))-1)
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
	end
	return sseqid_parse
end
logIO("sseqid_merge function initiated! - ")

logIO("Initiating parseIO function - ")
using IterTools
function parseIO() #DONE
	for file in keys(sseqid_parse)
		logIO("Parsing sequences for hits found in $file - ")
		file_bridge = SubString(file, 1, Integer(findfirst('_', file))-1)
		io = open("$file_bridge"*"_parsed.fa", "w")
		for id in keys(sseqid_parse[file])
			if findall(r"^Chromosome", id) === true
				for (i, (sstart, send)) in enumerate(partition(sseqid_parse[file][id], 2, 2))
					if sstart > send
						println("yes")
						temp_substring = SubString(reference[file][id], send:sstart)
						write(io, ">$id"*"_$i\n") + write(io, "$temp_substring\n")
					else
						println("no")
						temp_substring = SubString(reference[file][id], sstart:send)
						write(io, ">$id"*"_$i\n") + write(io, "$temp_substring\n")
					end
				end
			else
				continue
			end
		end
		close(io)
	end
end
logIO("parseIO function initiated! - ")

referenceIO(ARGS[1]) #upload subject sequences
logIO("referenceIO completed! - ")
blastIO(ARGS[1]) #upload blast output
logIO("blastIO completed! - ")
sseqid_merge(2000) #merge hits within n basepairs and extract hits
logIO("Successfully merged hits within the specified number of basepairs! - ")
parseIO() #output extracted hits
close(log)