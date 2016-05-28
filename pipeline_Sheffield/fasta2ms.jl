#!/usr/bin/env julia0.3.11

include("util.jl")


type Chunk
	id::String
	data::Vector{Char}
end

function Chunk(str)
	Chunk(str, Char[])
end

function add!(chunk::Chunk, str::Vector{Char})
	append!(chunk.data, str)
end

# read the entire fasta file
function read_fasta(input)
	# 1 chunk = 1 individual
	chunks::Vector{Chunk} = Chunk[]

	for line in eachline(input)
		# strip removes whitespace from beginning and end
		line = strip(line)

		# ignore empty lines
		if isempty(line)
			continue
		end

		# header
		if line[1] == '>'
			push!(chunks, Chunk(line))
		# data
		else
			add!(chunks[end], collect(line))
		end
	end
	
	# last expression -> return value of function
	chunks
end

type Pop
	name
	range
end

# assign names to populations
function get_pop_list(chunks::Vector{Chunk})
	cur_name = ""
	pop_start = 0
	pops = Pop[]

	for i in 1:length(chunks)
		header = chunks[i].id

		id = split(header, '#')[1]
		name = split(id, ':')[2]
		if isempty(name)
			println("ERROR: no ind name in $header!")
			exit()
		end
		pop_name = split(name, '_')[1]

		if pop_name ≠ cur_name
			if cur_name ≠ ""
				push!(pops, Pop(cur_name, pop_start:i-1))
			end
			cur_name = pop_name
			pop_start = i
		end
		# last one
		if i == length(chunks)
			push!(pops, Pop(pop_name, pop_start:i))
		end
	end

	pops
end

# translate population names into indices
# CAUTION: no tests for redundance!
function find_sel_og(pop_inds, pop_names, og_inds, og_names, pop_list)
	for name in pop_names
		pop = findfirst(p -> p.name == name, pop_list)
		@assert pop!=0 "no population named $(name)"
		
		push!(pop_inds, [pop_list[pop].range])
	end

	for name in og_names
		pop = findfirst(p -> p.name == name, pop_list)
		push!(og_inds, [pop_list[pop].range])
	end

	println("found $(length(og_inds)) outgroup pops:")
	println("sizes: ", join(map(o->length(o), og_inds), ", "))
	println("found $(length(pop_inds)) pops:")
	println("sizes: ", join(map(o->length(o), pop_inds), ", "))

	pop_inds, og_inds
end

is_valid(c::Char) = c == 'A' || c == 'T' || c == 'G' || c == 'C'

# generate outgroup sequence from a list of og pops
function gen_outgroup(ogs::Vector{Vector{Chunk}})
	data = fill('.', length(ogs[1][1].data))

	n_err = 0

	for site in 1:length(data)
		alleles = Set()

		for pop in ogs
			pop_n = length(pop)

			pop_n_1 = 0; pop_n_2 = 0
			pop_a_1 = '.'; pop_a_2 = '.'
			for chunk in pop
				c = chunk.data[site]

				# probably N
				if ! is_valid(c)
					@assert c == 'N' || c == '-' "unexpected allele: $(c)!"

					# treat '-' as 'N'
					data[site] = 'N'
					break
				end

				if c == pop_a_1 || pop_a_1 == '.'
					pop_a_1 = c
					pop_n_1 += 1
					continue
				end

				if c == pop_a_2 || pop_a_2 == '.'
					pop_a_2 = c
					pop_n_2 += 1
					continue
				end

				error("not a diallelic site!")
			end	# chunk

			# N overrides, so no need to check the other pops
			data[site] != 'N' || break
			
			# at this point only sites with either one or two valid alleles remain

			o = 'e'

			# find 90% majority
			if (pop_n > 10 && pop_n_1/pop_n > 0.9) || pop_n_1 ≥ pop_n-1
				o = pop_a_1
			elseif (pop_n > 10 && pop_n_2/pop_n > 0.9) || pop_n_2 ≥ pop_n-1
				o = pop_a_2
			end

			# no allele yet
			if data[site] == '.'
				data[site] = o
			# non-unanimous vote => error
			elseif data[site] ≠ o
				data[site] = 'e'
			end

			push!(alleles, pop_a_1)
			# diallelic
			if pop_a_2 ≠ '.'
				push!(alleles, pop_a_2)
			end
		end	# pop

		@assert(is_valid(data[site]) || data[site] == 'N' || data[site] == 'e',
			"found $(data[site]) in outgroup!")

		# ignore Ns
		if data[site] == 'N'
			continue
		end
		
		# assign random allele for undetermined sites
		if data[site] == 'e'
			@assert(length(alleles) == 2 && count(is_valid, alleles) == 2,
				"unknown situation in outgroup: $(data[site]) - " * join(alleles, ','))

			n_err += 1
			# not the most elegant way to get a random element out of the set
			# but at least it's short
			if rand(0:1) == 0
				pop!(alleles)
			end
			data[site] = pop!(alleles)
		end

	end	# site

	Chunk("outgroup", data), n_err
end


# convert one locus of fasta data into ms-like 0/1 data 
# plus list of Ns and segregating sites
# MOST IMPORTANT FUNCTION!
# chunks=haplos; outgroup=outgroup_seq
function conv_ms(chunks::Vector{Chunk}, outgroup::Chunk)
	@assert length(chunks[1].data) == length(outgroup.data)	
	
	# list of 0s; we'll need int later anyway
	hasN = fill(0, length(chunks[1].data))
	# get Ns from outgroup
	for i in 1:length(hasN)
		if outgroup.data[i] == 'N'
			hasN[i] = 1
		end
	end

	# we use booleans here since we won't print this list
	ref = chunks[1].data
	is_valid_snp = fill(false, length(ref))
	res = Vector{Char}[]

	for chunk in chunks
		# make sure all chunks are the same length
		@assert	length(chunk.data) == length(hasN)	

		push!(res, Char[])

		# check all sites
		for j in 1:length(chunk.data)
			if chunk.data[j] == 'N' || chunk.data[j] == '-' # previously the algorithm considered '-' as a valid value for conversion!
				hasN[j] = 1
				push!(res[end], ' ')	# we have to push *something*
			else
				# check the the ancestral/derived state of the allele
				state = chunk.data[j] == outgroup.data[j]
				push!(res[end], state ? '0' : '1')
				# check SNP validity in term of polymorphism WITHIN the INGROUP
				# is_valid_snp[j] = !same || is_valid_snp[j]
				same = chunk.data[j] == ref[j]
				is_valid_snp[j] = is_valid_snp[j] || !same
			end
		end
	end

	# invalidate SNPs with Ns
	for i in 1:length(hasN)
		if hasN[i] == 1 
			is_valid_snp[i] = false
		end
	end

	hasN, res, is_valid_snp
end


function print_ms(locus, is_valid_snp, out, name)
	println(out, "// $name")

	positions = Float16[]
	nsegs = 0
	for i in 1:length(is_valid_snp)
		if is_valid_snp[i]
			nsegs += 1
			push!(positions, i/length(is_valid_snp))
		end
	end
	println(out, "segsites: $nsegs")


	if nsegs!=0
		println(out, "positions: " * join(positions, " "))
		for haplo in locus 
			for i in 1:length(haplo)
				# only print heterozygous sites
				# also takes care of Ns
				if is_valid_snp[i]
					print(out, haplo[i])
				end
			end
			println(out)
		end
	else
		print(out, "\n")
		println("$name has 0 variable site")
	end

	println(out)
end


function print_mask(mask, out)
	println(out, join(mask))
end

function print_spinput(spout, loci, pops, msname)
	println(spout, length(loci), "\t\t\t// #loci")
	println(spout, length(pops), "\t\t\t// #populations")
	for locus in loci
		for pop in pops
			println(spout, length(pop))
		end
		println(spout, length(locus[1]), "\t\t\t// length of locus [bp]")
	end

	println(spout, 1)
	println(spout, msname)
end


# supposedly Julia is faster if everything is in a function...
function run()

# *** commandline args

	pop_inds = Vector{Int}[]
	pop_names = ASCIIString[]
	og_inds = Vector{Int}[]
	og_names = ASCIIString[]
	fasta_files = ASCIIString[]

	# TODO: maybe have repeated -p -o add groups, put all inds within one
	# -p -o in one group

	# println(ARGS)
	i = 2
	while i ≤ length(ARGS)
		a = ARGS[i]

		if a == "-p"
			str, i = get_arg(ARGS, i)
			pop_inds, pop_names = parse_group_spec(str)
		elseif a == "-o"
			str, i = get_arg(ARGS, i)
			og_inds, og_names = parse_group_spec(str)
		else
			fasta_files = ARGS[i:end]
			i = length(ARGS)
		end

		i += 1
	end


# *** reading fasta

	# read the first file to get header information
	chunks = open(read_fasta, fasta_files[1])
	# get populations (only used to be able to select by pop name)
	pop_list = get_pop_list(chunks)
	selected, outgroup = find_sel_og(pop_inds, pop_names, og_inds, og_names, pop_list)
	# for the ms file pop structure doesn't matter, just use a list of haplos
	sel_haplos = flatten(selected)

	haplos = Vector{Char}[]
	loci = Vector{Vector{Char}}[]
	mask = Vector{Int}[]
	is_valid_snp = Vector{Bool}[]

	out_fname = ARGS[1]


	statsf = open(out_fname * ".stats", "w")
	println(statsf, "file\tNb.UnknownState.Sites\tNb.Clean.Sites")
	
	for name in fasta_files
		print(".")

		# read raw sequence data
		chunks = open(read_fasta, name)
		# generate list of outgroup pops (pop structure important!)
		ogs = map(og -> chunks[og], outgroup)

		# get ancestral genotype
		outgroup_seq, og_errors = gen_outgroup(ogs)
		og_valid = count(is_valid, outgroup_seq.data)
		@assert og_errors ≤ og_valid 
		println(statsf, "$name\t$og_errors\t$og_valid")
		
		haplos = chunks[sel_haplos]

		# convert to ms data
		hasN, conv, valid = conv_ms(haplos, outgroup_seq)
		
		push!(loci, conv)
		push!(mask, hasN)
		push!(is_valid_snp, valid)
	end

	close(statsf)

	println()


# *** writing ms

	println("writing ms file...")

	open(out_fname * ".ms", "w") do ms_out
		println(ms_out, "Observed data as ms format")
		println(ms_out, "1 2 3 # dummy seed numbers\n")
		for i in 1:length(loci)
			print_ms(loci[i], is_valid_snp[i], ms_out, ARGS[5+i])
		end
	end


# *** writing spinput file

	open(out_fname * ".spi", "w") do sp_out
		print_spinput(sp_out, loci, selected, out_fname * ".ms")
	end


# *** writing mask

# function print_mask(mask, out)
# 	println(out, join(mask))
# end
	println("writing mask file...")

	open(out_fname * ".mask", "w") do mask_out
		for locus in mask
			print_mask(locus, mask_out)
		end
	end
end

run()

