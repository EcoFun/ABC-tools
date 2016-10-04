
function flatten{T}(a::Vector{Vector{T}})
	res = T[]
	for sub in a
		append!(res, sub)
	end
	res
end


function get_arg(args, i, fun = x->x)
	if length(args) < i+1
		error("expected argument after $(args[i])")
	end

	i += 1
	fun(args[i]), i
end

function parse_group_spec(str::String)
	# list of populations
	selected = Vector{Int}[]
	names = ASCIIString[]
	
	groups = split(str, ',')

	for g in groups
		if ismatch(r"^[0-9]+$", g)
			push!(selected, [int(g)])
		elseif ismatch(r"^[0-9]+-[0-9]+$", g)
			limits = split(g, '-')
			push!(selected, [int(limits[1]):int(limits[2])])
		else
			push!(names, g)
		end
	end

	return selected, names
end

