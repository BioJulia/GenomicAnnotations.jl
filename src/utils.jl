"""
	oneline(s)

Remove newlines, replacing them with spaces when it seems appropriate.
"""
function oneline(v::AbstractString)
	R = split(v, '\n')
	buf = IOBuffer()
	for (i, r) in enumerate(R)
		if i != lastindex(R) && (r[end] == '-' || occursin(' ', r))
			print(buf, r, " ")
		else
			print(buf, r)
		end
	end
	String(take!(buf))
end
oneline(v::Any) = v
