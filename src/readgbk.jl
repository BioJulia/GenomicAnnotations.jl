"""
Return the LOCUS entry of the header.
"""
function parseheader(header::String)
    lines = split(header, "\n")
    locus = split(lines[1], r" +")[2]
    return locus
end


"""
Parse lines encoding genomic position, returning the feature as a `Symbol`, and an instance of `Locus`.
"""
function parseposition(line::String)
    feature, posstring = split(strip(line), r" +")
    if occursin(r"(\.\.|\^)", posstring)
        position = UnitRange(parse.(Int, filter.(c->isnumeric(c), split(posstring, r"(\.\.|\^)(.*(\.\.|\^))?")))...)
        strand = occursin("complement", posstring) ? '-' : '+'
    else
        position = parse(Int, posstring)
        position = position:position
        strand = '.'
    end
    complete_left = !occursin('<', posstring)
    complete_right = !occursin('>', posstring)
    order = Vector{UnitRange{Int}}()
	join = occursin("join", posstring)
    if join || occursin("order", posstring)
        for m in eachmatch(r"\d+(\.\.|\^)\d+", posstring)
            r = Meta.parse.(split(m.match, r"(\.\.|\^)"))
            push!(order, r[1]:r[2])
        end
    end
    return Symbol(feature), Locus(position, strand, complete_left, complete_right, order, join)
end


"""
Parse footer (sequence) portion of a GenBank file, returning a `String`.
"""
function filterseq(io::IOBuffer)
    line = String(take!(io))
    String(replace(filter(x -> !(isspace(x) || isnumeric(x)), line), r"[^ATGCNatgcn]" => "n"))
end


"""
    parsechromosome(lines)

Parse and return one chromosome entry, and the line number that it ends at.
"""
function parsechromosome(lines, G::Type = Gene)
    genes = G[]
    iobuffer = IOBuffer()
	isheader = true
    isfooter = false
    spanning = false
    position_spanning = false
    qualifier = String("")
    content = String("")
    header = ""

    feature = :source
    locus = Locus()

    chromosome = Chromosome{G}()

    linecount = 0
    for line in lines
        linecount += 1

        # Catch cases where there's no header
        if linecount == 1 && occursin(r"gene", line)
            isheader = false
        end

        ### HEADER
		if isheader && occursin(r"FEATURES", line)
            chromosome.header = String(take!(iobuffer))
			isheader = false

        elseif isheader
            linecount == 1 ? print(iobuffer, line) : print(iobuffer, '\n', line)

        # Check if the footer has been reached
		elseif !isheader && !isfooter && (occursin(r"^BASE COUNT", line) || occursin(r"ORIGIN", line))
			# Stop parsing the file when the list of genes is over
            isfooter = true
			iobuffer = IOBuffer()

        ### BODY
		elseif !isheader && !isfooter
            if position_spanning && occursin(r"  /", line)
                position_spanning = false
                spanningline = filter(x -> x != ' ', String(take!(iobuffer)))
                try
                    feature, locus = parseposition(spanningline)
                catch
                    println(spanningline)
                    println(line)
                    @error "parseposition(spanningline) failed at line $linecount"
                end
            elseif position_spanning
                print(iobuffer, line)
            end
            if occursin(r"^ {5}\S", line)
                spanning = false
                try
                    feature, locus = parseposition(line)
                catch
                    println(line)
                    @error "parseposition(line) failed at line $linecount"
                end
                addgene!(chromosome, feature, locus)
            elseif !spanning && occursin(r"^ +/", line)
                if occursin(r"=", line)
                    if occursin("=\"", line)
                        (qualifier, content) = match(r"^ +/([^=]+)=\"?([^\"]*)\"?$", line).captures
                        content = String(content)
                    else
                        (qualifier, content) = match(r"^ +/(\S+)=(\S+)$", line).captures
                        try
                            tmpcontent = Meta.parse(content)
                            tmpcontent isa Expr && throw(Meta.ParseError)
							content = tmpcontent
                        catch
                            content = Symbol(content)
                        end
                    end

                    if occursin(r"\".*[^\"]$", line)
                        spanning = true
                    end

					isempty(names(chromosome.genedata)) ?
						(chromosome.genedata[!, Symbol(qualifier)] = Union{Missing, typeof(content)}[content]) :
	                    pushproperty!(chromosome.genes[end], Symbol(qualifier), content)

                else
                    # Qualifiers without a value assigned to them end up here
                    qualifier = split(line, '/')[2]
                    pushproperty!(chromosome.genes[end], Symbol(qualifier), true)
                end
            elseif spanning
                try
                    content = match(r" {21}([^\"]*)\"?$", line)[1]
                catch
                    @warn "Couldn't read content (line $linecount)"
                end
                if line[end] == '"'
                    spanning = false
                end
                if eltype(chromosome.genedata[!, Symbol(qualifier)]).b <: AbstractArray
                    i = index(chromosome.genes[end])
                    chromosome.genedata[!, Symbol(qualifier)][end][end] = Base.getproperty(chromosome.genes[end], Symbol(qualifier))[end] * "\n" * content
                else
                    Base.setproperty!(chromosome.genes[end], Symbol(qualifier), Base.getproperty(chromosome.genes[end], Symbol(qualifier)) * "\n" * content)
                end
            end

        ### FOOTER
        elseif isfooter
            if line == "//"
                break
            end
            if occursin(r"^ ", line)
                print(iobuffer, line)
            end
        end
    end
    chromosome.name = parseheader(chromosome.header)
    chromosome.sequence = LongDNASeq(filterseq(iobuffer))
    return linecount, chromosome
end


"""
    readgbk(input, G::Type = Gene; gunzip = false)

Parse GenBank-formatted file, returning a `Vector{Chromosome}`. `input` can be a file path or an `IOStream`. File names ending in ".gz" are assumed to be gzipped and are decompressed. Setting `gunzip` to `true` forces this behaviour.
The type of `AbstractGene` to be used can be specified with `G`, though currently the only option is `Gene`.
"""
function readgbk(filename::AbstractString, G::Type = Gene; gunzip = false)
    gz = filename[end-2:end] == ".gz"
    if gz || gunzip
        GZip.open(f -> readgbk(f, G), filename)
    else
        open(f -> readgbk(f, G), filename)
    end
end

function readgbk(input::IO, G::Type = Gene; gunzip = false)
	finished = false
	chrs = Chromosome{G}[]
	if gunzip
		lines = readlines(gzdopen(input))
	else
		lines = readlines(input)
	end
	currentline = 1
	while !finished
		if currentline >= length(lines)
			break
		end
		i, chr = parsechromosome(lines[currentline:end], G)
		currentline += i
		push!(chrs, chr)
	end
	chrs
end


## Printing


"""
    printgbk([io], chr)

Print `chr` in GenBank format.
"""
function printgbk(chrs::AbstractVector{C}) where {C <: Chromosome}
    io = IOBuffer()
    printgbk(io, chrs)
end
function printgbk(io::IO, chrs::AbstractVector{C}) where {C <: Chromosome}
    for chr in chrs
        printgbk(io, chr)
    end
    return io
end
function printgbk(chr::C) where {C <: Chromosome}
    io = IOBuffer()
    printgbk(io, chr)
end
function printgbk(io::IO, chr::C) where {C <: Chromosome}
    println(io, chr.header)
    println(io, rpad("FEATURES", 21, ' '), "Location/Qualifiers")
    println(io, chr.genes)
    println(io, "ORIGIN")
    formatsequence(chr.sequence, io)
    println(io)
    println(io, "//")
    return io
end


function formatsequence(sequence, io = IOBuffer)
    p = length(string(length(sequence))) + 2
    if length(sequence) > 60
        intervals = [i:i+60 for i in range(1; step = 60, stop = length(sequence)-60)]
        for interval in intervals
            println(io, lpad(string(first(interval)), p, ' '), " ", sequence[interval[1:10]],
                " ", sequence[interval[11:20]], " ", sequence[interval[21:30]],
                " ", sequence[interval[31:40]], " ", sequence[interval[41:50]],
                " ", sequence[interval[51:60]])
        end
    else
        intervals = [1:1]
    end
    i = intervals[end].stop
    if i <= length(sequence)
        print(io, lpad(i, p, ' '), " ")
        j = 0
        while i+j <= length(sequence)
            print(io, sequence[i+j])
            (j+1) % 10 == 0 && print(io, " ")
            j += 1
        end
    end
    return io
end
