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



function parsechromosome_gff(lines, G)
	chrs = Chromosome{G}[]
	chr = Chromosome{G}()
    genes = G[]
    iobuffer = IOBuffer()
	isheader = true
    isfooter = false
    qualifier = String("")
    content = String("")
    header = IOBuffer()
	currentchr = 0
	currentfasta = ""
	for (linecount, line) in enumerate(lines)
		# Store the header
		if line[1:2] == "##" && line != "##FASTA"
			isheader = true
			println(header, line)
		else
			isheader = false
		end
		if line == "##FASTA"
			isfooter = true
		elseif !isfooter && !isheader
			(seqid, source, feature, sstart, send, score, strand, phase, attributes) = split(line, '\t')
			if seqid ∉ [chr.name for chr in chrs]
				chr = Chromosome{G}()
				push!(chrs, chr)
				chr.name = seqid
				chr.genedata[!, :source] = Union{Missing, String}[]
				chr.genedata[!, :score] = Union{Missing, Float64}[]
				chr.genedata[!, :phase] = Union{Missing, Int}[]
			end
			locus = Locus(parse(Int, sstart):parse(Int, send), strand[1])
			addgene!(chr, Symbol(feature), locus)
			if source != "." && ismissing(chr.genes[end].source)
				chr.genes[end].source = source
			end
			if score != "." && ismissing(chr.genes[end].score)
				chr.genes[end].score = parse(Float64, score)
			end
			if phase != "."
				chr.genes[end].phase = parse(Int, phase)
			end
			if attributes != "."
				for attribute in split(attributes, ';')
					isempty(attribute) && continue
					qualifier, values = split(attribute, '=')
					values = split(values, ',')
					for value in values
						pushproperty!(chr.genes[end], Symbol(qualifier), String(value))
					end
				end
			end
		elseif isfooter
			if line[1] == '>'
				currentfasta = line[2:end]
				if currentchr > 0
					seq = String(take!(iobuffer))
					chrs[currentchr].sequence = LongDNASeq(seq)
				end
				currentchr = findfirst(chr -> chr.name == currentfasta, chrs)
			else
				if line[1] ∉ ['A', 'T', 'G', 'C']
					println(line)
				end
				print(iobuffer, line)
			end
		end
	end
	if !isempty(currentfasta)
		chrs[findfirst(chr -> chr.name == currentfasta, chrs)].sequence = LongDNASeq(String(take!(iobuffer)))
	end
	headerstring = String(take!(header))
	for chr in chrs
		chr.header = headerstring
	end
	return chrs
end


function readgff(filename::AbstractString, G::Type = Gene; gunzip = false)
	gz = filename[end-2:end] == "gz"
	if gz || gunzip
		GZip.open(f -> readgff(f, G), filename)
	else
		open(f -> readgff(f, G), filename)
	end
end

function readgff(input::IO, G::Type; gunzip = false)
	finished = false
	chrs = Chromosome{G}[]
	if gunzip
		lines = readlines(gzdopen(input))
	else
		lines = readlines(input)
	end
	currentline = 1
	parsechromosome_gff(lines, G)
end


function gffstring(gene::Gene)
	buf = IOBuffer()
	firstattribute = true
    for field in names(parent(gene).genedata)
        field in [:source, :score, :phase] && continue
        v = parent(gene).genedata[index(gene), field]
        if !ismissing(v)
			if firstattribute
				firstattribute = false
			else
				print(buf, ";")
			end
            if v isa AbstractVector
				print(buf, field, "=")
                for i in eachindex(v)
					print(buf, v[i])
					i == lastindex(v) ? print(buf, ";") : print(buf, ",")
                end
            else
				print(buf, field, "=", v)
            end
        end
    end
	join([parent(gene).name,
		get(gene, :source, "."),
		feature(gene),
		locus(gene).position.start,
		locus(gene).position.stop,
		get(gene, :score, "."),
		locus(gene).strand,
		get(gene, :phase, "."),
		String(take!(buf))], '\t')
end

printgff(filepath::AbstractString, chrs) = printgff(open(filepath, "w"), chrs)
printgff(io::IO, chr::Chromosome) = printgff(io, [chr])
function printgff(io::IO, chrs::AbstractVector{Chromosome{Gene}})
	iobuffer = IOBuffer()
	### Header
	print(iobuffer,  chrs[1].header)
	### Body
	for chr in chrs
		for gene in chr.genes
			println(iobuffer, gffstring(gene))
		end
	end
	### Footer
	if !all(isempty(chr.sequence) for chr in chrs)
		println(iobuffer, "##FASTA")
		for chr in chrs
			println(iobuffer, ">", chr.name)
			for s in Iterators.partition(chr.sequence, 80)
				println(iobuffer, join(s))
			end
		end
	end
	print(io, String(take!(iobuffer)))
end
