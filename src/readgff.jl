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

"""
    readgff(input, G::Type = Gene; gunzip = false)

Parse GFF3-formatted file, returning a `Vector{Chromosome}`. `input` can be a file path or an `IOStream`. File names ending in ".gz" are assumed to be gzipped and are decompressed. Setting `gunzip` to `true` forces this behaviour.
The type of `AbstractGene` to be used can be specified with `G`, though currently the only option is `Gene`.
"""
f
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


## Printing


"""
	_oneline(s)

Remove newlines, replacing them with spaces when it seems appropriate.
"""
function _oneline(v::AbstractString)
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
_oneline(v::Any) = v


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
					print(buf, _oneline(v[i]))
					i == lastindex(v) ? print(buf, ";") : print(buf, ",")
                end
            else
				print(buf, field, "=", _oneline(v))
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

"""
    printgff(io::IO, chr)
    printgff(path::AbstractString, chr)

Print `chr` in GFF3 format.
"""
printgff(filepath::AbstractString, chrs) = printgff(open(filepath, "w"), chrs)
printgff(io::IO, chr::Chromosome) = printgff(io, [chr])
function printgff(io::IO, chrs::AbstractVector{Chromosome{Gene}})
	iobuffer = IOBuffer()
	### Header
	if chrs[1].header[1] == '#'
		print(iobuffer, chrs[1].header)
	else
		println(iobuffer, "##gff-version 3")
	end
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
