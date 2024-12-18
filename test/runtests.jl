using GenomicAnnotations
using BioSequences
using Test

@testset "GenomicAnnotations" begin
    @testset "GenBank parsing" begin
        chrs = collect(open(GenBank.Reader, "example.gbk"))
        @test length(chrs) == 2
        @test chrs[2].name == "plasmid1"
    end

    chr = collect(open(GenBank.Reader, "example.gbk"))[1]

    @testset "GFF parsing" begin
        open(GFF.Writer, "example.gff") do w
            write(w, chr)
        end
        gff = collect(open(GFF.Reader, "example.gff"))[1]
        @test begin
            gbkbuf = IOBuffer()
            gffbuf = IOBuffer()
            print(gbkbuf, chr.genes[2:4])
            print(gffbuf, gff.genes[2:4])
            String(take!(gbkbuf)) == String(take!(gffbuf))
        end
    end

    @testset "EMBL parsing" begin
        open(EMBL.Writer, "example.embl") do w
            write(w, chr)
        end
        embl = collect(open(EMBL.Reader, "example.embl"))[1]
        @test begin
            gbkbuf = IOBuffer()
            emblbuf = IOBuffer()
            print(gbkbuf, chr.genes[2:4])
            print(emblbuf, embl.genes[2:4])
            String(take!(gbkbuf)) == String(take!(emblbuf))
        end
    end

    @testset "GTF parsing" begin
        gff = collect(open(GFF.Reader, "example.gff"))[1]
        open(GTF.Writer, "example.gtf") do w
            write(w, gff)
        end
        gtf = readgtf("example.gtf")[1]
        @test begin
            gff_data = map(row -> map(x -> string(coalesce(x, "")), row), eachrow(gff.genedata))
            gtf_data = map(row -> map(x -> string(coalesce(x, "")), row), eachrow(gtf.genedata))
            gff_data == gtf_data
        end
    end

    @testset "Extended methods" begin
        @test length(chr.genes[1]) == length(sequence(chr.genes[1]))
        @test length(chr.sequence) == length(sequence(chr))
    end

    @testset "Gene properties" begin
        @test length(propertynames(chr.genes[1])) == 12
        @test chr.genes[2].locus_tag == "tag01"
        @test (chr.genes[2].locus_tag = "tag01") == "tag01"
        @test begin
            chr.genes[2].test = "New column"
            chr.genes[2].test == "New column" && all(ismissing, chr.genes[[1,3,4,5,6]].test)
        end
        @test get(chr.genes[1], :locus_tag, "") == ""
        @test get(chr.genes[2], :locus_tag, "") == "tag01"
        @test get(chr.genes[2], :locustag, "") == ""
        @test begin
            GenomicAnnotations.pushproperty!(chr.genes[2], :db_xref, "GI:123")
            chr.genes[2].db_xref == ["GI:1293614", "GI:123"]
        end
        @test GenomicAnnotations.vectorise(Union{Missing,Int}[1,1,1]) == [[1],[1],[1]]
    end

    @testset "Iteration" begin
        @test length([g.locus_tag for g in chr.genes]) == 8
        @test [loc for loc in Locus("complement(join(1..3,7..9))")] ==
              [loc for loc in Locus("complement(order(1..3,7..9))")] == [Complement(ClosedSpan(7:9)), Complement(ClosedSpan(1:3))]
        @test [loc for loc in Locus("complement(1..3)")] == [Complement(ClosedSpan(1:3))]
        @test [loc for loc in Locus("join(complement(1..3),complement(7..9))")] ==
              [loc for loc in Locus("order(complement(1..3),complement(7..9))")] == [Complement(ClosedSpan(1:3)), Complement(ClosedSpan(7:9))]
        @test [loc for loc in Locus("1^2")] == [Locus("1^2")]
        @test [loc for loc in Locus("1..3")] == [Locus("1..3")]
        @test [loc for loc in Locus("1")] == [Locus("1")]
    end

    @testset "eachposition" begin
        @test collect(eachposition(Locus("1..3"))) == [1, 2, 3]
        @test collect(eachposition(Locus("join(1..3,7..9)"))) == [1, 2, 3, 7, 8, 9]
        @test collect(eachposition(Locus("complement(join(1..3,7..9))"))) == reverse([1, 2, 3, 7, 8, 9])
        @test collect(eachposition(Locus("order(complement(1..3),7..9)"))) == [3, 2, 1, 7, 8, 9]
    end

    @testset "Adding/removing genes" begin
        addgene!(chr, :CDS, ClosedSpan(300:390), locus_tag = "tag04")
        @test chr.genes[end].locus_tag == "tag04"
        delete!(chr.genes[end])
        @test chr.genes[end-1].locus_tag == "reg01"
    end

    @testset "@genes" begin
        @test length(union(@genes(chr, CDS), @genes(chr, !CDS))) == length(chr.genes)
        @test length(intersect(@genes(chr, CDS), @genes(chr, !CDS))) == 0
        @test length(union(@genes(chr, gene), @genes(chr, !gene))) == length(chr.genes)
        @test length(intersect(@genes(chr, gene), @genes(chr, !gene))) == 0
        @test length(@genes(chr)) == length(chr.genes)
        @test @genes(chr, feature(gene) == $:CDS) == chr.genes[[2,4,6]]
        @test @genes(chr, feature(gene) == $:CDS) == @genes(chr, CDS)
        @test @genes(chr, iscomplement(gene)) == chr.genes[[5,6,7]]
        @test @genes(chr, feature(gene) == $:CDS, !iscomplement(gene)) == chr.genes[[2,4]]
        @test @genes(chr, length(gene) < 300)[1] == chr.genes[2]
        @test length(@genes(chr, get(gene, :locus_tag, "") == "")) == 3
        gene = chr.genes[3]
        @test @genes(chr, gene == $gene)[1] == chr.genes[3]
        d = Dict(:a => "tag01")
        @test @genes(chr, :locus_tag == d[$:a]) == @genes(chr, :locus_tag == "tag01")
    end

    @testset "Broadcast" begin
        # Broadcasted assignment on existing property
        chr.genes[3:4].gene .= "AXL2P"
        @test all(chr.genes[3:4].gene .== "AXL2P")
        # Broadcasted assignment on previously missing property
        chr.genes[3:4].newproperty .= true
        @test all(chr.genes[3:4].newproperty .== true)
        # Broadcasted assignment with @genes
        @genes(chr, gene).newproperty .= false
        @test all(chr.genes[[3,5]].newproperty .== false)
    end

    @testset "Locus" begin
        @test locus(chr.genes[2]) < locus(chr.genes[4])
        @test locus(chr.genes[2]) == locus(chr.genes[2])
        @test iscomplement(chr.genes[2]) == false
        @test iscomplement(chr.genes[5]) == true
        @test sequence(chr.sequence, Locus("1..3")) == dna"aaa"
        @test sequence(chr.sequence, Locus("<1..>3")) == dna"aaa"
        @test sequence(chr.sequence, Locus("3..8"), translate = true) == aa"M"
        @test sequence(chr.sequence, Locus("1^2")) == dna"aa"
        @test sequence(chr.sequence, Locus("1")) == dna"a"
        @test sequence(chr.sequence, Locus("join(1..3,11..13)")) == dna"aaaata"
        @test sequence(chr.sequence, Locus("order(1..3,11..13)")) == [dna"aaa", dna"ata"]
        @test sequence(chr.sequence, Locus("complement(join(1..3,11..13))")) == dna"tatttt"
        @test sequence(chr.sequence, Locus("complement(order(1..3,11..13))")) == [dna"tat", dna"ttt"]
        @test sequence(chr.sequence, Locus("join(complement(1..3),complement(11..13))")) == dna"ttttat"
        @test sequence(chr.sequence, Locus("order(complement(1..3),complement(11..13))")) == [dna"ttt", dna"tat"]

        @test iscompound(Locus("join(1..3,7..9)")) == true
        @test iscompound(Locus("order(1..3,7..9)")) == true
        @test iscompound(Locus("complement(join(1..3,7..9))")) == true
        @test iscompound(Locus("join(complement(1..3),complement(7..9))")) == true
        @test iscompound(Locus("1..10")) == false
        @test iscompound(Locus("<1..>10")) == false
        @test iscompound(Locus("1^2")) == false
    end

    seq = dna"atgtccatatacaacggtatctccacctcaggtttagatctcaacaacggaaccattgccgacatgagacagttaggtatcgtcgagagttacaagctaaaacgagcagtagtcagctctgcatctgaagccgctgaagttctactaagggtggataacatcatccgtgcaagaccaagaaccgccaatagacaacatatgtaa"
    @test sequence(chr.genes[2]) == seq
    @test length(chr.genes[2]) == length(seq)

    @testset "Empty Record" begin
        chr = GenBank.Record()
        @test chr.name == ""
        @test chr.sequence == dna""
        @test chr.header == ""
        @test chr.genes == Gene[]
        @test names(chr.genedata) == []
    end

    @testset "Utils" begin
        seq = dna"a"^100
        @test relative_position(seq, Locus("1..10"), :start) == 1/100
        @test relative_position(seq, Locus("1..10"), :stop) == 10/100
        @test relative_position(seq, Locus("1..10"), :middle) ≈ 5.5/100
        @test relative_position(seq, Locus("join(100,1..10)"), :middle) == 5/100
    end
end