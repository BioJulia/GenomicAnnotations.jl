using Documenter, GenomicAnnotations

makedocs(sitename = "GenomicAnnotations.jl", authors = "Karl Dyrhage",
    pages = [
        "index.md",
        "I/O" => "io.md",
        "Accessing and modifying annotations" => "accessing.md",
        "Representing genomic loci" => "loci.md",
        "Filtering: the @genes macro" => "genes.md",
        "Examples" => "examples.md"
    ],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true")
)

deploydocs(
    repo = "github.com/BioJulia/GenomicAnnotations.jl.git",
)
