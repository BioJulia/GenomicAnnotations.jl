using Documenter, GenomicAnnotations

makedocs(sitename = "GenomicAnnotations.jl", authors = "Karl Dyrhage",
    pages = [
        "index.md",
        "I/O" => "io.md",
        "Accessing and modifying annotations" => "accessing.md",
        "Filtering: the @genes macro" => "genes.md"
    ],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true")
)

deploydocs(
    repo = "github.com/kdyrhage/GenomicAnnotations.jl.git",
)
