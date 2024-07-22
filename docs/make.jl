using ChromaProcesses
using Documenter

DocMeta.setdocmeta!(ChromaProcesses, :DocTestSetup, :(using ChromaProcesses); recursive=true)

makedocs(;
    modules=[ChromaProcesses],
    authors="Jonathan Miller jonathan.miller@fieldofnodes.com",
    sitename="ChromaProcesses.jl",
    format=Documenter.HTML(;
        canonical="https://fieldofnodes.github.io/ChromaProcesses.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/fieldofnodes/ChromaProcesses.jl",
    devbranch="main",
)
