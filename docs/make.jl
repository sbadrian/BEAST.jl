using Documenter, BEAST

DocMeta.setdocmeta!(BEAST, :DocTestSetup, :(using BEAST); recursive=true)

makedocs(;
    clean=false,
    modules=[BEAST],
    authors="Kristof Cools <kristof.cools@ugent.be> and contributors",
    repo="https://github.com/krcools/BEAST.jl/tree/master",
    sitename="BEAST.jl Documentation",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true", canonical="https://krcools.github.io/BEAST.jl", assets=String[]
    ),
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Assembly" => "assemble.md",
        "Quadrature Strategies" => "quadstrat.md",
        #"Manual" => "manual.md",
        #"Geometry" => Any["Coordinate System" => "coordinateSys.md", "Sphere Dimensions" => "scatterer.md"],
        #"Excitations" => Any[
        #    "Plane Wave" => "planeWave.md",
        #    "Dipoles" => "dipoles.md",
        #    "Ring Currents" => "ringCurrents.md",
        #    "Spherical Modes" => "sphModes.md",
        #    "Uniform Static Field" => "uniformStatic.md",
        #],
        #"Further Details" => "details.md",
        #"Contributing" => "contributing.md",
        #"API Reference" => "apiref.md",
    ],
)

deploydocs(; repo="github.com/krcools/BEAST.jl.git")