# Variograms

Commonly, natural phenomena such as grade mineralization show spatial interdependence. In this sense, the so-called variograms can be used in order to identify the mathematical law that describes the spatial structure of a regionalized variable. This notebook aims to present an intuitive and interactive introduction to a key tool in geostatistics: the **Variogram**.

## Instructions

1. Download and install [Julia ≥ v1.6.1](https://julialang.org/downloads/).

2. Install [Pluto](https://github.com/fonsp/Pluto.jl):
```julia
julia> using Pkg
julia> Pkg.add("Pluto")
```

3. Download and extract the repository folder:

![image](https://user-images.githubusercontent.com/63740520/123523125-28094000-d698-11eb-9213-267048409f49.png)

4. Open Julia REPL and instantiate required packages:
```julia
julia> using Pkg
julia> Pkg.activate("C:/repository/folder/local/path/.../Variograms-main")
julia> Pkg.instantiate()
```

5. Launch Pluto:
```julia
julia> using Pluto
julia> Pluto.run()
```

6. Copy `Variograms-main` local path, paste it on *Enter path or URL...* field and, finally, add `\notebook.jl` to path, so you can open `notebook.jl` file:

![image](https://user-images.githubusercontent.com/63740520/123523343-5dfaf400-d699-11eb-879f-eae4dced91ce.png)
