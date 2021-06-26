# Variograms

Commonly, natural phenomena such as grade mineralization show spatial interdependence. In this sense, the so-called variograms can be used in order to identify the mathematical law that describes the spatial structure of a regionalized variable. This notebook aims to present an intuitive and interactive introduction to a key tool in geostatistics: the **Variogram**.

## Instructions

1. Download and install [Julia ≥ v1.6.1](https://julialang.org/downloads/).

2. Open Julia REPL and install [Pluto](https://github.com/fonsp/Pluto.jl):
```julia
julia> using Pkg
julia> Pkg.add("Pluto")
```

3. Download and extract the repository folder:

![image](https://user-images.githubusercontent.com/63740520/123523125-28094000-d698-11eb-9213-267048409f49.png)

4. In Julia REPL, instantiate required packages:
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

6. Copy `Variograms-main` local path, paste it on *"Open from file"* field and, finally, add `\notebook.jl` to path, so you can open `notebook.jl` file:

![image](https://user-images.githubusercontent.com/63740520/123523580-cf877200-d69a-11eb-8f75-5f27d8c5dafa.png)


## References

Main references used (strongly recommended):

- [GeoStats Tutorials](https://github.com/JuliaEarth/GeoStatsTutorials)
- [An introduction to applied geostatistics (Isaaks & Srivastava (1990)](https://www.amazon.com.br/Introduction-Applied-Geostatistics-Edward-Isaaks/dp/0195050134)
