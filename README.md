# Variograms

Commonly, natural phenomena such as grade mineralization show spatial interdependence. In this sense, the so-called variograms can be used in order to identify the mathematical law that describes the spatial structure of a regionalized variable. This notebook aims to present an intuitive and interactive introduction to a key tool in geostatistics: the **Variogram**.

![image](https://user-images.githubusercontent.com/63740520/123526403-0915a880-d6ae-11eb-8a5d-f53d5d57c968.png)

## Instructions

1. Download and install [Julia ≥ v1.6.1](https://julialang.org/downloads/).

2. Open Julia REPL and install [Pluto](https://github.com/fonsp/Pluto.jl):
```julia
julia> using Pkg
julia> Pkg.add("Pluto")
```

3. Download and extract the repository folder:

![image](https://user-images.githubusercontent.com/63740520/123525964-01083980-d6ab-11eb-9c11-e02875ae1c71.png)

4. In Julia REPL, instantiate environment:
```julia
julia> using Pkg
julia> Pkg.activate("C:/downloaded/folder/local/path/.../Variograms-main")
julia> Pkg.instantiate()
```
**NOTE:** Replace all `\` with `/` (for Windows users)

5. Launch Pluto:
```julia
julia> using Pluto
julia> Pluto.run()
```

6. Copy `Variograms-main` local path, paste it on *"Open from file"* field and, finally, add `\variograms.jl` to path, so you can open `variograms.jl` notebook:

![image](https://user-images.githubusercontent.com/63740520/123525997-46c50200-d6ab-11eb-8752-19cb623badfb.png)

## References

Main references used (strongly recommended):

- [GeoStats Tutorials](https://github.com/JuliaEarth/GeoStatsTutorials)
- [An introduction to applied geostatistics](https://www.amazon.com.br/Introduction-Applied-Geostatistics-Edward-Isaaks/dp/0195050134)
