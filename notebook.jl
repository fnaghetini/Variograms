### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 1991a290-cfbf-11eb-07b6-7b3c8543dd28
begin
	# instantiate environment
    using Pkg; Pkg.activate(@__DIR__); Pkg.instantiate()
	
	# load packages used in this notebook
	using GeoStats, GeoStatsImages
	using CSV, DataFrames, Query
    using Statistics, StatsBase
	using Distributions, Random
	using PlutoUI
    using Plots, StatsPlots
	
	# default plot settings
	gr(format=:png)
end;

# ╔═╡ d4775d05-4943-4493-897e-4340f01475be
function sph2cart(azi)
	θ = deg2rad(azi)
	sin(θ), cos(θ)
end;

# ╔═╡ f8909bd5-9167-42ea-a302-a7a50bdc365c
html"""
<p style="background-color:lightgrey" xmlns:cc="http://creativecommons.org/ns#" xmlns:dct="http://purl.org/dc/terms/"><span property="dct:title">&nbsp&nbspVariograms</span> by <span property="cc:attributionName">Franco Naghetini</span> is licensed under <a href="http://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1"><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1"></a></p>
"""

# ╔═╡ 029c1951-054b-4f48-bc05-341250ce9f6a
md" # Variogramas"

# ╔═╡ 51107168-29ca-40b1-a658-9361199be3b1
md"""

## Função variograma e semivariograma

**Definição:** é uma **função matemática que mapeia o comportamento espacial de uma variável regionalizada** (e.g. Au, Cu, Co, Ag).

O variograma, quando existe, é único e válido para todo o domínio de estimativa.

A **função variograma pode ser anisotrópica**, sendo sensível à direção, mas não ao sentido (e.g. o variograma 000°/45° é igual ao variograma 180°/45°.

**Função variograma**:

```math
2\gamma(h) = \frac{1}{n} \sum_{i=1}^{n} [Z(x_i) - Z(x_i + h)]^2

```

**Função semivariograma**:

```math
\gamma(h) = \frac{1}{2n} \sum_{i=1}^{n} [Z(x_i) - Z(x_i + h)]^2

```

- **γ(h)**: valor do variograma $γ$ para uma distância $h$ entre dois pontos.

- **n**: número de pares de pontos.

- **Z(xᵢ)**: valor da variável $Z$ na posição $x_i$.

- **Z(xᵢ + h)**: valor da variável $Z$ na posição $x_i+h$.

Ao final da variografia, teremos em mãos um **modelo de variograma** representativo da continuidade espacial de uma variável e que será utilizado como **entrada no sistema linear de krigagem**.

"""

# ╔═╡ 0c00aee8-9db5-4fca-b92d-e19aa4fe5c1b
md"""
## Variogramas experimentais

A resolução da **função variograma** apresentada acima é representada graficamente e recebe o nome de **variograma experimental**.

Pode-se unir os pontos do variograma experimental para simplesmente facilitar a sua interpretação.

"""

# ╔═╡ b23b047e-1c02-40c5-ba88-825da85ba75c
md"""
Unir pontos do variograma: $(@bind join_points CheckBox())

"""

# ╔═╡ 4b12eecc-0645-4f46-b3be-8b8a095af599
begin
	# Amostras
	samples = geostatsimage("Gaussian30x10")

	# Variograma
	γ₁ = EmpiricalVariogram(samples, :Z, maxlag=60., nlags = 6)
	
	# Plotagem do variograma
	plot(γ₁, legend = false, ylims = (0,1.0), xlims = (0,60),
		 title = "Variograma Experimental", color = :orange,
	     line = join_points)
	
end

# ╔═╡ 5e623ea7-03f9-46a9-ba54-6d48d1a64057
md"""

Perceba que $h$ é a **distância cartesiana** entre duas amostras e $γ(h)$ seria a **distância estatística** entre essas amostras. Portanto, o **variograma experimental é a ferramenta que converte a distância geográfica em distância estatística**.

Note que a função variograma é uma **função discreta**, ou seja, os valores de $γ$ são calculados apenas para $h$ específicos.

Cada ponto no variograma experimental representa o valor $γ$ médio de um conjunto de pares de amostras separadas por uma distância $h$.

O número de pares de amostras é proporcional à altura das barras presentes no gráfico.

"""

# ╔═╡ 4b136ca1-f46f-43dc-9a1d-0659f1ef5e61
md""" ### Parâmetros para cálculo do variograma experimental

Para calcular variogramas experimentais, devemos definir alguns parâmetros:

- Direção
- Tamanho do passo
- Tolerância linear
- Tolerância angular
- Largura da banda

"""

# ╔═╡ c782a92c-cc4d-44bc-8521-2f70ad222bd5
md"""

#### Direção

Como visto anteriormente, os variogramas experimentais podem ser anisotrópicos, ou seja, variam de acordo com a direção. Nesse sentido, como os depósitos minerais são anisotrópicos, devemos escolher uma direção de cálculo para o variograma experimental.

- No contexto **2D**, informamos apenas o **azimute**.

- No contexto **3D**, informamos o **azimute** e o **mergulho**.

No exemplo 2D abaixo, percebemos que **quando variamos o azimute, o variograma experimental também sofre uma variação**:

"""

# ╔═╡ 43bc79ba-bb97-48bd-a8e4-c478bdc3a60b
md"""

Azimute: $(@bind azi Slider(0.:45.:135., default=0., show_value = true))°

"""

# ╔═╡ 3f39dcf2-055e-4aa8-8caa-09223175a5fa
begin
	walker_lake = CSV.File("data/walker_lake_proj.csv") |> DataFrame
	wl = walker_lake[:,[:X,:Y,:PB]]
	wl_georef = georef(wl, (:X,:Y))
	
	Random.seed!(1234)

    γ₂ = DirectionalVariogram(sph2cart(azi), wl_georef, :PB,
                              maxlag = 200, nlags = 7)
	
	plot(γ₂, legend = false, xlims = (0,200), ylims = (0,15),
		 title = "$(azi)°", color = :orange)
end

# ╔═╡ 5cb62b37-fe28-4816-b7ed-f5f40df255dc
md"""

#### Tamanho do passo

O **tamanho do passo (lag)** é a **distância média entre os furos na direção** em que o **variograma experimental** está sendo calculado.

Abaixo, tem-se um exemplo de cálculo de **variograma experimental na direção E-W para uma malha regular**. Adotou-se um **tamanho de passo** igual a **1 m**:


"""

# ╔═╡ f4e189ac-5d12-4de5-80e1-516103e5950f
md"""

Posição do vetor: $(@bind ix₁ Slider(1:1:4, default=1))
$(@bind iy₁ Slider(1:1:5, default=1))

Passo: $(@bind h₁ Slider(1:1:3, default=1, show_value = true)) m

"""

# ╔═╡ c3135efd-e69c-4270-8b45-b3f9f2dd586c
begin
	Random.seed!(42)
	
	values₁ = DataFrame(Au = rand(25))
	coords₁ = PointSet([(i,j) for i in 1:5 for j in 1:5])
	samples₁ = georef(values₁, coords₁)

	plot(samples₁, xlims = (0,6), ylims = (0,6), title = "Au (g/t)",
		 xlabel = "X(m)", ylabel = "Y(m)")
	
	plot!([(ix₁,iy₁),(ix₁+h₁,iy₁)], arrow = true, color = :red, lw = 2)
	
end

# ╔═╡ 3d25e2bc-9a6d-4712-92ed-de31dbdea3f2
md"""

#### Tolerância linear

Sabe-se que, na maioria dos casos, as malhas de sondagem são irregulares. Nesse caso, poucos pares de pontos serão buscados, uma vez que as amostras não se encontram equidistantes entre si.

"""

# ╔═╡ 874544a1-07af-4509-a34d-68d77558aaae
md"""

Posição do vetor: $(@bind ix₂ Slider(0.0:0.01:1.0, default=0.01))
$(@bind iy₂ Slider(0.0:0.01:1.0, default=0.17))

Tamanho do passo: $(@bind h₂ Slider(0.05:0.05:0.5, default=0.3, show_value = true)) m

"""

# ╔═╡ 2965ea0b-9b5e-4460-a735-e555733b2d83
begin
	Random.seed!(42)
	
	table₂ = georef(values₁, PointSet(rand(2,25)))
	
	plot(table₂, xlims = (-0.2,1.2), ylims = (-0.2,1.2), title = "Au (g/t)",
		 xlabel = "X(m)", ylabel = "Y(m)")
	
	plot!([(ix₂,iy₂),(ix₂+h₂,iy₂)], arrow = true, color = :red, lw = 2)
end

# ╔═╡ 7c00b7a2-5c09-46f5-ba8d-03786fd606b8
md"""

Uma alternativa comumente utilizada é a definição de uma **tolerância de passo** para que mais amostras sejam buscadas.

A **tolerância do passo (lagtol)** é definida como **metade do tamanho do passo**:

```math
lagtol = \frac{lag}{2} 

```

Essa abordagem permite que:

- Não haja *overlap* de amostras. As amostras são utilizadas apenas em um passo.

- Não haja *gap* de amostras. As amostras serão sempre utilizadas em algum dos passos.

"""

# ╔═╡ f1f163f7-eabc-4808-82e4-98ecfeddc730
md"""

O exemplo abaixo ilustra um perfil de uma **malha amostral irregular**. Note que, caso a tolerância de passo não fosse adotada, nenhum par de pontos seria buscado para o cálculo do variograma.

"""

# ╔═╡ 841ffdd2-16b4-4c19-8f03-70942a4ebb2e
md"""

Posição do vetor: $(@bind ix₃ Slider(1.:0.1:2.8, default=1))

"""

# ╔═╡ f738e56e-7725-4d52-a700-960ce372f025
begin
	Random.seed!(42)
	coords₃ = [(1.,1.),(1.6,1.),(1.9,1.),(2.2,1.),(2.8,1.),(3.6,1.),(3.8,1.)]
	values₃ = DataFrame(Au = rand(7))
	samples₃ = georef(values₃, coords₃)
	
	plot(samples₃, xlims = (0.,4.5), ylims = (0,2), title = "Au (g/t)",
		 xlabel = "X(m)", ylabel = "Y(m)")
	
	plot!([(ix₃,1.),(ix₃+1,1.)], arrow = true, color = :red, lw = 2)
	
	vline!([ix₃ + 0.5], color = :gray, ls = :dash)
	annotate!(ix₃ + 0.5, 2.1, text("lag - ½ lag", 7, :gray))
	
	vline!([ix₃+1], color = :red)
	annotate!(ix₃+1, 2.1, text("lag", 7, :red))
	
	vline!([ix₃ + 1.5], color = :gray, ls = :dash)
	annotate!(ix₃ + 1.5, 2.1, text("lag + ½ lag", 7, :gray))
	
end

# ╔═╡ ace40206-7ce6-4a64-b1ae-bd19d295158e
md"""

#### Tolerância angular

Além da tolerância linear (i.e. tolerância de passo), deve-se também definir uma **tolerância angular** que, por sua vez, é dividida em dois parâmetros:

- **Tolerância de azimute**

- **Tolerância de mergulho**

A **tolerância angular (angtol)** é definida como metade do incremento angular (*anginc*):

```math
angtol = \frac{anginc}{2} 

```

A convenção acima permite que:

- Não haja *overlap* de amostras.

- Não haja *gap* de amostras.


> A fórmula acima é válida tanto para a tolerância de azimute quanto para a tolerância de mergulho.

A figura abaixo ilustra um exemplo de tolerâncias angulares (azimute e mergulho) para um incremento angular de 45°:

"""

# ╔═╡ 728f75bd-0fc5-43c6-9551-4304925aa97b
md"""

Azimute: $(@bind dip_dir Slider(0:45:135, default=0, show_value = true))°

Mergulho: $(@bind dip Slider(0:45:135, default=0, show_value = true))°

"""

# ╔═╡ 4c5b95f2-5ad6-4f18-9cc0-9bd96eb3bf29
begin
	
	dip_dir_tol = plot([sph2cart(dip_dir+180),sph2cart(dip_dir)],
					   color = :red, ticks = false, xlims = (-1,1),
					   ylims = (-1,1), grid = false, lw = 3,
					   arrow = true, axis = false, legend = false,
					   title = "Tolerância de Azimute", size = (300,300))
	
	plot!([sph2cart(dip_dir+180-22.5),sph2cart(dip_dir-22.5)], color = :gray,
		  ls = :dash)
	
	plot!([sph2cart(dip_dir+180+22.5),sph2cart(dip_dir + 22.5)], color = :gray,
		  ls = :dash)
	
	vline!([0], color = :black)
	
	hline!([0], color = :black)
end;

# ╔═╡ 7fa3052f-52c8-48b5-ab3a-8401a6d8f93a
begin
	
	dip_tol = plot([sph2cart(dip+270),sph2cart(dip+90)], color = :red,
				   ticks = false, xlims = (-1,1), ylims = (-1,1),
				   grid = false, lw = 3, arrow = true, axis = false,
				   legend = false, title = "Tolerância de Mergulho",
				   size = (300,300))
	
	plot!([sph2cart(dip+270-22.5),sph2cart(dip+90-22.5)], color = :gray,
		  ls = :dash)
	
	plot!([sph2cart(dip+270+22.5),sph2cart(dip+90+22.5)], color = :gray,
		  ls = :dash)
	
	hline!([0], color = :black)
	
end;

# ╔═╡ 9709372d-3d5f-4bff-8ca1-adbb4dbeda23
plot(dip_dir_tol, dip_tol, layout = (1,2), size = (600,300))

# ╔═╡ 5e555810-f34d-402c-ac0a-17a423f420bc

md"""
#### Largura da banda

"""

# ╔═╡ 53979b80-406d-4eb9-aef8-5d3b626a5555
md"""

#### Resumo

- Variogramas **isotrópicos** são denominados **omnidirecionais**, ao passo que os variogramas **anisotrópicos** são chamados de **direcionais**.

- Como o variograma pode assumir anisotropia, devemos **definir uma ou mais direções para o cálculo do(s) variograma(s) experimental(is)**.

- Os **parâmetros de tolerância** (i.e. linear e angular) só são **necessários** em um contexto de malha **amostral irregular**.

- Em um **contexto 2D** de malha irregular, **apenas** as **tolerâncias de passo e de azimute** são necessárias.

- A **largura de banda** é um parâmetro **facultativo e restritivo**.

"""

# ╔═╡ 7723e35c-676b-465e-9a34-597745c3d284
md"""
Modelo Teórico: $(@bind model Select(["Gaussiano","Esférico","Pentaesférico","Exponencial"],
		default = "Esférico"))
"""

# ╔═╡ 6d0f5d99-f7e2-4f53-b835-c3b345613e4a
md"""

Ef. Pepita (c₀): $(@bind c₀ Slider(0.0:0.1:0.5, default=0.0, show_value=true))

Patamar (c₀+c₁): $(@bind cₜ Slider(0.5:0.1:1.0, default=1.0, show_value=true))

Alcance (a): $(@bind a Slider(3.0:2.0:10.0, default=5.0, show_value=true)) m

"""

# ╔═╡ 341ec3f6-c871-431f-8ffa-85f4c43ae138
if model == "Gaussiano"
	γ = GaussianVariogram(nugget = Float64(c₀),
						  sill = Float64(cₜ),
		 				  range = Float64(a))

elseif model == "Esférico"
	γ = SphericalVariogram(nugget = Float64(c₀),
						   sill = Float64(cₜ),
						   range = Float64(a))
	
elseif model == "Pentaesférico"
	γ = PentasphericalVariogram(nugget = Float64(c₀),
						   		sill = Float64(cₜ),
						   		range = Float64(a))

else
	γ = ExponentialVariogram(nugget = Float64(c₀),
						   sill = Float64(cₜ),
						   range = Float64(a))
end;

# ╔═╡ 61b8631b-8295-4dea-a5dd-189bf578bc8c
begin
	plot(γ, color = :black, lw = 2, label = model,
		 legend = :topleft, ylims = (0.,1.5), xlims = (0.,12.))
	
	hline!([c₀], color = :orange, ls = :dash, label = "Ef. Pepita")
	annotate!(11,c₀+0.05,text("Ef. Pepita",10,:orange))
	
	hline!([cₜ], color = :brown, ls = :dash, label = "Patamar")
	annotate!(11,cₜ+0.05,text("Patamar",10,:brown))
	
	vline!([a], color = :purple1, ls = :dash, label = "Alcance")
	annotate!(a,-0.05,text("Alcance",10,:purple1))
end

# ╔═╡ Cell order:
# ╟─1991a290-cfbf-11eb-07b6-7b3c8543dd28
# ╟─d4775d05-4943-4493-897e-4340f01475be
# ╟─f8909bd5-9167-42ea-a302-a7a50bdc365c
# ╟─029c1951-054b-4f48-bc05-341250ce9f6a
# ╟─51107168-29ca-40b1-a658-9361199be3b1
# ╟─0c00aee8-9db5-4fca-b92d-e19aa4fe5c1b
# ╟─b23b047e-1c02-40c5-ba88-825da85ba75c
# ╟─4b12eecc-0645-4f46-b3be-8b8a095af599
# ╟─5e623ea7-03f9-46a9-ba54-6d48d1a64057
# ╟─4b136ca1-f46f-43dc-9a1d-0659f1ef5e61
# ╟─c782a92c-cc4d-44bc-8521-2f70ad222bd5
# ╟─3f39dcf2-055e-4aa8-8caa-09223175a5fa
# ╟─43bc79ba-bb97-48bd-a8e4-c478bdc3a60b
# ╟─5cb62b37-fe28-4816-b7ed-f5f40df255dc
# ╟─c3135efd-e69c-4270-8b45-b3f9f2dd586c
# ╟─f4e189ac-5d12-4de5-80e1-516103e5950f
# ╟─3d25e2bc-9a6d-4712-92ed-de31dbdea3f2
# ╟─2965ea0b-9b5e-4460-a735-e555733b2d83
# ╟─874544a1-07af-4509-a34d-68d77558aaae
# ╟─7c00b7a2-5c09-46f5-ba8d-03786fd606b8
# ╟─f1f163f7-eabc-4808-82e4-98ecfeddc730
# ╟─f738e56e-7725-4d52-a700-960ce372f025
# ╟─841ffdd2-16b4-4c19-8f03-70942a4ebb2e
# ╟─ace40206-7ce6-4a64-b1ae-bd19d295158e
# ╟─4c5b95f2-5ad6-4f18-9cc0-9bd96eb3bf29
# ╟─7fa3052f-52c8-48b5-ab3a-8401a6d8f93a
# ╟─9709372d-3d5f-4bff-8ca1-adbb4dbeda23
# ╟─728f75bd-0fc5-43c6-9551-4304925aa97b
# ╟─5e555810-f34d-402c-ac0a-17a423f420bc
# ╟─53979b80-406d-4eb9-aef8-5d3b626a5555
# ╟─341ec3f6-c871-431f-8ffa-85f4c43ae138
# ╟─61b8631b-8295-4dea-a5dd-189bf578bc8c
# ╟─7723e35c-676b-465e-9a34-597745c3d284
# ╟─6d0f5d99-f7e2-4f53-b835-c3b345613e4a
