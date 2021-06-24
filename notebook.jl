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

**Definição:** é uma **função matemática que mapeia/descreve a continuidade espacial** (i.e. continuidade de teores) **de uma variável regionalizada** (e.g. Au, Cu, Co, Ag).

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

A resolução da **função variograma** apresentada acima é representada graficamente e recebe o nome de **variograma experimental** (Figura 1).

Pode-se unir os pontos do variograma experimental para simplesmente facilitar a sua interpretação.

"""

# ╔═╡ 4b12eecc-0645-4f46-b3be-8b8a095af599
begin
	# Amostras
	image = geostatsimage("Gaussian30x10")

	# Variograma
	γ₁ = EmpiricalVariogram(image, :Z, maxlag=60., nlags = 6)
	
end;

# ╔═╡ b23b047e-1c02-40c5-ba88-825da85ba75c
md"""
Unir pontos do variograma: $(@bind join_points CheckBox())

"""

# ╔═╡ 8cfef844-5e4d-44c8-817c-0021eecbcaa2
# Plotagem do variograma
plot(γ₁, legend = false, ylims = (0,1.0), xlims = (0,60),
	 title = "Variograma Experimental", color = :orange,
	 line = join_points)

# ╔═╡ 528f0bb5-4030-4006-a323-29f9cbc1efc0
html"""

<p align="center">
    <b>Figura 1</b>: Exemplo de variograma experimental.
</p>

"""

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

No exemplo 2D abaixo, percebemos que **quando variamos o azimute, o variograma experimental também sofre uma variação**.

Variogramas experimentais que assumem anisotropia são denominados **variogramas direcionais** (Figura 2).

"""

# ╔═╡ 43bc79ba-bb97-48bd-a8e4-c478bdc3a60b
md"""

Azimute: $(@bind azi Slider(0:45:135, default=0, show_value = true))°

"""

# ╔═╡ 3f39dcf2-055e-4aa8-8caa-09223175a5fa
begin
	walker_lake = CSV.File("data/walker_lake_proj.csv") |> DataFrame
	wl = walker_lake[:,[:X,:Y,:PB]]
	dropmissing!(wl)
	wl_georef = georef(wl, (:X,:Y))

    γ₂ = DirectionalVariogram(sph2cart(azi), wl_georef, :PB,
                              maxlag = 200, nlags = 7)
	
	plot(γ₂, legend = false, xlims = (0,200), ylims = (0,15),
		 title = "$(azi)°", color = :orange)
end

# ╔═╡ facdf937-4056-4699-b505-d9cada0c8ce3
html"""

<p align="center">
    <b>Figura 2</b>: Exemplo de variogramas direcionais.
</p>

"""

# ╔═╡ 5cb62b37-fe28-4816-b7ed-f5f40df255dc
md"""

#### Tamanho do passo

O **tamanho do passo (lag)** é a **distância média entre as amostras vizinhas** na direção em que o **variograma experimental** está sendo calculado.

Abaixo, tem-se um exemplo de cálculo de variograma experimental na direção E-W para uma malha regular (Figura 3). Adotou-se um tamanho de passo igual a 1 m.


"""

# ╔═╡ f4e189ac-5d12-4de5-80e1-516103e5950f
md"""

W-E: $(@bind ix₁ Slider(1:1:4, default=1))
N-S: $(@bind iy₁ Slider(1:1:5, default=1))

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

# ╔═╡ 8923e1c1-914d-47b5-a4b4-5f0c53c4e810
html"""

<p align="center">
    <b>Figura 3</b>: Exemplo de busca de pares de amostras em uma malha regular.
</p>

"""

# ╔═╡ 3d25e2bc-9a6d-4712-92ed-de31dbdea3f2
md"""

#### Tolerância linear

Sabe-se que, na maioria dos casos, as malhas de sondagem são irregulares. Nesse caso, poucos pares de pontos serão buscados, uma vez que as amostras não se encontram equidistantes entre si (Figura 4).

"""

# ╔═╡ 874544a1-07af-4509-a34d-68d77558aaae
md"""

W-E: $(@bind ix₂ Slider(0.0:0.001:1.0, default=0.015))
N-S: $(@bind iy₂ Slider(0.0:0.001:1.0, default=0.172))

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

# ╔═╡ 9d60b923-72e7-42c8-8ef3-d4a590e3f600
html"""

<p align="center">
    <b>Figura 4</b>: Exemplo de busca de pares de amostras em uma malha irregular.
</p>

"""

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

O exemplo abaixo ilustra amostras colineares de uma **malha amostral irregular** (Figura 5). Note que, caso a tolerância de passo não fosse adotada, apenas um par de pontos seria buscado para o cálculo do variograma.

"""

# ╔═╡ 8bd1b52b-b6a8-489e-ae74-be2931eef4ee
begin
	
	Random.seed!(42)
	coords₃ = [(1.,1.),(1.6,1.),(1.9,1.),(2.2,1.),(2.8,1.),(3.6,1.),(3.8,1.)]
	values₃ = DataFrame(Au = rand(7))
	samples₃ = georef(values₃, coords₃)
	
end;

# ╔═╡ 841ffdd2-16b4-4c19-8f03-70942a4ebb2e
md"""

Posição do vetor: $(@bind ix₃ Slider(1.:0.1:2.8, default=1))

Tolerância de passo: $(@bind lag_tol CheckBox())

"""

# ╔═╡ f738e56e-7725-4d52-a700-960ce372f025
begin
	
	if lag_tol
	
		plot(samples₃, xlims = (0.,4.5), ylims = (0,2), title = "Au (g/t)",
			 xlabel = "X(m)", ylabel = "Y(m)")

		plot!([(ix₃,1.),(ix₃+1,1.)], arrow = true, color = :red, lw = 2)

		vline!([ix₃ + 0.5], color = :gray, ls = :dash)
		annotate!(ix₃ + 0.5, 2.1, text("lag - ½ lag", 7, :gray))

		vline!([ix₃+1], color = :red)
		annotate!(ix₃+1, 2.1, text("lag", 7, :red))

		vline!([ix₃ + 1.5], color = :gray, ls = :dash)
		annotate!(ix₃ + 1.5, 2.1, text("lag + ½ lag", 7, :gray))
		
	else
		
		plot(samples₃, xlims = (0.,4.5), ylims = (0,2), title = "Au (g/t)",
			 xlabel = "X(m)", ylabel = "Y(m)")

		plot!([(ix₃,1.),(ix₃+1,1.)], arrow = true, color = :red, lw = 2)

		vline!([ix₃+1], color = :red)
		annotate!(ix₃+1, 2.1, text("lag", 7, :red))

	end
		
end

# ╔═╡ 650fc66a-3f8e-45d5-a898-5c783a8d12a1
html"""

<p align="center">
    <b>Figura 5</b>: Exemplo de amostras colineares não equidistantes.
</p>

"""

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

A Figura 6 ilustra um exemplo de tolerâncias angulares (azimute e mergulho) para um incremento angular de 45°.

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

# ╔═╡ b5bb3401-48d5-4d22-bbed-06427862062e
html"""

<p align="center">
    <b>Figura 6</b>: Tolerâncias angulares.
</p>

"""

# ╔═╡ 5e555810-f34d-402c-ac0a-17a423f420bc

md"""
#### Largura da banda

A **largura da banda** (largura máxima) é um **parâmetro de restrição facultativo** que pode ser utilizado em conjunto com a tolerância angular.

Esse parâmetro é definido pela distância entre a reta da direção de cálculo do variograma experimental e a linha de tolerância angular. A partir de uma distância igual à largura da banda a busca por amostras se paraleliza (Figura 7).

"""

# ╔═╡ 6433f0dc-04f8-450e-9a94-f8cfa8cda552
html"""

<p align="center">

    <img src="https://i.postimg.cc/7Zb0qPyc/Figure-01.png" style="height : 300px">

</p>

<p align="center">

    <b>Figura 7</b>: Largura da banda.

</p>

"""

# ╔═╡ 6f59e6db-0632-483a-89be-6c82dd188d60
md"""

### Malha regular x malha irregular

Para **malhas amostrais regulares**, os parâmetros de variograma experimental a serem definidos são:

- Número de direções de cálculo

    - Exemplo: 4 direções

- Direções de cálculo

    - Exemplo: 000°, 045°, 090°, 135°

- Tamanho do passo

    - Exemplo: 25 m

- Número de passos

    - Exemplo: 5 passos (25 m, 50 m, 75 m, 100 m, 125 m)

No caso de **malhas amostrais irregulares**, os parâmetros de variograma experimental a serem definidos são:

- Número de direções de cálculo

    - Exemplo: 4 direções

- Direções de cálculo e tolerância angular

    - Exemplo: 000° ± 22.5°, 045° ± 22.5°, 090° ± 22.5°, 135° ± 22.5°

- Tamanho do passo e tolerância do passo

    - Exemplo: 25 ± 12.5 m

- Número de passos

    - Exemplo: 5 passos (25 m, 50 m, 75 m, 100 m, 125 m)

- Largura da banda

    - Exemplo: 30 m

"""

# ╔═╡ 53979b80-406d-4eb9-aef8-5d3b626a5555
md"""

### Resumo

- Variogramas **isotrópicos** são denominados **omnidirecionais**, ao passo que os variogramas **anisotrópicos** são chamados de **direcionais**.

- Como o variograma pode assumir anisotropia, devemos **definir uma ou mais direções para o cálculo do(s) variograma(s) experimental(is)**.

- Os **parâmetros de tolerância** (i.e. linear e angular) só são **necessários** em um contexto de malha **amostral irregular**.

- Em um **contexto 2D de malha irregular**, apenas as **tolerâncias de passo e de azimute** são necessárias.

- A **largura de banda** é um parâmetro **facultativo e restritivo**.

"""

# ╔═╡ e80e5e16-59fb-4ec0-a9f0-6b8b97bc8d36
md"""

## Modelos teóricos de variograma

A partir dos variogramas experimentais só é possível obter valores médios de variograma (γ) para distâncias iguais a múltiplos do tamanho de passo (h) escolhido.

Portanto, é necessário o ajuste de um **modelo matemático contínuo**, de modo que saberemos o valor do variograma (γ) para qualquer distância entre pares de amostras (h).

O procedimento de se ajustar um modelo teórico contínuo ao variograma experimental é denominado **modelagem do variograma** (Figura 8).

"""

# ╔═╡ 9891913d-e735-4ec8-b09c-49b51417f18d
varmod₁ = fit(GaussianVariogram, γ₁);

# ╔═╡ 52c19e93-8534-4b59-a164-3f12d23d5440
md"""
Ajuste do modelo: $(@bind fit_model CheckBox())

"""

# ╔═╡ 7a8a49a2-007e-409a-9a45-c76f717f75f8
begin
	
	if fit_model
		plot(γ₁, color = :orange, line = false)
	
		plot!(varmod₁, ylims = (0,1.0), xlims = (0,50), label = "Ajuste Teórico",
		  	  legend = false)
	else
		plot(γ₁, color = :orange, line = false, ylims = (0,1.0), xlims = (0,50),
			 legend = false)
	end
	
end	

# ╔═╡ f92eedbf-097d-45f6-a550-ccc8c2f9841b
html"""

<p align="center">
    <b>Figura 8</b>: Exemplo de ajuste de um modelo teórico a um variograma experimental.
</p>

"""

# ╔═╡ ea67e941-496f-45a3-8b0f-058d573291d8
md"""

### Elementos do variograma

Para se modelar o variograma, deve-se configurar quatro elementos:

- Alcance (*Range*)

- Efeito Pepita (*Nugget Effect*)

- Patamar (*Sill*)

- Tipo de Modelo (*Model*)

"""

# ╔═╡ 03a46b21-c522-4036-a710-bd6ce0a26a1b
md"""

##### Alcance (a)

O **alcance (range)** consiste na distância máxima até onde se consegue estabelecer alguma interdependência espacial entre pares de amostras.

Em outras palavras, o alcance define até qual distância $h$ existe correlação espacial entre pares de amostras. Portanto, **para distâncias $h$ superiores ao alcance, não há interdependência espacial entre amostras**.

"""

# ╔═╡ 510759f5-3838-4db7-b683-e39677a2551b
md"""

##### Efeito Pepita (C₀)

O **efeito pepita (nugget effect)** é definido como a **descontinuidade próxima a origem do variograma**. É o valor de $γ(h)$ quando $h$ **tende** a zero. Note que:

$γ(h) = c₀, h → 0$
$γ(h) = 0, h = 0$

"""

# ╔═╡ c75b37b1-b13a-4d3b-949c-639e4e5dc01b
md"""

##### Patamar (C + C₀)

O **patamar (sill)** é definido como o máximo valor de γ(h) que as amostras podem apresentar.

Na prática, uma abordagem muito utilizada é considerar o patamar como o valor da **variância à priori** (i.e. variância amostral) da variável de interesse.

"""

# ╔═╡ 83593f8e-8dd2-40b1-903b-8712bb9eb048
md"""

##### Tipo de modelo

Nem todas as funções matemáticas podem ser utilizadas para ajustar variogramas experimentais. Para ser utilizada, a função determinística deve ser:

- Contínua

- Monotônica crescente

Existe um pouco mais de uma dezena de funções que podem ser utilizadas como ajustes teóricos. Dessas funções, há três que ajustam a grande maioria dos fenômenos naturais:

- Modelo Gaussiano

- Modelo Esférico

- Modelo Exponencial

O **modelo teórico controla a variabilidade a pequenas distâncias**, ou seja, o comportamento próximo à origem do variograma.


"""

# ╔═╡ a6802bda-7b7a-4d98-bb08-bcbe8a990e01
md"""

###### Modelo Gaussiano

- Apresenta **comportamento próximo à origem parabólico**.

- Esse tipo de modelo teórico é normalmente utilizado para ajustar variogramas experimentais de fenômenos de **baixa heterogeneidade**.

- Sua equação é descrita como:

``` math

γ(h) = C₀ + C \left[ 1 - exp \left[- \left(\frac{h}{a} \right)^2 \right]  \right] 

```

"""


# ╔═╡ bd8328a0-9cd3-4eec-8c68-f6f15f296f4b
md"""

###### Modelo Esférico

- Apresenta **comportamento próximo à origem linear**.

- Esse tipo de modelo teórico é normalmente utilizado para ajustar variogramas experimentais de fenômenos de **intermediária heterogeneidade**.

- É o **modelo mais utilizado na indústria da mineração**.

- Sua equação é descrita como:

``` math
γ(h) = C₀ + C \left[\frac{3h}{2a} - \frac{1}{2}\left(\frac{h}{a}\right)^3 \right], ∀ h < a)
```

``` math
γ(h) = C₀ + C, ∀ h ≥ a
```

"""

# ╔═╡ ed5011d1-8d0a-4ace-8baa-0af219b3e3e2
md"""

###### Modelo Exponencial

- Apresenta **comportamento próximo à origem linear**. Entretanto, a inclinação desse modelo nas proximidades da origem é maior do que a inclinação do Modelo Esférico.

- Esse tipo de modelo teórico é normalmente utilizado para ajustar variogramas experimentais de fenômenos de **elevada heterogeneidade**.

- Sua equação é descrita como:

``` math
γ(h) = C₀ + C \left[1 - exp \left[-\left(\frac{h}{a} \right) \right] \right]
```

"""

# ╔═╡ d1e9297b-04c0-47ad-91e8-77648b205d10
md"""

Termos das equações:

- **γ(h)** = variograma
- **C₀** = efeito pepita
- **h** = passo
- **C** = contribuição ao patamar
- **a** = alcance

"""

# ╔═╡ e7157b79-368d-4ce1-97d3-22110e3359da
md"""

A Figura 9 ilustra graficamente os quatro elementos que constituem o variograma. 

"""

# ╔═╡ 6d0f5d99-f7e2-4f53-b835-c3b345613e4a
md"""

Modelo Teórico: $(@bind model Select(["Gaussiano","Esférico","Pentaesférico","Exponencial"],
		default = "Gaussiano"))

Ef. Pepita (C₀): $(@bind c₀ Slider(0.0:0.1:0.5, default=0.0, show_value=true))

Patamar (C₀+C): $(@bind cₜ Slider(0.5:0.1:1.0, default=1.0, show_value=true))

Alcance (a): $(@bind a Slider(5.0:10:45.0, default=25.0, show_value=true)) m

"""

# ╔═╡ 341ec3f6-c871-431f-8ffa-85f4c43ae138
if model == "Gaussiano"
	γ₃ = GaussianVariogram(nugget = Float64(c₀),
						  sill = Float64(cₜ),
		 				  range = Float64(a))

elseif model == "Esférico"
	γ₃ = SphericalVariogram(nugget = Float64(c₀),
						   sill = Float64(cₜ),
						   range = Float64(a))
	
elseif model == "Pentaesférico"
	γ₃ = PentasphericalVariogram(nugget = Float64(c₀),
						   		sill = Float64(cₜ),
						   		range = Float64(a))

else
	γ₃ = ExponentialVariogram(nugget = Float64(c₀),
						   sill = Float64(cₜ),
						   range = Float64(a))
end;

# ╔═╡ 61b8631b-8295-4dea-a5dd-189bf578bc8c
begin
	plot(γ₃, color = :black, lw = 2, label = model,
		 legend = :topleft, ylims = (0.,1.5), xlims = (0.,65.))
	
	hline!([c₀], color = :red, ls = :dash, label = "Ef. Pepita")
	annotate!(55,c₀+0.05,text("Ef. Pepita",10,:red))
	
	hline!([cₜ], color = :green, ls = :dash, label = "Patamar")
	annotate!(55,cₜ+0.05,text("Patamar",10,:green))
	
	vline!([a], color = :blue, ls = :dash, label = "Alcance")
	annotate!(a,-0.05,text("Alcance",10,:blue))
end

# ╔═╡ ab5e6c19-789b-4944-ba8e-f983a9a2652c
html"""

<p align="center">
    <b>Figura 9</b>: Elementos do variograma.
</p>

"""

# ╔═╡ 8b4ee7b2-2a01-44ae-8539-27f1815fe634
md"""

## Tipos de anisotropia

Na geoestatística, a **anisotropia** existe quando um ou mais elementos do variograma variam com a mudança da direção. Existem três tipos (Figura 10):

- **Anisotropia Zonal**: patamar varia de acordo com a mudança de direção.

- **Anisotropia Geométrica**: alcance varia de acordo com a mudança de direção.

- **Anisotropia mista**: patamar e alcance variam de acordo com a mudança de direção.

> Embora existam três tipos de anisotropia, é comum considerar apenas a **anisotropia geométrica** para a modelagem do variograma.

> **Não existe anisotropia de efeito pepita**, uma vez que esse elemento é, por definição, isotópico.

"""

# ╔═╡ 83d6c4fe-bcd6-4d7f-93ef-2e093b1284fa
md"""

Tipo de anisotropia: $(@bind aniso Select(["Zonal","Geométrica","Mista"],
										  default = "Geométrica"))

"""

# ╔═╡ 187c01ca-053e-4994-a748-cf9b16683a50
if aniso == "Zonal"
	γ_aniso₁ = SphericalVariogram(nugget = 0.1, range = 50.0, sill = 1.0)
	γ_aniso₂ = SphericalVariogram(nugget = 0.1, range = 50.0, sill = 1.5)
	
elseif aniso == "Geométrica"
	γ_aniso₁ = SphericalVariogram(nugget = 0.1, range = 50.0, sill = 1.0)
	γ_aniso₂ = SphericalVariogram(nugget = 0.1, range = 30.0, sill = 1.0)
	
else
	γ_aniso₁ = SphericalVariogram(nugget = 0.1, range = 50.0, sill = 1.5)
	γ_aniso₂ = SphericalVariogram(nugget = 0.1, range = 30.0, sill = 1.0)
end;

# ╔═╡ b9634a1e-f225-4986-867f-fd36f56882df
begin
	
	# Variograma 1
	plot(γ_aniso₁, color = :red, lw = 2, legend = false,
		 title = "Anisotropia $aniso")
	
	# Variograma 2
	plot!(γ_aniso₂, color = :blue, lw = 2, xlims = (0,80), ylims = (0,2))
	
end

# ╔═╡ ee09dcab-2298-444c-ad9f-f79268c9056c
html"""

<p align="center">
    <b>Figura 10</b>: Tipos de anisotropia.
</p>

"""

# ╔═╡ 0f28a997-4945-47fe-83b9-058726bc8041
md"""

## Estruturas aninhadas

A **estrutura do variograma** é a porção da equação do ajuste teórico em que o valor de $C$ cresce com o aumento da distância $h$.

``` math
γ(h) = C_0 +
\underbrace{C \left[\frac{3h}{2a} - \frac{1}{2}\left(\frac{h}{a}\right)^3 \right]}_\text{estrutura do variograma}
```

> O efeito pepita $C₀$ não pertence à estrutura do variograma.

O **imbricamento das estruturas** é definido como a soma de $n$ estruturas do variograma:

``` math
γ(h) = C_0 +
C_1 \left[\frac{3h}{2a_1} - \frac{1}{2}\left(\frac{h}{a_1}\right)^3 \right] +
C_2 \left[\frac{3h}{2a_2} - \frac{1}{2}\left(\frac{h}{a_2}\right)^3 \right] + ... + C_n \left[\frac{3h}{2a_n} - \frac{1}{2}\left(\frac{h}{a_n}\right)^3 \right]
```

O imbricamento de estruturas permite uma **maior flexibilidade na modelagem do variograma**. A Figura 11 ilustra, graficamente, um exemplo de modelo de variograma construído a partir do imbricamento de duas estruturas.

"""

# ╔═╡ f95ffa70-f924-404a-8cec-fc281b8588e2
md"""

Efeito pepita: $(@bind nested_cₒ Slider(0.00:0.1:4.0,
										default=3.0, show_value=true))

Contrib. 1ª estrutura: $(@bind nested_c₁ Slider(0.0:0.1:10.0,
													default=2.6, show_value=true))

Contrib. 2ª estrutura: $(@bind nested_c₂ Slider(0.0:0.1:10.0,
													default=2.7, show_value=true))

Alcance 1ª estrutura: $(@bind nested_r₁ Slider(10.0:1.0:156.0, default=83.0, show_value=true)) m

Alcance 2ª estrutura: $(@bind nested_r₂ Slider(10.0:1.0:156.0, default=101.0, show_value=true)) m

"""

# ╔═╡ 3f0465bc-444c-4026-a677-a182366790ae
begin
	
	# Efeito pepita
	γ_nested₀ = NuggetEffect(Float64(nested_cₒ))
	
	# 1ª estrutura
    γ_nested₁ = SphericalVariogram(sill = Float64(nested_c₁),
								   range = Float64(nested_r₁))

	# 2ª estrutura
    γ_nested₂ = SphericalVariogram(sill = Float64(nested_c₂),
								   range = Float64(nested_r₂))

	# Modelo aninhado
    γ_nested  = γ_nested₀ + γ_nested₁ + γ_nested₂
	
	# Plotagem do variograma experimental
	plot(γ₂, color = :orange, legend = false, line = false)
	
	# Plotagem do ajuste teórico
	plot!(γ_nested, color = :black, lw = 2, xlims = (0,220), ylims = (0,15))
	
	# Alcance
	vline!([nested_r₂], color = :gray, ls = :dash)

end

# ╔═╡ 864c9c06-e52b-4de8-bc16-d053fa3c0346
html"""

<p align="center">
    <b>Figura 11</b>: Exemplo de modelo de variograma imbricado.
</p>

"""

# ╔═╡ Cell order:
# ╟─1991a290-cfbf-11eb-07b6-7b3c8543dd28
# ╟─d4775d05-4943-4493-897e-4340f01475be
# ╟─f8909bd5-9167-42ea-a302-a7a50bdc365c
# ╟─029c1951-054b-4f48-bc05-341250ce9f6a
# ╟─51107168-29ca-40b1-a658-9361199be3b1
# ╟─0c00aee8-9db5-4fca-b92d-e19aa4fe5c1b
# ╟─4b12eecc-0645-4f46-b3be-8b8a095af599
# ╟─b23b047e-1c02-40c5-ba88-825da85ba75c
# ╟─8cfef844-5e4d-44c8-817c-0021eecbcaa2
# ╟─528f0bb5-4030-4006-a323-29f9cbc1efc0
# ╟─5e623ea7-03f9-46a9-ba54-6d48d1a64057
# ╟─4b136ca1-f46f-43dc-9a1d-0659f1ef5e61
# ╟─c782a92c-cc4d-44bc-8521-2f70ad222bd5
# ╟─43bc79ba-bb97-48bd-a8e4-c478bdc3a60b
# ╟─3f39dcf2-055e-4aa8-8caa-09223175a5fa
# ╟─facdf937-4056-4699-b505-d9cada0c8ce3
# ╟─5cb62b37-fe28-4816-b7ed-f5f40df255dc
# ╟─f4e189ac-5d12-4de5-80e1-516103e5950f
# ╟─c3135efd-e69c-4270-8b45-b3f9f2dd586c
# ╟─8923e1c1-914d-47b5-a4b4-5f0c53c4e810
# ╟─3d25e2bc-9a6d-4712-92ed-de31dbdea3f2
# ╟─874544a1-07af-4509-a34d-68d77558aaae
# ╟─2965ea0b-9b5e-4460-a735-e555733b2d83
# ╟─9d60b923-72e7-42c8-8ef3-d4a590e3f600
# ╟─7c00b7a2-5c09-46f5-ba8d-03786fd606b8
# ╟─f1f163f7-eabc-4808-82e4-98ecfeddc730
# ╟─8bd1b52b-b6a8-489e-ae74-be2931eef4ee
# ╟─841ffdd2-16b4-4c19-8f03-70942a4ebb2e
# ╟─f738e56e-7725-4d52-a700-960ce372f025
# ╟─650fc66a-3f8e-45d5-a898-5c783a8d12a1
# ╟─ace40206-7ce6-4a64-b1ae-bd19d295158e
# ╟─4c5b95f2-5ad6-4f18-9cc0-9bd96eb3bf29
# ╟─7fa3052f-52c8-48b5-ab3a-8401a6d8f93a
# ╟─728f75bd-0fc5-43c6-9551-4304925aa97b
# ╟─9709372d-3d5f-4bff-8ca1-adbb4dbeda23
# ╟─b5bb3401-48d5-4d22-bbed-06427862062e
# ╟─5e555810-f34d-402c-ac0a-17a423f420bc
# ╟─6433f0dc-04f8-450e-9a94-f8cfa8cda552
# ╟─6f59e6db-0632-483a-89be-6c82dd188d60
# ╟─53979b80-406d-4eb9-aef8-5d3b626a5555
# ╟─e80e5e16-59fb-4ec0-a9f0-6b8b97bc8d36
# ╟─9891913d-e735-4ec8-b09c-49b51417f18d
# ╟─52c19e93-8534-4b59-a164-3f12d23d5440
# ╟─7a8a49a2-007e-409a-9a45-c76f717f75f8
# ╟─f92eedbf-097d-45f6-a550-ccc8c2f9841b
# ╟─ea67e941-496f-45a3-8b0f-058d573291d8
# ╟─03a46b21-c522-4036-a710-bd6ce0a26a1b
# ╟─510759f5-3838-4db7-b683-e39677a2551b
# ╟─c75b37b1-b13a-4d3b-949c-639e4e5dc01b
# ╟─83593f8e-8dd2-40b1-903b-8712bb9eb048
# ╟─a6802bda-7b7a-4d98-bb08-bcbe8a990e01
# ╟─bd8328a0-9cd3-4eec-8c68-f6f15f296f4b
# ╟─ed5011d1-8d0a-4ace-8baa-0af219b3e3e2
# ╟─d1e9297b-04c0-47ad-91e8-77648b205d10
# ╟─e7157b79-368d-4ce1-97d3-22110e3359da
# ╟─341ec3f6-c871-431f-8ffa-85f4c43ae138
# ╟─6d0f5d99-f7e2-4f53-b835-c3b345613e4a
# ╟─61b8631b-8295-4dea-a5dd-189bf578bc8c
# ╟─ab5e6c19-789b-4944-ba8e-f983a9a2652c
# ╟─8b4ee7b2-2a01-44ae-8539-27f1815fe634
# ╟─187c01ca-053e-4994-a748-cf9b16683a50
# ╟─83d6c4fe-bcd6-4d7f-93ef-2e093b1284fa
# ╟─b9634a1e-f225-4986-867f-fd36f56882df
# ╟─ee09dcab-2298-444c-ad9f-f79268c9056c
# ╟─0f28a997-4945-47fe-83b9-058726bc8041
# ╟─f95ffa70-f924-404a-8cec-fc281b8588e2
# ╟─3f0465bc-444c-4026-a677-a182366790ae
# ╟─864c9c06-e52b-4de8-bc16-d053fa3c0346
