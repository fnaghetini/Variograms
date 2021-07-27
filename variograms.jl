### A Pluto.jl notebook ###
# v0.15.1

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
	# load packages used in this notebook
	using GeoStats, GeoStatsImages
	using CSV, DataFrames, Query
    using Statistics, Random
	using PlutoUI
    using Plots
	
	# default plot settings
	gr(format=:png)
end;

# ╔═╡ 3bd915e1-2f58-451c-a0fb-8aec6d6f75d9
PlutoUI.TableOfContents(aside=true, title="Sumário",
						indent=true, depth=3)

# ╔═╡ f8909bd5-9167-42ea-a302-a7a50bdc365c
html"""
<p style="background-color:lightgrey" xmlns:cc="http://creativecommons.org/ns#" xmlns:dct="http://purl.org/dc/terms/"><span property="dct:title">&nbsp&nbsp<b>Variogramas</b></span> por <span property="cc:attributionName">Franco Naghetini</span> é licenciado sob <a href="http://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1"><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1"></a></p>
"""

# ╔═╡ 029c1951-054b-4f48-bc05-341250ce9f6a
md" # Variogramas"

# ╔═╡ 51107168-29ca-40b1-a658-9361199be3b1
md"""

## Função variograma

O **variograma** é uma função matemática que mapeia/descreve a continuidade espacial de uma variável regionalizada (VR).
- *Exemplos de VR:* Au(g/t), Cu(%), Pb(%), Ag(ppb).

O variograma, quando existe, é único e válido para todo o domínio de estimativa.

A **função variograma** pode ser anisotrópica, sendo sensível à direção, mas não ao sentido.
- *Exemplo:* $γ(000°)$ ≠ $γ(045°)$.
- *Exemplo:* $γ(000°)$ = $γ(180°)$.

Função variograma:

```math
2\gamma(h) = \frac{1}{n} \sum_{i=1}^{n} [Z(x_i) - Z(x_i + h)]^2

```

Função semivariograma:

```math
\gamma(h) = \frac{1}{2n} \sum_{i=1}^{n} [Z(x_i) - Z(x_i + h)]^2

```

-  $γ(h)$: valor do variograma $γ(h)$ para uma distância $h$ entre dois pontos.

-  $n$: número de pares de pontos.

-  $Z(x_i)$: valor da variável $Z$ na posição $(x_i)$.

-  $Z(x_i + h)$: valor da variável $Z$ na posição $(x_i+h)$.

> **Nota:** o termo semivariograma foi cunhado para enfatizar o termo $\frac{1}{2n}$ da função. Entretanto, atualmente, ele é considerado obsoleto e, por isso, o termo variograma tende a ser mais utilizado.

Ao final da variografia, teremos em mãos um **modelo de variograma** representativo da continuidade espacial de uma variável e que será utilizado como entrada no sistema linear de krigagem.

"""

# ╔═╡ d4775d05-4943-4493-897e-4340f01475be
# Convert geologic orientation to cartesian orientation
function sph2cart(azi)
	θ = deg2rad(azi)
	sin(θ), cos(θ)
end;

# ╔═╡ 0c00aee8-9db5-4fca-b92d-e19aa4fe5c1b
md"""
## Variogramas experimentais

A resolução da **função variograma** apresentada acima é representada graficamente e recebe o nome de **variograma experimental** (*Figura 1*).

Pode-se unir os pontos do variograma experimental para simplesmente facilitar a sua interpretação.

"""

# ╔═╡ 4b12eecc-0645-4f46-b3be-8b8a095af599
begin
	# Sample image
	image = geostatsimage("Gaussian30x10")

	# Calculating experimental variogram
	γ₁ = EmpiricalVariogram(image, :Z, maxlag=60., nlags = 6)
end;

# ╔═╡ b23b047e-1c02-40c5-ba88-825da85ba75c
md"""

Unir pontos do variograma: $(@bind join_points CheckBox())

"""

# ╔═╡ 8cfef844-5e4d-44c8-817c-0021eecbcaa2
# Ploting experimental variogram
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

Perceba que $h$ é a **distância cartesiana** entre duas amostras, enquanto $γ(h)$ é a **distância estatística** entre essas amostras. Portanto, o variograma experimental é a ferramenta que converte a distância geográfica $h$ em distância estatística $γ(h)$.

Note que a função variograma é uma **função discreta**, ou seja, os valores de $γ(h)$ são calculados apenas para $h$ específicos.

Cada ponto no variograma experimental representa o valor $γ(h)$ médio de um conjunto de pares de amostras separadas por uma distância $h$.

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

- No contexto 2D, informamos apenas o azimute.

- No contexto 3D, informamos o azimute e o mergulho.

No exemplo 2D abaixo, percebemos que quando variamos o azimute, o variograma experimental também sofre uma variação.

Variogramas experimentais que assumem anisotropia são denominados **variogramas direcionais** (*Figura 2*).

"""

# ╔═╡ bdf7046f-f955-446c-8437-f889be9e22c5
begin
	# Defining data URL
	url = "https://raw.githubusercontent.com/fnaghetini/Variograms/main/data/walker_lake_proj.csv"
	
	# Downloading data
	data = download(url)
	
	# Importing Walker Lake dataset
	walker_lake = CSV.File(data, type = Float64) |> DataFrame
	
	# Select just necessary columns
	wl = walker_lake[:,[:X,:Y,:PB]]
	
	# Droping missing values
	dropmissing!(wl)
	
	# Georeferencing data
	wl_georef = georef(wl, (:X,:Y))
end;

# ╔═╡ 43bc79ba-bb97-48bd-a8e4-c478bdc3a60b
md"""

Azimute: $(@bind azm Slider(0:45:180, default=0, show_value = true))°

"""

# ╔═╡ 3f39dcf2-055e-4aa8-8caa-09223175a5fa
begin
	# Directional variogram
    γ₂ = DirectionalVariogram(sph2cart(azm), wl_georef, :PB,
                              maxlag = 200, nlags = 7)
	
	# Ploting directional variogram
	plot(γ₂, legend = false, xlims = (0,200), ylims = (0,15),
		 title = "$(azm)°", color = :orange)
end

# ╔═╡ facdf937-4056-4699-b505-d9cada0c8ce3
html"""

<p align="center">
    <b>Figura 2</b>: Exemplo de variogramas direcionais.
</p>

"""

# ╔═╡ 1704bcbf-004c-4e30-a09a-6325f116b53c
md"""

> **Nota:** repare que $γ(000°) = γ(180°)$.

"""

# ╔═╡ 5cb62b37-fe28-4816-b7ed-f5f40df255dc
md"""

#### Tamanho do passo

O **tamanho do passo (lag)** é a distância média entre as amostras vizinhas na direção em que o variograma experimental está sendo calculado.

Abaixo, tem-se um exemplo de cálculo de variograma experimental na direção E-W para uma **malha amostral regular** (*Figura 3*). Adotou-se um tamanho de passo igual a 1 m.


"""

# ╔═╡ f4e189ac-5d12-4de5-80e1-516103e5950f
md"""

W-E: $(@bind ix₁ Slider(1:1:4, default=1))
N-S: $(@bind iy₁ Slider(1:1:5, default=1))

Passo: $(@bind h Slider(1:1:3, default=1, show_value = true)) m

"""

# ╔═╡ c3135efd-e69c-4270-8b45-b3f9f2dd586c
begin
	# Random seed
	Random.seed!(42)
	
	# Generating samples
	values₁ = DataFrame(Au = rand(25))
	coords₁ = PointSet([(i,j) for i in 1:5 for j in 1:5])
	
	# Georeferencing samples
	samples₁ = georef(values₁, coords₁)

	# Ploting samples
	plot(samples₁, xlims = (0,6), ylims = (0,6), title = "Au (g/t)",
		 xlabel = "X(m)", ylabel = "Y(m)", clims = (0,1))
	
	# Ploting h vector
	plot!([(ix₁,iy₁),(ix₁+h,iy₁)], arrow = true, color = :red, lw = 2)
end

# ╔═╡ 8923e1c1-914d-47b5-a4b4-5f0c53c4e810
html"""

<p align="center">
    <b>Figura 3</b>: Exemplo de busca de pares de amostras em uma malha regular.
</p>

"""

# ╔═╡ 4d6946ec-3d4b-40bd-9ba6-40fb68c80142
md"""

Note que, quanto maior é o número de passos e menor é o tamanho do passo definido, mais o variograma experimental apresenta aspecto de "eletrocardiograma" (*Figura 4*).

"""

# ╔═╡ 472e42d3-00ba-49f9-93ae-438f6cccae08
md"""

N° de passos: $(@bind num_lags Slider(2:1:35, default=6, show_value=true))

"""

# ╔═╡ d0455ea3-9ab0-4545-bb06-ae8809d99290
begin
	# Directional variogram
    γ_lag = DirectionalVariogram(sph2cart(0), wl_georef, :PB,
                              maxlag = 200, nlags = num_lags)
	
	# Ploting directional variogram
	plot(γ_lag, legend = false, xlims = (0,200), ylims = (0,15),
		 title = "Variograma Experimental", color = :orange)
end

# ╔═╡ 1bdc6f70-b3b3-4f9d-aacd-3038fbe3dfa8
html"""

<p align="center">
    <b>Figura 4</b>: Influência do número de passos no variograma experimental.
</p>

"""

# ╔═╡ 3d25e2bc-9a6d-4712-92ed-de31dbdea3f2
md"""

#### Tolerância linear

Sabe-se que, na maioria dos casos, as malhas de sondagem são irregulares. Nesse caso, poucos pares de pontos serão buscados, já que as amostras não se encontram equidistantes entre si (*Figura 5*).

"""

# ╔═╡ 874544a1-07af-4509-a34d-68d77558aaae
md"""

W-E: $(@bind ix₂ Slider(0.0:0.001:1.0, default=0.015))
N-S: $(@bind iy₂ Slider(0.0:0.001:1.0, default=0.172))

Tamanho do passo: $(@bind lag_size Slider(0.05:0.05:1.0, default=0.3,
										  show_value = true)) m

"""

# ╔═╡ 2965ea0b-9b5e-4460-a735-e555733b2d83
begin
	# Random seed
	Random.seed!(42)
	
	# Generating georeferenced samples
	table₂ = georef(values₁, PointSet(rand(2,25)))
	
	# Ploting samples
	plot(table₂, xlims = (-0.2,1.2), ylims = (-0.2,1.2), title = "Au (g/t)",
		 xlabel = "X(m)", ylabel = "Y(m)", clims = (0,1))
	
	# Ploting h vector
	plot!([(ix₂,iy₂),(ix₂+lag_size,iy₂)], arrow = true, color = :red, lw = 2)
end

# ╔═╡ 9d60b923-72e7-42c8-8ef3-d4a590e3f600
html"""

<p align="center">
    <b>Figura 5</b>: Exemplo de busca de pares de amostras em uma malha irregular.
</p>

"""

# ╔═╡ 7c00b7a2-5c09-46f5-ba8d-03786fd606b8
md"""

Uma alternativa comumente utilizada é a definição de uma **tolerância de passo (lagtol)** para que mais amostras sejam buscadas.

A tolerância do passo é definida como metade do tamanho do passo:

```math
lagtol = \frac{lag}{2} 

```

Essa abordagem permite que:

- Não haja overlap de amostras. As amostras são utilizadas apenas em um passo.

- Não haja gap de amostras. As amostras serão sempre utilizadas em algum dos passos.

"""

# ╔═╡ f1f163f7-eabc-4808-82e4-98ecfeddc730
md"""

O exemplo abaixo ilustra amostras colineares de uma **malha amostral irregular** (*Figura 6*). Note que, caso a tolerância de passo não fosse adotada, apenas um par de pontos seria buscado para o cálculo do variograma.

"""

# ╔═╡ 8bd1b52b-b6a8-489e-ae74-be2931eef4ee
begin
	# Random seed
	Random.seed!(42)
	
	# Generating samples
	coords₃ = [(1.,1.),(1.6,1.),(1.9,1.),(2.2,1.),(2.8,1.),(3.6,1.),(3.8,1.)]
	values₃ = DataFrame(Au = rand(7))
	
	# Georeferencing samples
	samples₃ = georef(values₃, coords₃)
end;

# ╔═╡ 841ffdd2-16b4-4c19-8f03-70942a4ebb2e
md"""

W-E: $(@bind ix₃ Slider(1.:0.1:2.8, default=1))

Tolerância de passo: $(@bind lag_tol CheckBox())

"""

# ╔═╡ f738e56e-7725-4d52-a700-960ce372f025
begin
	
	if lag_tol
		# Ploting samples
		plot(samples₃, xlims = (0.,4.5), ylims = (0,2), title = "Au (g/t)",
			 xlabel = "X(m)", ylabel = "Y(m)", clims = (0,1))
		
		# Ploting h vector
		plot!([(ix₃,1.),(ix₃+1,1.)], arrow = true, color = :red, lw = 2)

		# Ploting lag - ½ lag dashed line
		vline!([ix₃ + 0.5], color = :gray, ls = :dash)
		annotate!(ix₃ + 0.5, 2.1, text("lag - ½ lag", 7, :gray))

		# Ploting lag solid line
		vline!([ix₃+1], color = :red)
		annotate!(ix₃+1, 2.1, text("lag", 7, :red))

		# Ploting lag + ½ lag dashed line
		vline!([ix₃ + 1.5], color = :gray, ls = :dash)
		annotate!(ix₃ + 1.5, 2.1, text("lag + ½ lag", 7, :gray))
		
	else
		# Ploting samples
		plot(samples₃, xlims = (0.,4.5), ylims = (0,2), title = "Au (g/t)",
			 xlabel = "X(m)", ylabel = "Y(m)", clims = (0,1))

		# Ploting h vector
		plot!([(ix₃,1.),(ix₃+1,1.)], arrow = true, color = :red, lw = 2)

		# Ploting lag solid line
		vline!([ix₃+1], color = :red)
		annotate!(ix₃+1, 2.1, text("lag", 7, :red))
	end
		
end

# ╔═╡ 650fc66a-3f8e-45d5-a898-5c783a8d12a1
html"""

<p align="center">
    <b>Figura 6</b>: Exemplo de amostras colineares não equidistantes.
</p>

"""

# ╔═╡ ace40206-7ce6-4a64-b1ae-bd19d295158e
md"""

#### Tolerância angular

Além da tolerância linear (i.e. tolerância de passo), deve-se também definir uma **tolerância angular** que, por sua vez, é dividida em dois parâmetros:

- Tolerância de azimute

- Tolerância de mergulho

A **tolerância angular (angtol)** é definida como metade do **incremento angular (anginc)**:

```math
angtol = \frac{anginc}{2} 

```

A convenção acima permite que:

- Não haja overlap de amostras.

- Não haja gap de amostras.


> **Nota:** a fórmula acima é válida tanto para a tolerância de azimute quanto para a tolerância de mergulho.

> **Nota:** para malhas irregulares, é necessária a definição das tolerâncias de azimute e de mergulho em um contexto 3D, mas apenas a definição da tolerância de azimute em um contexto 2D.

"""

# ╔═╡ 049568a8-0d02-403d-9d87-ce9a5bf5e242
md"""

A *Figura 7* ilustra um exemplo de tolerâncias angulares (i.e. azimute e mergulho) para um incremento angular de 45°.

"""

# ╔═╡ 728f75bd-0fc5-43c6-9551-4304925aa97b
md"""

Azimute: $(@bind dip_dir Slider(0:45:180, default=0, show_value = true))°

Mergulho: $(@bind dip Slider(0:45:180, default=0, show_value = true))°

"""

# ╔═╡ 4c5b95f2-5ad6-4f18-9cc0-9bd96eb3bf29
begin
	# Ploting dip direction solid line
	dip_dir_tol = plot([sph2cart(dip_dir+180),sph2cart(dip_dir)],
					   color = :red, ticks = false, xlims = (-1,1),
					   ylims = (-1,1), grid = false, lw = 3,
					   arrow = true, axis = false, legend = false,
					   title = "Tolerância de Azimute", size = (300,300))
	
	# Ploting dashed dip direction tolerance lines
	plot!([sph2cart(dip_dir+180-22.5),sph2cart(dip_dir-22.5)], color = :gray,
		  ls = :dash)
	
	plot!([sph2cart(dip_dir+180+22.5),sph2cart(dip_dir + 22.5)], color = :gray,
		  ls = :dash)
	
	# Ploting N-S line
	vline!([0], color = :black)
	
	# Ploting W-E line
	hline!([0], color = :black)
end;

# ╔═╡ 7fa3052f-52c8-48b5-ab3a-8401a6d8f93a
begin
	# Ploting dip solid line
	dip_tol = plot([sph2cart(dip+270),sph2cart(dip+90)], color = :red,
				   ticks = false, xlims = (-1,1), ylims = (-1,1),
				   grid = false, lw = 3, arrow = true, axis = false,
				   legend = false, title = "Tolerância de Mergulho",
				   size = (300,300))
	
	# Ploting dashed dip tolerance lines
	plot!([sph2cart(dip+270-22.5),sph2cart(dip+90-22.5)], color = :gray,
		  ls = :dash)
	
	plot!([sph2cart(dip+270+22.5),sph2cart(dip+90+22.5)], color = :gray,
		  ls = :dash)
	
	# Ploting W-E line
	hline!([0], color = :black)
end;

# ╔═╡ 9709372d-3d5f-4bff-8ca1-adbb4dbeda23
# Ploting dip direction and dip
plot(dip_dir_tol, dip_tol, layout = (1,2), size = (600,300))

# ╔═╡ b5bb3401-48d5-4d22-bbed-06427862062e
html"""

<p align="center">
    <b>Figura 7</b>: Tolerâncias angulares.
</p>

"""

# ╔═╡ 5e555810-f34d-402c-ac0a-17a423f420bc

md"""
#### Largura da banda

A **largura da banda** é um parâmetro de restrição facultativo que pode ser utilizado em conjunto com a tolerância angular.

Esse parâmetro é definido pela distância entre a reta da direção de cálculo do variograma experimental e a linha de tolerância angular. A partir de uma distância igual à largura da banda a busca por amostras se paraleliza (*Figura 8*).

"""

# ╔═╡ 6433f0dc-04f8-450e-9a94-f8cfa8cda552
html"""

<p align="center">
    <img src="https://i.postimg.cc/7Zb0qPyc/Figure-01.png" style="height : 300px">
</p>

<p align="center">
    <b>Figura 8</b>: Largura da banda.
</p>

"""

# ╔═╡ 6f59e6db-0632-483a-89be-6c82dd188d60
md"""

### Malhas regulares x malhas irregulares

Para **malhas amostrais regulares**, os parâmetros de variograma experimental a serem definidos são:

- Número de direções de cálculo

    - *Exemplo:* 4 direções

- Direções de cálculo

    - *Exemplo:* 000°, 045°, 090°, 135°

- Tamanho do passo

    - *Exemplo:* 25 m

- Número de passos

    - *Exemplo:* 5 passos (25 m, 50 m, 75 m, 100 m, 125 m)

No caso de **malhas amostrais irregulares**, os parâmetros de variograma experimental a serem definidos são:

- Número de direções de cálculo

    - *Exemplo:* 4 direções

- Direções de cálculo e tolerância angular

    - *Exemplo:* 000° ± 22.5°, 045° ± 22.5°, 090° ± 22.5°, 135° ± 22.5°

- Tamanho do passo e tolerância do passo

    - *Exemplo:* 25 ± 12.5 m

- Número de passos

    - *Exemplo:* 5 passos (25 m, 50 m, 75 m, 100 m, 125 m)

- Largura da banda

    - *Exemplo:* 30 m

"""

# ╔═╡ e80e5e16-59fb-4ec0-a9f0-6b8b97bc8d36
md"""

## Modelos teóricos

A partir dos variogramas experimentais só é possível obter valores médios de variograma $γ(h)$ para distâncias iguais a múltiplos do tamanho de passo $h$ escolhido.

Portanto, é necessário o ajuste de um **modelo matemático contínuo**, de modo que saberemos o valor do variograma $γ(h)$ para qualquer distância entre pares de amostras $h$.

O procedimento de se ajustar um modelo teórico contínuo ao variograma experimental é denominado **modelagem do variograma** (*Figura 9*).

"""

# ╔═╡ 9891913d-e735-4ec8-b09c-49b51417f18d
# Fitting experimental variogram 
varmod = fit(GaussianVariogram, γ₁);

# ╔═╡ 52c19e93-8534-4b59-a164-3f12d23d5440
md"""
Ajuste do modelo: $(@bind fit_model CheckBox())

"""

# ╔═╡ 7a8a49a2-007e-409a-9a45-c76f717f75f8
begin
	if fit_model
		# Ploting experimental variogram
		plot(γ₁, color = :orange, line = false)
	
		# Ploting variogram model
		plot!(varmod, ylims = (0,1.0), xlims = (0,50), label = "Ajuste Teórico",
		  	  legend = false)

	else
		# Ploting experimental variogram
		plot(γ₁, color = :orange, line = false, ylims = (0,1.0), xlims = (0,50),
			 legend = false)
	end
end	

# ╔═╡ f92eedbf-097d-45f6-a550-ccc8c2f9841b
html"""

<p align="center">
    <b>Figura 9</b>: Exemplo de ajuste de um modelo teórico a um variograma experimental.
</p>

"""

# ╔═╡ ea67e941-496f-45a3-8b0f-058d573291d8
md"""

### Elementos do variograma

Para se modelar o variograma, deve-se configurar quatro elementos:

- Alcance

- Efeito Pepita

- Patamar

- Tipo de Modelo

"""

# ╔═╡ 03a46b21-c522-4036-a710-bd6ce0a26a1b
md"""

##### Alcance (a)

O **alcance (_range_)** consiste na distância máxima até onde se consegue estabelecer alguma interdependência espacial entre pares de amostras.

Em outras palavras, o alcance define até qual distância $h$ existe correlação espacial entre pares de amostras. Portanto:

> **Nota:** para distâncias $h$ superiores ao alcance, não há interdependência espacial entre amostras.

"""

# ╔═╡ 510759f5-3838-4db7-b683-e39677a2551b
md"""

##### Efeito Pepita (C₀)

O **efeito pepita (_nugget effect_)** é definido como a descontinuidade próxima a origem do variograma. É o valor de $γ(h)$ quando $h$ tende a zero. Perceba que:

> **Nota:** $γ(h) = C_0, \forall h → 0$

> **Nota:** $γ(h) = 0, \forall h = 0$

"""

# ╔═╡ c75b37b1-b13a-4d3b-949c-639e4e5dc01b
md"""

##### Patamar (C + C₀)

O **patamar (_sill_)** é definido como o máximo valor de $γ(h)$ que as amostras podem apresentar.

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

O **tipo de modelo (_model_)** controla a variabilidade a pequenas distâncias, ou seja, o comportamento próximo à origem do variograma.


"""

# ╔═╡ a6802bda-7b7a-4d98-bb08-bcbe8a990e01
md"""

###### Modelo Gaussiano

- Apresenta comportamento próximo à origem parabólico.

- Esse tipo de modelo teórico é normalmente utilizado para ajustar variogramas experimentais de fenômenos de baixa heterogeneidade.

- Sua equação é descrita como:

``` math

γ(h) = C_0 + C \left[ 1 - exp \left[- \left(\frac{h}{a} \right)^2 \right]  \right] 

```

"""


# ╔═╡ bd8328a0-9cd3-4eec-8c68-f6f15f296f4b
md"""

###### Modelo Esférico

- Apresenta comportamento próximo à origem linear.

- Esse tipo de modelo teórico é normalmente utilizado para ajustar variogramas experimentais de fenômenos de intermediária heterogeneidade.

- É o modelo mais utilizado na indústria da mineração.

- Sua equação é descrita como:

``` math
γ(h) = C_0 + C \left[\frac{3h}{2a} - \frac{1}{2}\left(\frac{h}{a}\right)^3 \right], ∀ h < a)
```

``` math
γ(h) = C_0 + C, ∀ h ≥ a
```

"""

# ╔═╡ ed5011d1-8d0a-4ace-8baa-0af219b3e3e2
md"""

###### Modelo Exponencial

- Apresenta comportamento próximo à origem linear. Entretanto, a inclinação desse modelo nas proximidades da origem é maior do que a inclinação do Modelo Esférico.

- Esse tipo de modelo teórico é normalmente utilizado para ajustar variogramas experimentais de fenômenos de elevada heterogeneidade.

- Sua equação é descrita como:

``` math
γ(h) = C_0 + C \left[1 - exp \left[-\left(\frac{h}{a} \right) \right] \right]
```

"""

# ╔═╡ d1e9297b-04c0-47ad-91e8-77648b205d10
md"""

Termos das equações:

-  $γ(h)$ = valor do variograma para a distância $h$
-  $C_0$ = efeito pepita
-  $h$ = passo
-  $C$ = contribuição ao patamar
-  $a$ = alcance

"""

# ╔═╡ e7157b79-368d-4ce1-97d3-22110e3359da
md"""

A *Figura 10* ilustra graficamente os quatro elementos que constituem o variograma. 

"""

# ╔═╡ 6d0f5d99-f7e2-4f53-b835-c3b345613e4a
md"""

Modelo Teórico: $(@bind model Select(["Gaussiano","Esférico","Pentaesférico","Exponencial"],
		default = "Gaussiano"))

Ef. Pepita: $(@bind c₀ Slider(0.0:0.1:0.5, default=0.0, show_value=true))

Patamar: $(@bind cₜ Slider(0.5:0.1:1.0, default=1.0, show_value=true))

Alcance: $(@bind a Slider(5.0:10:45.0, default=25.0, show_value=true)) m

"""

# ╔═╡ 341ec3f6-c871-431f-8ffa-85f4c43ae138
# Gaussian variogram model
if model == "Gaussiano"
	γ₃ = GaussianVariogram(nugget = Float64(c₀),
						   sill = Float64(cₜ),
		 				   range = Float64(a))

# Spherical variogram model
elseif model == "Esférico"
	γ₃ = SphericalVariogram(nugget = Float64(c₀),
						    sill = Float64(cₜ),
						    range = Float64(a))

# Pentaspherical variogram model
elseif model == "Pentaesférico"
	γ₃ = PentasphericalVariogram(nugget = Float64(c₀),
						   		 sill = Float64(cₜ),
						   		 range = Float64(a))

# Exponential variogram model
else
	γ₃ = ExponentialVariogram(nugget = Float64(c₀),
						      sill = Float64(cₜ),
						      range = Float64(a))

end;

# ╔═╡ 61b8631b-8295-4dea-a5dd-189bf578bc8c
begin
	# Ploting variogram model
	plot(γ₃, color = :black, lw = 2, label = model,
		 legend = :topleft, ylims = (0.,1.5), xlims = (0.,65.))
	
	# Ploting nugget effect dashed line
	hline!([c₀], color = :red, ls = :dash, label = "Ef. Pepita")
	annotate!(55,c₀+0.05,text("Ef. Pepita",10,:red))
	
	# Ploting sill dashed line
	hline!([cₜ], color = :green, ls = :dash, label = "Patamar")
	annotate!(55,cₜ+0.05,text("Patamar",10,:green))
	
	# Ploting range dashed line
	vline!([a], color = :blue, ls = :dash, label = "Alcance")
	annotate!(a,-0.05,text("Alcance",10,:blue))
end

# ╔═╡ ab5e6c19-789b-4944-ba8e-f983a9a2652c
html"""

<p align="center">
    <b>Figura 10</b>: Elementos do variograma.
</p>

"""

# ╔═╡ 8b4ee7b2-2a01-44ae-8539-27f1815fe634
md"""

## Tipos de anisotropia

Na geoestatística, a **anisotropia** existe quando um ou mais elementos do variograma variam com a mudança da direção. Existem três tipos (*Figura 11*):

- *Anisotropia Zonal*: patamar varia de acordo com a mudança de direção.

- *Anisotropia Geométrica*: alcance varia de acordo com a mudança de direção.

- *Anisotropia Mista*: patamar e alcance variam de acordo com a mudança de direção.

> **Nota:** embora existam três tipos de anisotropia, é comum considerar apenas a anisotropia geométrica para a modelagem do variograma.

> **Nota:** não existe anisotropia de efeito pepita, uma vez que esse elemento é, por definição, isotópico.

"""

# ╔═╡ 83d6c4fe-bcd6-4d7f-93ef-2e093b1284fa
md"""

Tipo de anisotropia: $(@bind aniso Select(["Zonal","Geométrica","Mista"],
										  default = "Geométrica"))

"""

# ╔═╡ 187c01ca-053e-4994-a748-cf9b16683a50
# Zonal anisotropy
if aniso == "Zonal"
	γ_aniso₁ = SphericalVariogram(nugget = 0.1, range = 50.0, sill = 1.0)
	γ_aniso₂ = SphericalVariogram(nugget = 0.1, range = 50.0, sill = 1.5)

# Geometric anisotropy
elseif aniso == "Geométrica"
	γ_aniso₁ = SphericalVariogram(nugget = 0.1, range = 50.0, sill = 1.0)
	γ_aniso₂ = SphericalVariogram(nugget = 0.1, range = 30.0, sill = 1.0)

# Mixed anisotropy
else
	γ_aniso₁ = SphericalVariogram(nugget = 0.1, range = 50.0, sill = 1.5)
	γ_aniso₂ = SphericalVariogram(nugget = 0.1, range = 30.0, sill = 1.0)

end;

# ╔═╡ b9634a1e-f225-4986-867f-fd36f56882df
begin
	# Ploting red variogram model
	plot(γ_aniso₁, color = :red, lw = 2, legend = false,
		 title = "Anisotropia $aniso")
	
	# Ploting blue variogram model
	plot!(γ_aniso₂, color = :blue, lw = 2, xlims = (0,80), ylims = (0,2))
end

# ╔═╡ ee09dcab-2298-444c-ad9f-f79268c9056c
html"""

<p align="center">
    <b>Figura 11</b>: Tipos de anisotropia.
</p>

"""

# ╔═╡ 0f28a997-4945-47fe-83b9-058726bc8041
md"""

## Estruturas imbricadas

A **estrutura do variograma** é a porção da equação do ajuste teórico em que o valor de $C$ cresce com o aumento da distância $h$.

``` math
γ(h) = C_0 +
\underbrace{C \left[\frac{3h}{2a} - \frac{1}{2}\left(\frac{h}{a}\right)^3 \right]}_\text{estrutura do variograma}
```

> **Nota:** o efeito pepita $C_0$ não pertence à estrutura do variograma.

O **imbricamento/aninhamento das estruturas** é definido como a soma de $n$ estruturas do variograma. A equação abaixo ilustra um imbricamento de $n$ estruturas para um modelo esférico:

``` math
γ(h) = C_0 +
\underbrace{C_1 \left[\frac{3h}{2a_1} - \frac{1}{2}\left(\frac{h}{a_1}\right)^3 \right]}_\text{1ª estrutura} +
\underbrace{C_2 \left[\frac{3h}{2a_2} - \frac{1}{2}\left(\frac{h}{a_2}\right)^3 \right]}_\text{2ª estrutura} + ... +
\underbrace{C_n \left[\frac{3h}{2a_n} - \frac{1}{2}\left(\frac{h}{a_n}\right)^3 \right]}_\text{n-ésima estrutura}
```

O patamar ($C$) consiste na soma entre todas as contribuições ao patamar e o efeito pepita:

```math
C = C_0 + C_1 + C_2 + ... + C_n
```

> **Nota:** normalmente, utiliza-se, no máximo, três estruturas imbricadas em um modelo de variograma.

O imbricamento de estruturas permite uma maior flexibilidade na modelagem do variograma. A *Figura 12* ilustra, graficamente, um exemplo de modelo de variograma construído a partir do imbricamento de duas estruturas.

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
	# Nugget structure
	γ_nested₀ = NuggetEffect(Float64(nested_cₒ))
	
	# 1ª structure
    γ_nested₁ = SphericalVariogram(sill = Float64(nested_c₁),
								   range = Float64(nested_r₁))

	# 2ª structure
    γ_nested₂ = SphericalVariogram(sill = Float64(nested_c₂),
								   range = Float64(nested_r₂))

	# Nested model
    γ_nested  = γ_nested₀ + γ_nested₁ + γ_nested₂
	
	# Ploting experimental variogram
	plot(γ₂, color = :orange, legend = false, line = false)
	
	# Ploting variogram model
	plot!(γ_nested, color = :black, lw = 2, xlims = (0,220), ylims = (0,15),
		  title = "Variograma Imbricado")
	
	# Ploting range dashed line
	vline!([nested_r₂], color = :gray, ls = :dash)
end

# ╔═╡ 864c9c06-e52b-4de8-bc16-d053fa3c0346
html"""

<p align="center">
    <b>Figura 12</b>: Exemplo de modelo de variograma imbricado.
</p>

"""

# ╔═╡ 538bf67b-33c6-45c3-b5bf-328922debb26
md"""

## Variograma anisotrópico

Como a continuidade espacial de fenômenos naturais tende a ser anisotrópica e o objetivo da variografia é justamente descrever a continuidade espacial desses fenômenos, é plausível que o variograma seja anisotrópico.

A forma mais simples e coerente para se representar anisotropia é por meios de elipses (2D) ou elipsoides (3D).

Em um contexto 3D, assumindo condições de **anisotropia geométrica**, para representar a continuidade espacial de um fenômeno, basta encontrarmos os eixos principais do elipsoide, de modo que:

- O efeito pepita será isotrópico.

- O patamar será assumido como isotrópico.

- O alcance será anisotrópico.

Portanto, os eixos do elipsoide representam justamente a variação do alcance para diferentes direções.

A equação de um **modelo esférico anisotrópico** é descrita como:

``` math
γ(h) = C_0 + C \left[\frac{3h}{2(a_x,a_y,a_z)} - \frac{1}{2}\left(\frac{h}{(a_x,a_y,a_z)}\right)^3 \right]
```

A *Figura 13* ilustra graficamente um exemplo de modelo de variograma anisotrópico.

"""

# ╔═╡ 18282939-e7ef-4da4-aade-72e7b01886de
md"""

Alcance em Y: $(@bind range_y Slider(10.0:2.0:120.0, default=100.0, show_value=true)) m

Alcance em X: $(@bind range_x Slider(10.0:2.0:120.0, default=66.0, show_value=true)) m

Alcance em Z: $(@bind range_z Slider(10.0:2.0:120.0, default=26.0, show_value=true)) m

"""

# ╔═╡ dc47965d-e732-44e4-875c-b4922ff4bd1f
begin
	# Primary model variogram
	γ_1st = SphericalVariogram(nugget = 0.1, range = Float64(range_y), sill = 5.0)
	
	# Secondary model variogram
	γ_2nd = SphericalVariogram(nugget = 0.1, range = Float64(range_x), sill = 5.0)
	
	# Tertiary model variogram
	γ_3rd = SphericalVariogram(nugget = 0.1, range = Float64(range_z), sill = 5.0)
end;

# ╔═╡ d6a4e6dd-7ace-4406-be57-804b4c2537e5
md"""

| Estrutura |Ef. Pepita | Alcance em X | Alcance em Y | Alcance em Z | Variância |
|:---------:|:---------:|:------------:|:------------:|:------------:| -------|
| 0 | 0.01 | - | - | - | - |
| 1 | - | $(range_x) m | $(range_y) m | $(range_z) m | 5.0 |

"""

# ╔═╡ b2ea2e47-4fa5-4d17-8341-889069a717c7
begin
	# Ploting primary model variogram
	plot(γ_1st, color = :red, lw = 2, label = "Primário",
		 title = "Variograma Anisotrópico")
	
	# Ploting secondary model variogram
	plot!(γ_2nd, color = :green, lw = 2, label = "Secundário")
	
	# Ploting tertiary model variogram
	plot!(γ_3rd, color = :blue, lw = 2, label = "Terciário",
		  xlims = (0,120), ylims = (0,8))
end

# ╔═╡ 7e05a32f-44ba-45ec-8db2-6d23a966a298
html"""

<p align="center">
    <b>Figura 13</b>: Exemplo de modelo de variograma anisotrópico.
</p>

"""

# ╔═╡ 6feb0cb4-7bff-4635-ae38-4400affe89f3
md"""

## Modelo de variograma x estimativas

Sabe-se que o modelo de variograma é utilizado como entrada na estimativa por krigagem. Nesse sentido, cada um de seus parâmetros e elementos exerce uma influência no modelo de teores estimados:

- A **direção** indica a orientação da continuidade espacial de teores.

- O **alcance** controla o comprimento de continuidade espacial médio ("elipses" na imagem).

- O **patamar** define a "altura" das "elipses".

- O **efeito pepita** define uma variabilidade adicional para escalas menores.

- O **tipo de modelo** define o comportamento próximo a origem.

O exemplo abaixo auxilia na compreensão da influência de cada um desses parâmetros nas estimativas resultantes. A *Figura 14* mostra o modelo de variograma anisotrópico utilizado na estimativa por krigagem. A *Figura 15* representa o mapa da localização das amostras.

"""

# ╔═╡ 8079a74c-005d-4654-8e44-d763a12aefd8
md"""

Direção de maior continuidade: $(@bind azi₂ Slider(0.0:45.0:90.0, default=0.0, show_value=true))°

Modelo Teórico: $(@bind m Select(["Gaussiano","Esférico","Exponencial"],
							 default = "Esférico"))

Efeito pepita: $(@bind Cₒ Slider(0.00:0.1:5.0, default=3.0, show_value=true))

Alcance primário (Y): $(@bind ry Slider(10.0:1.0:156.0, default=101.0, show_value=true)) m

Alcance secundário (X): $(@bind rx Slider(10.0:1.0:156.0, default=32.0, show_value=true)) m

"""

# ╔═╡ 39e7cb17-7813-4103-880d-64803c636039
begin
	# Theoretical variogram model
	model_type = Dict("Gaussiano" => GaussianVariogram,
					  "Esférico" => SphericalVariogram,
					  "Exponencial" => ExponentialVariogram)
	
	# Sample variance (defining variogram sill)
	σ² = var(wl_georef[:PB])
	
	# Calculating primary experimental variogram
	γexp_pri = DirectionalVariogram(sph2cart(azi₂), wl_georef, :PB,
                                    maxlag = 200, nlags = 5)
	
	# Calculating secondary experimental variogram
	γexp_sec = DirectionalVariogram(sph2cart(azi₂+90), wl_georef, :PB,
                                    maxlag = 200, nlags = 5)
	
	# Fitting primary experimental variogram 
	γm_pri = model_type[m](nugget = Float64(Cₒ),
						   sill = Float64(σ²),
						   range = Float64(ry))
	
	# Fitting secondary experimental variogram 
	γm_sec = model_type[m](nugget = Float64(Cₒ),
						   sill = Float64(σ²),
						   range = Float64(rx))
end;

# ╔═╡ 308abd53-d536-4ff0-8e1d-9ac118742d93
begin
	# Graphic parameters
	col_pri = :red
	col_sec = :blue
	xlim    = (0,200)
	ylim    = (0,15)
	
	# Ploting primary experimental variogram
	plot(γexp_pri, color = col_pri, label = false)
	
	# Ploting secondary experimental variogram
	plot!(γexp_sec, color = col_sec, label = false)
	
	# Ploting primary variogram model
	plot!(γm_pri, color = col_pri, lw = 2, legend = false)
	
	# Ploting secondary variogram model
	plot!(γm_sec, color = col_sec, lw = 2, xlims = xlim, ylims = ylim)
	
	# Ploting sill dashed line
	hline!([σ²], color = :gray, ls = :dash)
	
	# Ploting primary range dashed line
	vline!([ry], color = col_pri, ls = :dash)
	
	# Ploting secondary range dashed line
	vline!([rx], color = col_sec, ls = :dash)
end

# ╔═╡ a0b3b930-5f2a-47a1-bc81-c70c2ff595e6
html"""

<p align="center">
    <b>Figura 14</b>: Modelo de variograma anisotrópico utilizado na estimativa.
</p>

"""

# ╔═╡ ee0d2529-a6bc-4ee3-bd74-27a38a585c14
md"""

Quando a checkbox `Filtrar apenas high grades` é marcada, apenas os high grades de Pb (%) são apresentados. Esses valores correspondem a todos os teores de Pb maiores que o P90. Isso facilita a visualização espacial das amostras.

Além disso, ao ativar a checkbox `Visualizar estimativas`, é possível visualizar as estimativas realizadas a partir do modelo de variograma configurado (*Figura 13*). As estimativas serão mostradas na *Figura 15*.

"""

# ╔═╡ fb99bba7-e81b-4653-a7dc-3558f6fc7e2c
md"""

Filtrar apenas high grades: $(@bind filter_hg CheckBox())

Visualizar estimativas: $(@bind show_model CheckBox())

"""

# ╔═╡ cd5c821d-247e-4d18-93cf-065197b79f1b
begin
	if show_model
		# Defining anisotropy ellipsoid
		ellip = Ellipsoid([ry,rx],[azi₂], convention = GSLIB)

		# Anisotropic variogram model
		γ = model_type[m](nugget = Float64(Cₒ),
						  sill = Float64(σ²),
						  distance = metric(ellip))

		# Estimation domain
		dom = CartesianGrid((243,283),(8.,8.),(1.,1.))

		# Defining estimation problem
		problem = EstimationProblem(wl_georef, dom, :PB)

		# Defining solver
		OK = Kriging(:PB => (variogram = γ,
						     neighborhood = ellip,
						     minneighbors = 8,
							 maxneighbors = 16)
				    )

		# Solving estimation problem
		sol = solve(problem, OK)
		
		# Manipulating estimates
		estimates = sol |> @map({PB = _.PB, geometry = _.geometry}) |> GeoData
	end
end;

# ╔═╡ c90bdb75-1918-4d57-93b1-6cda3e8fbb0c
begin
	if show_model
		# Ploting estimates		
		plot(estimates, color=:coolwarm, xlabel="X", ylabel="Y",
			 xlims=(8,251), ylims=(8,291),clims = (0,12),
			 marker=(:square,1.2), markerstrokewidth=0,
			 size=(500,500))
		
		# Ploting samples
		plot!(wl_georef, color=:coolwarm, marker=(:square,2),
			  markerstrokecolor=:black, markerstrokewidth=0.3,
		      title="Pb (%)")
		
	else
		
		if filter_hg
			# Filtering high grades
			wl_filt = wl_georef |> @filter(_.PB > quantile(wl.PB, 0.9)) |> GeoData
			
			# Ploting high grades
			plot(wl_filt, color=:coolwarm, marker=(:square,2),
				 markerstrokecolor=:black, markerstrokewidth=0.3,
				 xlims=(8,251), ylims=(8,291),clims = (0,12),
				 size=(500,500),title="Pb (%)", xlabel="X", ylabel="Y")

		else
			# Ploting samples
			plot(wl_georef, color=:coolwarm, marker=(:square,2),
				 markerstrokecolor=:black, markerstrokewidth=0.3,
				 xlims=(8,251), ylims=(8,291),clims = (0,12),
				 size=(500,500),title="Pb (%)", xlabel="X", ylabel="Y")
		end
		
	end
end

# ╔═╡ 2f1d77a0-e5cd-4d77-8031-cff161f67a45
html"""

<p align="center">
    <b>Figura 15</b>: Mapa de localização das amostras de Pb (%). Ative a primeira checkbox para visualizar apenas os high grades ou a segunda para visualizar as estimativas.
</p>

"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
GeoStats = "dcc97b0b-8ce5-5539-9008-bb190f959ef6"
GeoStatsImages = "7cd16168-b42c-5e7d-a585-4f59d326662d"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Query = "1a8c2f83-1ff3-5112-b086-8aa67b057ba1"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
CSV = "~0.8.5"
DataFrames = "~1.2.1"
GeoStats = "~0.25.2"
GeoStatsImages = "~0.6.0"
Plots = "~1.19.3"
PlutoUI = "~0.7.9"
Query = "~1.0.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "485ee0867925449198280d4af84bdb46a2a404d0"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.0.1"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgCheck]]
git-tree-sha1 = "dedbbb2ddb876f899585c4ec4433265e3017215a"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.1.0"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArrayInterface]]
deps = ["IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "a71d224f61475b93c9e196e83c17c6ac4dedacfa"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.1.18"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[AverageShiftedHistograms]]
deps = ["LinearAlgebra", "RecipesBase", "Statistics", "StatsBase", "UnicodePlots"]
git-tree-sha1 = "2e16019c442fd57011a440f6c7693aef81fadc8d"
uuid = "77b51b56-6f8f-5c3a-9cb4-d71f9594ea6e"
version = "0.8.6"

[[BangBang]]
deps = ["Compat", "ConstructionBase", "Future", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables", "ZygoteRules"]
git-tree-sha1 = "e239020994123f08905052b9603b4ca14f8c5807"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.31"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c3598e525718abcc440f69cc6d5f60dda0a1b61e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.6+5"

[[CSV]]
deps = ["Dates", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode"]
git-tree-sha1 = "b83aa3f513be680454437a0eee21001607e5d983"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.8.5"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "e2f47f6d8337369411569fd45ae5753ca10394c6"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.0+6"

[[CategoricalArrays]]
deps = ["DataAPI", "Future", "JSON", "Missings", "Printf", "RecipesBase", "Statistics", "StructTypes", "Unicode"]
git-tree-sha1 = "1562002780515d2573a4fb0c3715e4e57481075e"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.10.0"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "ad613c934ec3a3aa0ff19b91f15a16d56ed404b5"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.0.2"

[[CircularArrays]]
deps = ["OffsetArrays"]
git-tree-sha1 = "786e067a47d43952f483f4eb107c02f59ea6f93e"
uuid = "7a955b69-7140-5f4e-a0ed-f168c5e2e749"
version = "1.2.0"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random", "StaticArrays"]
git-tree-sha1 = "ed268efe58512df8c7e224d2e170afd76dd6a417"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.13.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "dc7dedc2c2aa9faf59a55c622760a25cbefbe941"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.31.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[CompositionsBase]]
git-tree-sha1 = "455419f7e328a1a2493cabc6428d79e951349769"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.1"

[[ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "32d125af0fb8ec3f8935896122c5e345709909e5"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.0"

[[Crayons]]
git-tree-sha1 = "3f71217b538d7aaee0b69ab47d9b7724ca8afa0d"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.0.4"

[[DataAPI]]
git-tree-sha1 = "ee400abb2298bd13bfc3df1c412ed228061a2385"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.7.0"

[[DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "a19645616f37a2c2c3077a44bc0d3e73e13441d7"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.2.1"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4437b64df1e0adccc3e5d1adbc3ac741095e4677"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.9"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[DataValues]]
deps = ["DataValueInterfaces", "Dates"]
git-tree-sha1 = "d88a19299eba280a6d062e135a43f00323ae70bf"
uuid = "e7dc6d0d-1eca-5fa6-8ad6-5aecde8b7ea5"
version = "0.4.13"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DefineSingletons]]
git-tree-sha1 = "77b4ca280084423b728662fe040e5ff8819347c5"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.1"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DensityRatioEstimation]]
deps = ["LinearAlgebra", "Parameters", "Requires", "Statistics", "StatsBase"]
git-tree-sha1 = "9fd052a2fab80e2e033c9c28dc014a35aaf7da4b"
uuid = "ab46fb84-d57c-11e9-2f65-6f72e4a7229f"
version = "0.4.3"

[[Dictionaries]]
deps = ["Indexing", "Random"]
git-tree-sha1 = "011d10589449681ad631c0f866a2ce72771648c6"
uuid = "85a47980-9c8c-11e8-2b9f-f7ca1fa99fb4"
version = "0.3.10"

[[DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[DiffRules]]
deps = ["NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "214c3fcac57755cfda163d91c58893a8723f93e9"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.0.2"

[[Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "abe4ad222b26af3337262b8afb28fab8d215e9f8"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.3"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns"]
git-tree-sha1 = "3889f646423ce91dd1055a76317e9a1d3a23fff1"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.11"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "a32185f5428d3986f47c2ab78b1f216d5e6cc96f"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.5"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "92d8f9f208637e8d2d28c664051a00569c01493d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.1.5+1"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "LibVPX_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "3cc57ad0a213808473eafef4845a74766242e05f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.3.1+4"

[[FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "f985af3b9f4e278b1d24434cbb546d6092fca661"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.3"

[[FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3676abafff7e4ff07bbd2c42b3d8201f31653dcc"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.9+8"

[[FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "256d8e6188f3f1ebfa1a5d17e072a0efafa8c5bf"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.10.1"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "693210145367e7685d8604aee33d9bfb85db8b31"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.11.9"

[[FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "8b3c09b56acaf3c0e581c66638b85c8650ee9dca"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.8.1"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "35895cf184ceaab11fd778b4590144034a167a2f"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.1+14"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "NaNMath", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "e2af66012e08966366a43251e1fd421522908be6"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.18"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "cbd58c9deb1d304f5a245a0b7eb841a2560cfec6"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.1+5"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "dba1e8614e98949abfa60480b13653813d8f0157"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+0"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "9f473cdf6e2eb360c576f9822e7c765dd9d26dbc"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.58.0"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "eaf96e05a880f3db5ded5a5a8a7817ecba3c7392"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.58.0+0"

[[GeoClustering]]
deps = ["Distances", "LinearAlgebra", "Meshes", "Statistics", "TableOperations", "Tables"]
git-tree-sha1 = "30199d291102f123faf8086df6d9ab08abec5a99"
uuid = "7472b188-6dde-460e-bd07-96c4bc049f7e"
version = "0.1.5"

[[GeoEstimation]]
deps = ["Distances", "GeoStatsBase", "KrigingEstimators", "LinearAlgebra", "Meshes", "NearestNeighbors", "Variography"]
git-tree-sha1 = "38b9d31eaa49265a9deb75d068c2a1dbd0d3dd4c"
uuid = "a4aa24f8-9f24-4d1a-b848-66d123bfa54d"
version = "0.7.2"

[[GeoLearning]]
deps = ["Distributions", "GeoStatsBase", "MLJModelInterface", "Meshes", "TableOperations", "Tables"]
git-tree-sha1 = "af31350de7c005e643cc13e3f253e6619c43db58"
uuid = "90c4468e-a93e-43b4-8fb5-87d804bc629f"
version = "0.1.4"

[[GeoSimulation]]
deps = ["CpuId", "Distributions", "FFTW", "GeoStatsBase", "KrigingEstimators", "LinearAlgebra", "Meshes", "Statistics", "Variography"]
git-tree-sha1 = "d2ce0390e030299af5d7c18c39eff8d43de39efa"
uuid = "220efe8a-9139-4e14-a4fa-f683d572f4c5"
version = "0.3.6"

[[GeoStats]]
deps = ["DensityRatioEstimation", "Distances", "GeoClustering", "GeoEstimation", "GeoLearning", "GeoSimulation", "GeoStatsBase", "KrigingEstimators", "LossFunctions", "Meshes", "PointPatterns", "Reexport", "Variography"]
git-tree-sha1 = "3a31406b580b2b435940d894edf72d48380586df"
uuid = "dcc97b0b-8ce5-5539-9008-bb190f959ef6"
version = "0.25.2"

[[GeoStatsBase]]
deps = ["AverageShiftedHistograms", "CSV", "CategoricalArrays", "Combinatorics", "DensityRatioEstimation", "Distances", "Distributed", "Distributions", "LinearAlgebra", "LossFunctions", "MLJModelInterface", "Meshes", "Optim", "Parameters", "RecipesBase", "ScientificTypes", "StaticArrays", "Statistics", "StatsBase", "TableOperations", "Tables", "Transducers", "TypedTables"]
git-tree-sha1 = "929d3c7fe92ac9f61383d097d311188bfc0b30b9"
uuid = "323cb8eb-fbf6-51c0-afd0-f8fba70507b2"
version = "0.21.9"

[[GeoStatsImages]]
deps = ["FileIO", "GslibIO"]
git-tree-sha1 = "1e1c08019257940fd78e6aed4644299c143991f9"
uuid = "7cd16168-b42c-5e7d-a585-4f59d326662d"
version = "0.6.0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "15ff9a14b9e1218958d3530cc288cf31465d9ae2"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.3.13"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "47ce50b742921377301e15005c96e979574e130b"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.1+0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[GslibIO]]
deps = ["DelimitedFiles", "FileIO", "GeoStatsBase", "Meshes", "Printf", "Tables"]
git-tree-sha1 = "e5193f014eead4e5334e34212831e8523c36e7d0"
uuid = "4610876b-9b01-57c8-9ad9-06315f1a66a5"
version = "0.7.7"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "c6a1fff2fd4b1da29d3dccaffb1e1001244d844e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.12"

[[IfElse]]
git-tree-sha1 = "28e837ff3e7a6c3cdb252ce49fb412c8eb3caeef"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.0"

[[Indexing]]
git-tree-sha1 = "ce1566720fd6b19ff3411404d4b977acd4814f9f"
uuid = "313cdc1a-70c2-5d6a-ae34-0150d3930a38"
version = "1.1.1"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[InitialValues]]
git-tree-sha1 = "26c8832afd63ac558b98a823265856670d898b6c"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.2.10"

[[IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[InvertedIndices]]
deps = ["Test"]
git-tree-sha1 = "15732c475062348b0165684ffe28e85ea8396afc"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.0.0"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IterableTables]]
deps = ["DataValues", "IteratorInterfaceExtensions", "Requires", "TableTraits", "TableTraitsUtils"]
git-tree-sha1 = "70300b876b2cebde43ebc0df42bc8c94a144e1b4"
uuid = "1c8ee90f-4401-5389-894e-7a04a3dc0f4d"
version = "1.0.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "81690084b6198a2e1da36fcfda16eeca9f9f24e4"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.1"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[KrigingEstimators]]
deps = ["Combinatorics", "GeoStatsBase", "LinearAlgebra", "Meshes", "Statistics", "Variography"]
git-tree-sha1 = "979f8ce302f0c2fa06403efbcdd4b785082363ff"
uuid = "d293930c-a38c-56c5-8ebb-12008647b47a"
version = "0.7.5"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "c7f1c695e06c01b95a67f0cd1d34994f3e7db104"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.2.1"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a4b12a1bd2ebade87891ab7e36fdbce582301a92"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.6"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LearnBase]]
git-tree-sha1 = "a0d90569edd490b82fdc4dc078ea54a5a800d30a"
uuid = "7f8f8fb0-2700-5f03-b4bd-41f8cfc144b6"
version = "0.4.1"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[LibVPX_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "12ee7e23fa4d18361e7c2cde8f8337d4c3101bc7"
uuid = "dd192d2f-8180-539f-9fb4-cc70b1dcf69a"
version = "1.10.0+0"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "761a393aeccd6aa92ec3515e428c26bf99575b3b"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+0"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "f27132e551e959b3667d8c93eae90973225032dd"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.1.1"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["DocStringExtensions", "LinearAlgebra"]
git-tree-sha1 = "7bd5f6565d80b6bf753738d2bc40a5dfea072070"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.2.5"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[LossFunctions]]
deps = ["InteractiveUtils", "LearnBase", "Markdown", "RecipesBase", "StatsBase"]
git-tree-sha1 = "0f057f6ea90a84e73a8ef6eebb4dc7b5c330020f"
uuid = "30fc2ffe-d236-52d8-8643-a9d8f7c094a7"
version = "0.7.2"

[[MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "c253236b0ed414624b083e6b72bfe891fbd2c7af"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2021.1.1+1"

[[MLJModelInterface]]
deps = ["Random", "ScientificTypesBase", "StatisticalTraits"]
git-tree-sha1 = "55c785a68d71c5fd7b64b490e0d9ab18cf13a04c"
uuid = "e80e1ace-859a-464e-9ed9-23947d8ae3ea"
version = "1.1.1"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "6a8a2a625ab0dea913aba95c11370589e0239ff0"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.6"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Meshes]]
deps = ["CategoricalArrays", "CircularArrays", "Distances", "IterTools", "IteratorInterfaceExtensions", "LinearAlgebra", "NearestNeighbors", "Random", "RecipesBase", "ReferenceFrameRotations", "SimpleTraits", "SparseArrays", "SpecialFunctions", "StaticArrays", "StatsBase", "TableTraits", "Tables"]
git-tree-sha1 = "8ca246dd3d13cc5187ba75ae8f0253c3740b3038"
uuid = "eacbb407-ea5a-433e-ab97-5258b1ca43fa"
version = "0.16.11"

[[MicroCollections]]
deps = ["BangBang", "Setfield"]
git-tree-sha1 = "e991b6a9d38091c4a0d7cd051fcb57c05f98ac03"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.0"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f8c673ccc215eb50fcadb285f522420e29e69e1c"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "0.4.5"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50608f411a1e178e0129eab4110bd56efd08816f"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.0"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "16baacfdc8758bc374882566c9187e785e85c2f0"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.9"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "4f825c6da64aebaa22cc058ecfceed1ab9af1c7e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.3"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7937eda4681660b4d6aeeecc2f7e1c81c8ee4e2f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+0"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Optim]]
deps = ["Compat", "FillArrays", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "d34366a3abc25c41f88820762ef7dfdfe9306711"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.3.0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "4dd403333bcf0909341cfe57ec115152f937d7d8"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.1"

[[Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "2276ac65f1e236e0a6ea70baff3f62ad4c625345"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.2"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "c8abc88faa3f7a3950832ac5d6e690881590d6dc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "1.1.0"

[[PersistenceDiagramsBase]]
deps = ["Compat", "Tables"]
git-tree-sha1 = "ec6eecbfae1c740621b5d903a69ec10e30f3f4bc"
uuid = "b1ad91c1-539c-4ace-90bd-ea06abc420fa"
version = "0.1.1"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "501c20a63a34ac1d015d5304da0e645f42d91c9f"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.0.11"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs"]
git-tree-sha1 = "1bbbb5670223d48e124b388dee62477480e23234"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.19.3"

[[PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "Suppressor"]
git-tree-sha1 = "44e225d5837e2a2345e69a1d1e01ac2443ff9fcb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.9"

[[PointPatterns]]
deps = ["Distributions", "GeoStatsBase", "Meshes"]
git-tree-sha1 = "aa6c09ea2efe8d782b3e6b64b68013103e52bbb2"
uuid = "e61b41b6-3414-4803-863f-2b69057479eb"
version = "0.3.8"

[[PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "cde4ce9d6f33219465b55162811d8de8139c0414"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.2.1"

[[PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "0d1245a357cc61c8cd61934c07447aa569ff22e6"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.1.0"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "12fbe86da16df6679be7521dfb39fbc861e1dc7b"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.1"

[[Query]]
deps = ["DataValues", "IterableTables", "MacroTools", "QueryOperators", "Statistics"]
git-tree-sha1 = "a66aa7ca6f5c29f0e303ccef5c8bd55067df9bbe"
uuid = "1a8c2f83-1ff3-5112-b086-8aa67b057ba1"
version = "1.0.0"

[[QueryOperators]]
deps = ["DataStructures", "DataValues", "IteratorInterfaceExtensions", "TableShowUtils"]
git-tree-sha1 = "911c64c204e7ecabfd1872eb93c49b4e7c701f02"
uuid = "2aef5ad7-51ca-5a8f-8e88-e75cf067b44b"
version = "0.9.3"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "b3fb709f3c97bfc6e948be68beeecb55a0b340ae"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.1"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "2a7a2469ed5d94a98dea0e85c46fa653d76be0cd"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.3.4"

[[Reexport]]
git-tree-sha1 = "5f6c21241f0f655da3952fd60aa18477cf96c220"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.1.0"

[[ReferenceFrameRotations]]
deps = ["Crayons", "LinearAlgebra", "Printf", "StaticArrays"]
git-tree-sha1 = "fecac02781f5c475c957d8088c4b43a0a44316b5"
uuid = "74f56ac7-18b3-5285-802d-d4bd4f104033"
version = "1.0.0"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[ScientificTypes]]
deps = ["CategoricalArrays", "ColorTypes", "Dates", "PersistenceDiagramsBase", "PrettyTables", "ScientificTypesBase", "StatisticalTraits", "Tables"]
git-tree-sha1 = "345e33061ad7c49c6e860e42a04c62ecbea3eabf"
uuid = "321657f4-b219-11e9-178b-2701a2544e81"
version = "2.0.0"

[[ScientificTypesBase]]
git-tree-sha1 = "3f7ddb0cf0c3a4cff06d9df6f01135fa5442c99b"
uuid = "30f210dd-8aff-4c5f-94ba-8e64358c1161"
version = "1.0.0"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "35927c2c11da0a86bcd482464b93dadd09ce420f"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.5"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "d5640fc570fb1b6c54512f0bd3853866bd298b3e"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.7.0"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "LogExpFunctions", "OpenSpecFun_jll"]
git-tree-sha1 = "508822dca004bf62e210609148511ad03ce8f1d8"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.6.0"

[[SplitApplyCombine]]
deps = ["Dictionaries", "Indexing"]
git-tree-sha1 = "3cdd86a92737fa0c8af19aecb1141e71057dc2db"
uuid = "03a91e81-4c3e-53e1-a0a4-9c0c8f19dd66"
version = "1.1.4"

[[SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "edef25a158db82f4940720ebada14a60ef6c4232"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.13"

[[Static]]
deps = ["IfElse"]
git-tree-sha1 = "62701892d172a2fa41a1f829f66d2b0db94a9a63"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.3.0"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "5b2f81eeb66bcfe379947c500aae773c85c31033"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.8"

[[StatisticalTraits]]
deps = ["ScientificTypesBase"]
git-tree-sha1 = "5114841829816649ecc957f07f6a621671e4a951"
uuid = "64bff920-2084-43da-a3e6-9bb72801c0c9"
version = "2.0.0"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "2f6792d523d7448bbe2fec99eca9218f06cc746d"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.8"

[[StatsFuns]]
deps = ["LogExpFunctions", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "30cd8c360c54081f806b1ee14d2eecbef3c04c49"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.8"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "000e168f5cc9aded17b6999a560b7c11dda69095"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.0"

[[StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "e36adc471280e8b346ea24c5c87ba0571204be7a"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.7.2"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableOperations]]
deps = ["SentinelArrays", "Tables", "Test"]
git-tree-sha1 = "a7cf690d0ac3f5b53dd09b5d613540b230233647"
uuid = "ab02a1b2-a7df-11e8-156e-fb1833f50b87"
version = "1.0.0"

[[TableShowUtils]]
deps = ["DataValues", "Dates", "JSON", "Markdown", "Test"]
git-tree-sha1 = "14c54e1e96431fb87f0d2f5983f090f1b9d06457"
uuid = "5e66a065-1f0a-5976-b372-e0b8c017ca10"
version = "0.2.5"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[TableTraitsUtils]]
deps = ["DataValues", "IteratorInterfaceExtensions", "Missings", "TableTraits"]
git-tree-sha1 = "8fc12ae66deac83e44454e61b02c37b326493233"
uuid = "382cd787-c1b6-5bf2-a167-d5b971a19bda"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "8ed4a3ea724dac32670b062be3ef1c1de6773ae8"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.4.4"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "34f27ac221cb53317ab6df196f9ed145077231ff"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.65"

[[TypedTables]]
deps = ["Adapt", "SplitApplyCombine", "Tables", "Unicode"]
git-tree-sha1 = "19f9fca769be841fdf42241a4daa2ab52a11f89b"
uuid = "9d95f2ec-7b3d-5a63-8d20-e2491e220bb9"
version = "1.3.1"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[UnicodePlots]]
deps = ["Crayons", "Dates", "SparseArrays", "StatsBase"]
git-tree-sha1 = "1a63e6eea76b291378ff9f95801f8b6d96213208"
uuid = "b8865327-cd53-5732-bb35-84acbb429228"
version = "1.3.0"

[[Variography]]
deps = ["Distances", "GeoStatsBase", "InteractiveUtils", "LinearAlgebra", "Meshes", "NearestNeighbors", "Optim", "Parameters", "Printf", "RecipesBase", "Setfield", "SpecialFunctions", "Transducers"]
git-tree-sha1 = "7353264bababd8616c09ea4ee101f022bb5fb473"
uuid = "04a0146e-e6df-5636-8d7f-62fa9eb0b20c"
version = "0.12.9"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll"]
git-tree-sha1 = "2839f1c1296940218e35df0bbb220f2a79686670"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.18.0+4"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "9e7a1e8ca60b742e508a315c17eef5211e7fbfd7"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.1"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "acc685bcf777b2202a904cdcb49ad34c2fa1880c"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.14.0+4"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7a5780a0d9c6864184b3a2eeeb833a0c871f00ab"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "0.1.6+4"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "c45f4e40e7aafe9d086379e5578947ec8b95a8fb"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d713c1ce4deac133e3334ee12f4adff07f81778f"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2020.7.14+2"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "487da2f8f2f0c8ee0e83f39d13037d6bbf0a45ab"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.0.0+3"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─1991a290-cfbf-11eb-07b6-7b3c8543dd28
# ╟─3bd915e1-2f58-451c-a0fb-8aec6d6f75d9
# ╟─f8909bd5-9167-42ea-a302-a7a50bdc365c
# ╟─029c1951-054b-4f48-bc05-341250ce9f6a
# ╟─51107168-29ca-40b1-a658-9361199be3b1
# ╟─d4775d05-4943-4493-897e-4340f01475be
# ╟─0c00aee8-9db5-4fca-b92d-e19aa4fe5c1b
# ╟─4b12eecc-0645-4f46-b3be-8b8a095af599
# ╟─b23b047e-1c02-40c5-ba88-825da85ba75c
# ╟─8cfef844-5e4d-44c8-817c-0021eecbcaa2
# ╟─528f0bb5-4030-4006-a323-29f9cbc1efc0
# ╟─5e623ea7-03f9-46a9-ba54-6d48d1a64057
# ╟─4b136ca1-f46f-43dc-9a1d-0659f1ef5e61
# ╟─c782a92c-cc4d-44bc-8521-2f70ad222bd5
# ╟─bdf7046f-f955-446c-8437-f889be9e22c5
# ╟─43bc79ba-bb97-48bd-a8e4-c478bdc3a60b
# ╟─3f39dcf2-055e-4aa8-8caa-09223175a5fa
# ╟─facdf937-4056-4699-b505-d9cada0c8ce3
# ╟─1704bcbf-004c-4e30-a09a-6325f116b53c
# ╟─5cb62b37-fe28-4816-b7ed-f5f40df255dc
# ╟─f4e189ac-5d12-4de5-80e1-516103e5950f
# ╟─c3135efd-e69c-4270-8b45-b3f9f2dd586c
# ╟─8923e1c1-914d-47b5-a4b4-5f0c53c4e810
# ╟─4d6946ec-3d4b-40bd-9ba6-40fb68c80142
# ╟─d0455ea3-9ab0-4545-bb06-ae8809d99290
# ╟─472e42d3-00ba-49f9-93ae-438f6cccae08
# ╟─1bdc6f70-b3b3-4f9d-aacd-3038fbe3dfa8
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
# ╟─049568a8-0d02-403d-9d87-ce9a5bf5e242
# ╟─7fa3052f-52c8-48b5-ab3a-8401a6d8f93a
# ╟─728f75bd-0fc5-43c6-9551-4304925aa97b
# ╟─9709372d-3d5f-4bff-8ca1-adbb4dbeda23
# ╟─b5bb3401-48d5-4d22-bbed-06427862062e
# ╟─5e555810-f34d-402c-ac0a-17a423f420bc
# ╟─6433f0dc-04f8-450e-9a94-f8cfa8cda552
# ╟─6f59e6db-0632-483a-89be-6c82dd188d60
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
# ╟─538bf67b-33c6-45c3-b5bf-328922debb26
# ╟─dc47965d-e732-44e4-875c-b4922ff4bd1f
# ╟─18282939-e7ef-4da4-aade-72e7b01886de
# ╟─d6a4e6dd-7ace-4406-be57-804b4c2537e5
# ╟─b2ea2e47-4fa5-4d17-8341-889069a717c7
# ╟─7e05a32f-44ba-45ec-8db2-6d23a966a298
# ╟─6feb0cb4-7bff-4635-ae38-4400affe89f3
# ╟─39e7cb17-7813-4103-880d-64803c636039
# ╟─8079a74c-005d-4654-8e44-d763a12aefd8
# ╟─308abd53-d536-4ff0-8e1d-9ac118742d93
# ╟─a0b3b930-5f2a-47a1-bc81-c70c2ff595e6
# ╟─ee0d2529-a6bc-4ee3-bd74-27a38a585c14
# ╟─cd5c821d-247e-4d18-93cf-065197b79f1b
# ╟─fb99bba7-e81b-4653-a7dc-3558f6fc7e2c
# ╟─c90bdb75-1918-4d57-93b1-6cda3e8fbb0c
# ╟─2f1d77a0-e5cd-4d77-8031-cff161f67a45
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
