### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 76823222-2137-466c-9be2-d251a1132de3
begin
	using Combinatorics
	using LinearAlgebra
	using Optim
	using Random
end

# ╔═╡ 6d5700d0-22ec-11f0-3b35-ab9c2a365742
md""" # Final de Optimización
## Distribución de máxima entropía usando dualidad

Para resolver el problema de encontrar la distribución de máxima entropía, donde $\cal{M} = \{ A \in \cal{P}([n]) : \# A=k \}$ 
y $p:\cal{M} \rightarrow \mathbb{R}$. Es decir, nuestro espacio de búsqueda serán las distribuciones de probabilidad sobre los subconjuntos de $[n]$ de cardinal $k$, que a su vez satisfacen las condiciones de distribución marginal: la probabilidad de que un elemento $e \in [n]$ sea seleccionado es $\theta_{e}$. El hecho de que sean distribuciones de probabilidad y que verifiquen las condiciones marginales se visualiza en las restricciones del problema.

#### Problema primal
$$max \sum_{M \in \cal{M}}{p_M \space ln(\frac{1}{p_M})}$$
$$s.a.: \sum_{\substack{M \in \cal{M} \\ e \in M}}{p_M} = \theta_{e} \space \space \forall e \in [n]$$
$$\sum_{M \in \cal{M}}{p_{M}} = 1$$
$$p_{M} \geq 0 \space \space \forall M \in \cal{M}$$


#### Problema dual
$$min \sum_{e \in [n]}{\theta_{e} \lambda_{e}} + ln(\sum_{M \in \cal{M}}{e^{\lambda(M)}})$$
$$s.a.: \lambda_{e} \in \mathbb{R} \space \space \forall e \in [m]$$
"""

# ╔═╡ abbb0b74-6d5d-4e93-8ba4-207d0ecc95ef
md"""
Realizo corridas sobre vectores $\theta$, $\theta'$ $\in \mathbb{R}^{n}$ tales que suman k, $\theta_i \leq \theta'_i \space$ $\forall i \in \{1, \dots, k\} \space$ y $\space \theta_i \geq \theta'_i$ $\forall i \in \{k+1, \dots, n\}$, para encontrar un contraejemplo a $\textit{selection monotonicity}$.
"""

# ╔═╡ 21307643-4656-4e3d-a9bf-2fde973dad66
md"""Defino como p1 el contraejemplo encontrado por Jamie para $n = 6$ y $k = 3$
"""

# ╔═╡ c895cd78-6fd0-472d-a51f-eafd0891ba57
begin
	p1 = [
    Rational{BigInt}(BigInt(50725394825993519278096800656275153334288384969061806535166055966652175553215110626112154355363727),
                     BigInt(820344544708777686628383921246921807685840301005725685964582767932498825215270734395260460182489056)),

    Rational{BigInt}(BigInt(16995635448467104183757178014765702738718813669542372343274191569647329181122211511135709477229515),
                     BigInt(820344544708777686628383921246921807685840301005725685964582767932498825215270734395260460182489056)),

    Rational{BigInt}(BigInt(16995635448467104183757178014765702738718813669542372343274191569647329181122211511135709477229515),
                     BigInt(820344544708777686628383921246921807685840301005725685964582767932498825215270734395260460182489056)),

    Rational{BigInt}(BigInt(203732522694899282198523379819523750032272292460663481324189224126535182285019873726167269517595483),
                     BigInt(205086136177194421657095980311730451921460075251431421491145691983124706303817683598815115045622264)),

    Rational{BigInt}(BigInt(203732522694899282198523379819523750032272292460663481324189224126535182285019873726167269517595483),
                     BigInt(205086136177194421657095980311730451921460075251431421491145691983124706303817683598815115045622264)),

    Rational{BigInt}(BigInt(746456786844211074651353568498768863987616551023722656078520071679268183450193679728059651096880547),
                     BigInt(820344544708777686628383921246921807685840301005725685964582767932498825215270734395260460182489056))
]

end

# ╔═╡ 07aa2842-28b2-4dd7-a4b4-6d2b9f9fd0ad
begin
	Random.seed!(1234)
	n2 = 6
	k2 = 3
	M2 = collect(CoolLexCombinations(n2,k2))
	M_bool2 = []
	
	for A in M2
		a = zeros(n2)
		a[A] .= 1
		push!(M_bool2, a)
	end	
end

# ╔═╡ 089ae76f-a0c7-4213-93ed-76b4a6532e1a
begin
	# Función para perturbar el vector θ de manera aleatoria y generar θ'
	function perturb_vector_rand(p::Vector{Float64})
	    admisible = false
		d = length(p)
		delta = rand()/2
		
		while !admisible
			p_aux = copy(p)
			# Total amount to shift
		    shift_total = delta * k2  # increase k2 elements by delta
		    shift_per_rest = shift_total / (d - k2)
		
		    for i in 1:k2
		        p_aux[i] += delta
		    end
		
		    for i in (k2+1):d
		        p_aux[i] -= shift_per_rest
		    end
		
		    # Check if the result is still valid
		    if all(x -> 0.0 ≤ x < 1.0, p_aux) && abs(sum(p_aux) - sum(p)) < 1e-8
				admisible = true
				return(p_aux)
		    else
		        delta = delta*0.85
		    end
		end
		return p_new
	end
end

# ╔═╡ 0821fac1-0b9f-4311-81e1-03c6531db730
begin
	function f2(λ, p)
	    suma = 0
		for A in M_bool2
		suma += ℯ^(-dot(λ, A))
		end
		logaritmo = log(suma)
		return dot(λ, p) + logaritmo
	end

	function ∇f2(λ, p)
		suma1 = 0
		for A in M_bool2
			suma1 += ℯ^(-dot(λ, A))
		end
		logaritmo = log(suma1)

		g = zeros(n2)
		for i in 1:n2
			suma = 0
			for A in M_bool2
				if A[i] == 1
					suma += ℯ^(-dot(λ, A))
				end
			end
			g[i] = suma
		end
		
		return p - g*suma1
	end
end

# ╔═╡ a2794b94-0d64-4fd6-b528-bf5752e8c588
md"""Calculo la distribución de máxima entropía para $p1$ y la probabilidad de que Conditional Poisson Rounding seleccione B con esa distribución.

"""

# ╔═╡ 82b47801-51bc-4651-bdbc-7ebc5d620735
begin
	p2 = Float64.(p1)
	res_p = optimize(λ -> f2(λ, p2),-50ones(n2),50ones(n2), zeros(n2),SAMIN(),Optim.Options(iterations=1000000))
	λ2_opt = Optim.minimizer(res_p)
	z2 = 0    # Cte de normalización para p
	
	for A in M_bool2
		z2 += ℯ^(-dot(λ2_opt, A))
	end

	B = zeros(n2)
	B[1:k2] .= 1
	p_B = ℯ^(-dot(λ2_opt, B)) * (1/z2)
end

# ╔═╡ c67600a0-5544-48ea-97a7-9db11ef31e79
begin
	p2_prima = zeros(n2)
	encontrado = false
	intentos = 0
	λ2_prima_opt = zeros(n2)

	while !encontrado
		intentos += 1
		println("Intentos: ", intentos)
		
		# Defino p'
		p2_prima = perturb_vector_rand(p2)
	
		# Resuelvo para p'
		res_p_prima = optimize(λ -> f2(λ, p2_prima),-50ones(n2),50ones(n2), zeros(n2),SAMIN(),Optim.Options(iterations=1000000))
		λ2_prima_opt = Optim.minimizer(res_p_prima)
		#println("Distribución óptima λ para p:", λ2_opt)
		#println("Distribución óptima λ' para p':", λ2_prima_opt)
	
		# Calculo P(S={1,...,k}) para p'
		local z2_prima = 0     # Cte de normalización para p'
		for A in M_bool2
			z2_prima += ℯ^(-dot(λ2_prima_opt, A))
		end
	
		p_B_prima = ℯ^(-dot(λ2_prima_opt, B)) * (1/z2_prima)
	
		encontrado = p_B > p_B_prima
	end

	println("Contraejemplo encontrado:")
	println("Vector p:", p2)
	println("Vector p':", p2_prima)	
end

# ╔═╡ a00ed0dc-1fe7-45f4-92cd-650b61604226
begin
	zz2_prima = 0     # Cte de normalización para p'

	for A in M_bool2
		zz2_prima += ℯ^(-dot(λ2_prima_opt, A))
	end

	p_B_prima = ℯ^(-dot(λ2_prima_opt, B)) * (1/zz2_prima)

	m = length(M2)
	probabilidades2 = zeros(m)
	probabilidades2_prima = zeros(m)
	for i in 1:m
		probabilidades2[i] = ℯ^(-dot(λ2_opt, M_bool2[i]))
		probabilidades2_prima[i] = ℯ^(-dot(λ2_prima_opt, M_bool2[i]))
	end

	probabilidades2 = probabilidades2*(1/z2)
	probabilidades2_prima = probabilidades2*(1/zz2_prima)

	println("El vector de probabilidades es ", probabilidades2)
	println("La suma de probabilidades para p es ", sum(probabilidades2))
	println("La suma de probabilidades para p' es ", sum(probabilidades2_prima))
end

# ╔═╡ b2711614-5fa4-4451-89a7-e220fa084d9b
begin
	println("Contraejemplo encontrado:")
	println("N = ", n2, "   k = ", k2)
	println("Vector p:", p2)
	println("Vector p':", p2_prima)	
	println("Distribución óptima λ para p:", λ2_opt)
	println("Distribución óptima λ' para p':", λ2_prima_opt)
	println("Probabilidades para B = {1,...,k} según λ y λ':", p_B, ";  ", p_B_prima)
end

# ╔═╡ 701ed612-cce5-49b9-850b-f124d2de3457
begin
	marginales2 = zeros(n2)
	marginales2_prima = zeros(n2)
	for i in 1:n2
		for j in 1:m
			if (M_bool2[j])[i] == 1
				marginales2[i] += probabilidades2[j]
				marginales2_prima[i] += probabilidades2[j]
			end
		end
	end

	println("Las marginales dadas por los residuos son ", p2)
	println(p2_prima)
	println("\nLas marginales calculadas con la distribución óptima son ", marginales2)
	println(marginales2_prima)
	println("\nLa norma 2 de la diferencia entre las marginales es ", norm(p2-marginales2))
	println(norm(p2_prima-marginales2_prima))
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Combinatorics = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Optim = "429524aa-4258-5aef-a3af-852621145aeb"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[compat]
Combinatorics = "~1.0.2"
Optim = "~1.9.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.2"
manifest_format = "2.0"
project_hash = "0c56d13732ac6a6035bf41c1be1731cc677bae87"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "50c3c56a52972d78e8be9fd135bfb91c9574c140"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.1.1"

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

    [deps.Adapt.weakdeps]
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "ed2ec3c9b483842ae59cd273834e5b46206d6dda"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.11.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.0+0"

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

    [deps.FillArrays.weakdeps]
    PDMats = "90014a1f-27ba-587c-ab20-58faa44d9150"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "2de436b72c3422940cbe1367611d137008af7ec3"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.23.1"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

    [deps.ForwardDiff.weakdeps]
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "be3dc50a92e5a386872a493a10050136d4703f9b"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.6.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "d9b79c4eed437421ac4285148fcadf42e0700e89"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.9.4"

    [deps.Optim.extensions]
    OptimMOIExt = "MathOptInterface"

    [deps.Optim.weakdeps]
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.OrderedCollections]]
git-tree-sha1 = "12f1439c4f986bb868acda6ea33ebc78e19b95ad"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.7.0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─6d5700d0-22ec-11f0-3b35-ab9c2a365742
# ╠═76823222-2137-466c-9be2-d251a1132de3
# ╟─abbb0b74-6d5d-4e93-8ba4-207d0ecc95ef
# ╠═089ae76f-a0c7-4213-93ed-76b4a6532e1a
# ╟─21307643-4656-4e3d-a9bf-2fde973dad66
# ╠═c895cd78-6fd0-472d-a51f-eafd0891ba57
# ╠═07aa2842-28b2-4dd7-a4b4-6d2b9f9fd0ad
# ╠═0821fac1-0b9f-4311-81e1-03c6531db730
# ╟─a2794b94-0d64-4fd6-b528-bf5752e8c588
# ╠═82b47801-51bc-4651-bdbc-7ebc5d620735
# ╠═c67600a0-5544-48ea-97a7-9db11ef31e79
# ╠═a00ed0dc-1fe7-45f4-92cd-650b61604226
# ╠═b2711614-5fa4-4451-89a7-e220fa084d9b
# ╠═701ed612-cce5-49b9-850b-f124d2de3457
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
