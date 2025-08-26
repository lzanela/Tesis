### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ eb417b8f-8824-4263-bcf3-7e722a2291f6
begin
	using Random
	using StatsBase: Weights
	using Combinatorics
	using LinearAlgebra
end

# ╔═╡ a29751ad-5eab-4988-ae35-0c6b098014ba
begin
	using StatsBase
	# Función para samplear de N con probabilidades dadas por q	
	function sample_from_subset(N::Set{Int}, q::Vector{Float64})
	    # Aseguramos que los índices estén en el rango válido
	    elements = collect(N)
	    weights = q[elements]  # Extraemos las probabilidades correspondientes solo a N
	    return StatsBase.sample(elements, Weights(weights))
	end
end

# ╔═╡ 040604b0-3b6e-11f0-35ff-c10fb2433828
md""" # Algoritmo de Brewer
## Búsqueda de contraejemplo para selection monotonicity
"""

# ╔═╡ 57af78d6-c62c-41a4-90dd-6f66cd91122a
begin
	n = 7
	k = 3
end

# ╔═╡ 1d7b2b4a-2a8b-40bd-9048-5d7f3224f63b
begin
	# Función para generar un vector θ que cumpla con los requisitos enunciados
	function rand_vector_sum(dim::Int=8, total::Float64=3.0)
	    while true
	        # Generate dim-1 sorted uniform random numbers in (0, total)
	        cuts = sort(rand(dim - 1) .* total)
	        # Add 0 and total to make full partition
	        parts = diff([0.0; cuts; total])
	        # Retry if any component ≥ 1
	        if all(x -> x < 1.0, parts)
	            return parts
	        end
	    end
	end

	# Función para perturbar el vector θ y generar θ'
	function perturb_vector(p::Vector{Float64}, delta::Float64=0.01)
	  admisible = false
		d = length(p)

		while !admisible
			p_aux = copy(p)
			# Total amount to shift
		  shift_total = delta * k  # increase 3 elements by delta
		  shift_per_rest = shift_total / (d - k)
		  for i in 1:k
		      p_aux[i] += delta
		  end
		  for i in (k+1):d
		      p_aux[i] -= shift_per_rest
		  end
		  # Check if the result is still valid
		  if all(x -> 0.0 ≤ x < 1.0, p_aux) && abs(sum(p_aux) - sum(p)) < 1e-10
			admisible = true
		    return(p_aux)
		  else
		    delta = delta*0.95
		  end
		end
		return p_new
	end
end

# ╔═╡ 1c84a304-7e2a-4038-abbb-f06c0585754d
md"""
### Algoritmo de Brewer
"""

# ╔═╡ bbeb909b-c55e-4232-ba91-d5e8f3a06266
begin
	function brewer(k, p)
		N = Set(1:length(p))
		elegidos = Set{Int}()
		for t in k:-1:1
			suma = sum(p[i] for i in N)
			q = p * (1/suma)
			probabilidades = q .* (1 .- q) ./ (1 .- t .* q)
			for i in 1:n
			    if i ∉ N
			        probabilidades[i] = 0.0
			    end
			end
			probabilidades = probabilidades * (1/sum(probabilidades))
			elemento = sample_from_subset(N, probabilidades)
			delete!(N, elemento)
			push!(elegidos, elemento)
		end
		return elegidos
	end
end

# ╔═╡ 1c8f6a7d-3308-473f-86ec-ab3189d7b77a
begin
	p1 = rand_vector_sum(n,float(k))
	elegidos = brewer(k, p1)
	println(collect(elegidos))
end

# ╔═╡ 28b2efaa-bd96-429b-8a38-042b71f7e9a8
md"""
### Probabilidad de que Brewer seleccione a un conjunto A particular de cardinal k:
$P_B(S=A)$
"""

# ╔═╡ 140ec3dc-8169-4211-8b46-e5bee4397aa8
# Esta función no puede ser ejecutada con vectores p que tengan 0 o 1 en alguna coordenada (tampoco tendría sentido, pues en esos casos, simplemente se elimina de la repartición al partido que tenga coordenada 0 y se le asigna una banca al partido que tenga proba 1)

function selection_probability(A, p)
    N = collect(1:length(p))
    k = length(A)
    permutaciones = collect(permutations(A))
    total = 0.0

    for perm in permutaciones
        prod = 1.0

        for (idx, j) in enumerate(perm)
            N_j = setdiff(N, perm[1:idx-1])
            sum_2 = sum(p[N_j])
            
            q_j = zeros(length(p))
            q_j[N_j] .= p[N_j] ./ sum_2

            z_j = 0.0
            for i in N_j
                z_j += q_j[i] * (1 - q_j[i]) / (1 - (k - idx + 1) * q_j[i])
            end

            prod *= q_j[j] * (1 - q_j[j]) / (z_j * (1 - (k - idx + 1) * q_j[j]))
        end

        total += prod
    end

    return total
end

# ╔═╡ f6294aad-4a04-4fea-aabf-cf9561d30b7e
md"""
    generate_hyperplane_grid(n_dimensions::Int, target_sum::Float64, points_per_dim::Int)

Generates a grid of points for testing apportionment methods.

The function finds points `p = (p₁, p₂, ..., pₙ)` that satisfy three conditions:
1.  The vector has `n_dimensions` coordinates.
2.  Each coordinate `pᵢ` is strictly between 0 and 1 (i.e., `0 < pᵢ < 1`).
3.  The coordinates sum to a specific value: `sum(p) = target_sum`.

### Arguments
- `n_dimensions::Int`: The number of coordinates in each vector (`n`). Must be >= 1.
- `target_sum::Float64`: The value `k` that the coordinates must sum to. For solutions to exist, `0 < target_sum < n_dimensions`.
- `points_per_dim::Int`: The number of grid points to generate along each of the first `n-1` independent dimensions. A larger number creates a denser grid.

### Returns
- `Vector{Vector{Float64}}`: A vector of vectors, where each inner vector is a point on the grid that satisfies all conditions. Returns an empty vector if no such points are found on the grid or if the inputs are invalid.

### Methodology
The problem has `n-1` degrees of freedom. We generate a grid for the first `n-1` coordinates and then calculate the `n`-th coordinate as `pₙ = target_sum - sum(p₁, ..., pₙ₋₁)`. We then filter out any points where any coordinate (including the calculated `pₙ`) falls outside the `(0, 1)` range.

A small epsilon is used to ensure coordinates are *strictly* greater than 0 and less than 1.
"""

# ╔═╡ f1ffd09a-da9c-4c4a-8997-28747f92f210
function generate_hyperplane_grid(n_dimensions::Int, target_sum::Float64, points_per_dim::Int)
    # --- Input Validation ---
    if n_dimensions < 1
        @error "Number of dimensions must be at least 1."
        return Vector{Float64}[]
    end
    if !(0 < target_sum < n_dimensions)
        # If the sum is outside this range, no solutions are possible.
        @warn "target_sum is outside the valid range (0, n_dimensions). No solutions exist."
        return Vector{Float64}[]
    end

    # A small value to ensure strict inequality (0 < p_i < 1)
    epsilon = 1e-9

    # This will hold the valid grid points
    grid_points = Vector{Float64}[]

    # Handle the simple n=1 case separately for clarity
    if n_dimensions == 1
        if 0 < target_sum < 1
            push!(grid_points, [target_sum])
        end
        return grid_points
    end

    # --- Grid Generation ---

    # 1. Create the grid ticks for one dimension
    # We create a grid for each of the first n-1 dimensions
    ticks = range(epsilon, 1.0 - epsilon, length=points_per_dim)

    # 2. Create a Cartesian product iterator for the n-1 dimensions
    # ntuple creates a tuple of `ticks` of length `n_dimensions - 1`
    # The splat `...` operator passes them as separate arguments to product()
    prefix_iterator = Iterators.product(ntuple(i -> ticks, n_dimensions - 1)...)

    # 3. Iterate, calculate the last component, and filter
    for p_prefix in prefix_iterator
        # p_prefix is a tuple, e.g., (0.1, 0.25)
        prefix_sum = sum(p_prefix)
        last_component = target_sum - prefix_sum

        # 4. Check if the generated point is valid
        # We already know the prefix components are in (0,1).
        # We only need to check the last one.
        if epsilon < last_component < 1.0 - epsilon
            # Also, ensure the prefix sum itself doesn't make a valid solution impossible.
            # This check is implicitly handled by the check on last_component, but it's
            # good to be aware of the constraint: target_sum - 1 < prefix_sum < target_sum
            
            # Construct the full valid vector and add it to our list
            full_vector = [p_prefix..., last_component]
            push!(grid_points, full_vector)
        end
    end

    return grid_points
end

# ╔═╡ e57e9672-a42c-42e9-8854-2bdece4f4475
begin
	encontrado = false
	intentos = 0
	points_per_dim = 10
	p = zeros(n)
	p_prima = zeros(n)
	grilla = generate_hyperplane_grid(n, float(k), points_per_dim)
	A = collect(1:k)

	for vect in grilla
		if encontrado
			break
		end
		intentos += 1
		println("Intentos: ", intentos)
		
		# Defino p y p'
		p = vect
		p_prima = perturb_vector(p, 0.1)

		# Calculo P(S=A) para p y p', con A={1,...,k}
		proba_p = selection_probability(A, p)
		proba_p_prima = selection_probability(A, p_prima)

		encontrado = proba_p > proba_p_prima
	end

	if encontrado
		println("Contraejemplo encontrado:")
		println("Vector p:", p)
		println("Vector p':", p_prima)
	else
		println("Grilla de n = ", n, " y k = ", k, " completada. Sin contraejemplos encontrados.")
	end
end

# ╔═╡ 64a3499b-5e78-442f-9178-7277fe88531c
begin
	tuplas = []

	for i in 3:15
		for j in 2:i-1
			tupla = (i, j)
			push!(tuplas, tupla)
		end
	end

	println(tuplas)
end

# ╔═╡ d133bac7-271b-44e8-8350-7515eb05e66a
begin
	for tupla in tuplas
		encontrado = false
		intentos = 0
		n = tupla[1]
		k = tupla[2]
		points_per_dim = 10
		tupla
		p = zeros(n)
		p_prima = zeros(n)
		grilla = generate_hyperplane_grid(n, float(k), points_per_dim)
		A = collect(1:k)
	
		for vect in grilla
			if encontrado
				break
			end
			intentos += 1
			
			# Defino p y p'
			p = vect
			p_prima = perturb_vector(p, 0.1)
	
			# Calculo P(S=A) para p y p', con A={1,...,k}
			proba_p = selection_probability(A, p)
			proba_p_prima = selection_probability(A, p_prima)
	
			encontrado = proba_p > proba_p_prima
		end
	
		if encontrado
			println("Contraejemplo encontrado:")
			println("Vector p:", p)
			println("Vector p':", p_prima)
		else
			println("Grilla de n = ", n, " y k = ", k, " completada. Sin contraejemplos encontrados.")
		end		
	end
end

# ╔═╡ 4123a03b-b86b-4427-b59a-a65bd9249a4d
# ╠═╡ disabled = true
#=╠═╡
begin
	encontrado = false
	intentos = 0
	p = zeros(n)
	p_prima = zeros(n)
	A = collect(1:k)

	while !encontrado
		intentos += 1
		println("Intentos: ", intentos)
		
		# Defino p y p'
		p = rand_vector_sum(n, float(k))
		p_prima = perturb_vector(p, 0.1)

		# Calculo P(S=A) para p y p', con A={1,...,k}
		proba_p = selection_probability(A, p)
		proba_p_prima = selection_probability(A, p_prima)

		encontrado = proba_p > proba_p_prima
	end

	println("Contraejemplo encontrado:")
	println("Vector p:", p)
	println("Vector p':", p_prima)
end
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Combinatorics = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
Combinatorics = "~1.0.2"
StatsBase = "~0.34.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.2"
manifest_format = "2.0"
project_hash = "4e1107f5e957a2938b6cf7faf3d93357feddff18"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

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

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

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

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OrderedCollections]]
git-tree-sha1 = "12f1439c4f986bb868acda6ea33ebc78e19b95ad"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.7.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

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

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"
"""

# ╔═╡ Cell order:
# ╟─040604b0-3b6e-11f0-35ff-c10fb2433828
# ╠═eb417b8f-8824-4263-bcf3-7e722a2291f6
# ╠═57af78d6-c62c-41a4-90dd-6f66cd91122a
# ╠═a29751ad-5eab-4988-ae35-0c6b098014ba
# ╠═1d7b2b4a-2a8b-40bd-9048-5d7f3224f63b
# ╟─1c84a304-7e2a-4038-abbb-f06c0585754d
# ╠═bbeb909b-c55e-4232-ba91-d5e8f3a06266
# ╠═1c8f6a7d-3308-473f-86ec-ab3189d7b77a
# ╟─28b2efaa-bd96-429b-8a38-042b71f7e9a8
# ╠═140ec3dc-8169-4211-8b46-e5bee4397aa8
# ╟─f6294aad-4a04-4fea-aabf-cf9561d30b7e
# ╠═f1ffd09a-da9c-4c4a-8997-28747f92f210
# ╠═e57e9672-a42c-42e9-8854-2bdece4f4475
# ╠═64a3499b-5e78-442f-9178-7277fe88531c
# ╠═d133bac7-271b-44e8-8350-7515eb05e66a
# ╠═4123a03b-b86b-4427-b59a-a65bd9249a4d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
