# Imports
import itertools
import math
import numpy as np
import gurobipy as gp
from gurobipy import GRB

# Global parameters
N = 5
K = 2
M = math.comb(N, K)
tol = 5e-6
delta = 0.1

# Aux functions
def indicator_vector(S):
    """
    Returns the indicator vector of the set S.
    """
    vec = np.zeros(N)
    for i in S:
        vec[i] = 1
    return vec

def build_matrix(sets):
    """
    Builds the matrix A where each column is the indicator vector of a set in sets.
    """
    indicator_vectors = [indicator_vector(S) for S in sets]
    A = np.column_stack(indicator_vectors)
    return A

def build_bases(A: np.ndarray) -> list:
    """
    Builds the bases of the matrix A. 
    Returns a list of tuples, each tuple containing the column indices of a base.
    """
    m, n = A.shape
    bases = []
    for cols in itertools.combinations(range(0, n), m):
        submatrix = A[:, cols]
        if np.linalg.matrix_rank(submatrix) == m:
            bases.append(cols)
    return bases
    
def chequear_contraejemplo(B, B_, indice_seleccionado, modo='interseccion'):
    """
    Busca vectores de probabilidades p y p'(p, alpha) tales que p y p' sean factibles
    para las bases B y B' respectivamente, y que a su vez provoquen que las distribuciones 
    óptimas no sean selection monotone.
    """
    
    # Model initialization
    model = gp.Model("contraejemplo_finder")
    model.Params.TimeLimit = 100            # (seconds) Give Gurobi a time budget
    model.Params.SolutionLimit = 2          # Stop after second feasible solution
    model.Params.OutputFlag = 1             # Show solver output
    model.Params.FeasibilityTol = 1e-09     # Tolerance for feasibility
    model.Params.PoolSolutions = 4          # Store up to 4 solutions
    model.Params.PoolSearchMode = 2         # Find the n best solutions, where n is given by PoolSolutions
    model.Params.NumericFocus = 3           # Focus on numerical accuracy

    # Define the matrices and their inverses
    B_matrix = A[:, B]
    B_matrix_ = A[:, B_]
    B_inv = np.linalg.inv(B_matrix)
    B_inv_ = np.linalg.inv(B_matrix_)

    # Determine the index of the selected set in B and B'
    indice_seleccionado_en_B = B.index(indice_seleccionado)
    indice_seleccionado_en_B_ = B_.index(indice_seleccionado)
    
    # Variables alpha_1,...,alpha_N determining the marginal probabilities vector.
    alpha = model.addVars(N, lb=0.0, ub=1, name="alpha")
    p = model.addVars(N, lb=delta, ub=1-delta, name="p")

    # Conjuntos de indices seleccionados para verificar selection monotonicity
    selected_set = set(all_sets[indice_seleccionado])
    unselected_set = set(range(N)).difference(selected_set)

    ## Definition of p(eps) and p'(eps, alpha)
    #p_prime = list(p) + [alpha[0], -alpha[1], -alpha[2], alpha[3]]

    # Valid range for p'
    for i in selected_set:
        model.addConstr(p[i] + alpha[i] <= 1, name=f"p_prime_{i}_leq_1")
        model.addConstr(p[i] + alpha[i] >= 0, name=f"p_prime_{i}_geq_0")
    for i in unselected_set:
        model.addConstr(p[i] - alpha[i] <= 1, name=f"p_prime_{i}_leq_1")
        model.addConstr(p[i] - alpha[i] >= 0, name=f"p_prime_{i}_geq_0")
    
    # p suma K
    model.addConstr(gp.quicksum(p[i] for i in range(N)) == K, name="p_suma_K")

    # p' suma 2 (alpha suma 0)
    model.addConstr(gp.quicksum(alpha[i] for i in selected_set) - 
                    gp.quicksum(alpha[i] for i in unselected_set) == 0, name="alpha_suma_0")

    # Feasibility of B and B' for p and p' resp
    for i in range(N):
        model.addConstr(gp.quicksum(B_inv[i, j] * p[j] for j in range(N)) >= 0, 
                        name=f"feasibility_p_geq_0_row_{i}")
        model.addConstr(gp.quicksum(B_inv[i, j] * p[j] for j in range(N)) <= 1, 
                        name=f"feasibility_p_leq_1_row_{i}")
        model.addConstr(gp.quicksum(B_inv_[i, j] * (p[j] + alpha[j]) for j in selected_set) +
                        gp.quicksum(B_inv_[i, j] * (p[j] - alpha[j]) for j in unselected_set) >= 0, 
                        name=f"feasibility_p_prime_geq_0_row_{i}")
        model.addConstr(gp.quicksum(B_inv_[i, j] * (p[j] + alpha[j]) for j in selected_set) +
                        gp.quicksum(B_inv_[i, j] * (p[j] - alpha[j]) for j in unselected_set) <= 1, 
                        name=f"feasibility_p_prime_leq_1_row_{i}")
        
    # Unfeasibility of all other bases for p and p', introducing auxiliary variables and constraints
    for i in range(len(bases)):
        base_aux = bases[i]
        if base_aux == B or base_aux == B_:
            continue
        base_aux_matrix = A[:, base_aux]
        base_aux_inv = np.linalg.inv(base_aux_matrix)
        
        # Auxiliary variables for infeasibility
        z_p = model.addVars(N, vtype=GRB.BINARY, name=f"z_p_base_{i}")
        z_p_ = model.addVars(N, vtype=GRB.BINARY, name=f"z_p_prime_base_{i}")

        # At least one constraint violated for p
        for row in range(N):
            model.addConstr((gp.quicksum(base_aux_inv[row, j] * p[j] for j in range(N)) + tol) * z_p[row] <= 0,
                            name=f"infeasibility_p_base_{i}_row_{row}_leq_0")
        model.addConstr(gp.quicksum(z_p[row] for row in range(N)) == 1, 
                        name=f"infeasibility_p_base_{i}_at_least_one_violated")
        
        # At least one constraint violated for p'
        for row in range(N):
            model.addConstr((gp.quicksum(base_aux_inv[row, j] * (p[j] + (alpha[j] if j in selected_set else -alpha[j])) 
                                        for j in range(N)) + tol)*z_p_[row] <= 0,
                            name=f"infeasibility_p_prime_base_{i}_row_{row}_leq_0")
        model.addConstr(gp.quicksum(z_p_[row] for row in range(N)) == 1, 
                        name=f"infeasibility_p_prime_base_{i}_at_least_one_violated")
        
    # Unfeasibility of B for p' and of B' for p
    # Auxiliary variables for infeasibility
    z_p_B_prime = model.addVars(N, vtype=GRB.BINARY, name=f"z_p_base_B_prime")
    z_p_prime_B = model.addVars(N, vtype=GRB.BINARY, name=f"z_p_prime_base_B")

    # At least one constraint violated for p
    for row in range(N):
        model.addConstr((gp.quicksum(B_inv_[row, j] * p[j] for j in range(N)) + tol) * z_p_B_prime[row] <= 0,
                        name=f"infeasibility_p_base_B_prime_row_{row}_leq_0")
    model.addConstr(gp.quicksum(z_p_B_prime[row] for row in range(N)) >= 1, 
                    name=f"infeasibility_p_base_B_prime_at_least_one_violated")
    
    # At least one constraint violated for p'
    for row in range(N):
        model.addConstr((gp.quicksum(B_inv[row, j] * (p[j] + (alpha[j] if j in selected_set else -alpha[j])) 
                                    for j in range(N)) + tol)*z_p_prime_B[row] <= 0,
                        name=f"infeasibility_p_prime_base_B_row_{row}_leq_0")
    model.addConstr(gp.quicksum(z_p_prime_B[row] for row in range(N)) >= 1, 
                    name=f"infeasibility_p_prime_base_B_at_least_one_violated")
    
    # Selection monotonicity violation
    if modo=='interseccion':
        model.addConstr(gp.quicksum(B_inv_[indice_seleccionado_en_B_, j] * (p[j] + alpha[j]) for j in selected_set) +
                        gp.quicksum(B_inv_[indice_seleccionado_en_B_, j] * (p[j] - alpha[j]) for j in unselected_set) <= 
                        gp.quicksum(B_inv[indice_seleccionado_en_B, j] * p[j] for j in range(N)) - tol,
                        name="selection_monotonicity_violation")
    else:
        model.addConstr(gp.quicksum(B_inv[indice_seleccionado_en_B, j] * p[j] for j in range(N)) >= tol,
                        name="selection_monotonicity_violation")

    # Objective function
    model.setObjective(gp.quicksum(alpha[j] for j in selected_set), GRB.MAXIMIZE)
    
    # Solve the problem
    model.write("bases_contraejemplo_finder.lp")
    model.optimize()

    n_sols = model.SolCount
    print(f"Optimization ended with status {model.Status} and found {n_sols} solutions.")

    # Retrieve results
    if model.Status == GRB.OPTIMAL or model.Status == GRB.SUBOPTIMAL:
        p_sol = {i: p[i].X for i in range(N)}
        alpha_sol = {i: alpha[i].X for i in range(N)}
        return p_sol, alpha_sol
        
    if model.Status == GRB.INFEASIBLE:
        print("Model is infeasible. Computing IIS...")
        model.computeIIS()
        print("Constraints in IIS:")
        for c in model.getConstrs():
            if c.IISConstr:
                print(f"  {c.ConstrName}: {model.getRow(c)} {c.Sense} {c.RHS}")

    else:
        print(f"Optimization ended with status {model.Status}")
        return None

if __name__ == "__main__":
    # Global setup
    all_sets = list(itertools.combinations(range(N), K))
    print(all_sets)
    
    A = build_matrix(all_sets)
    bases = build_bases(A)
    found = False
    N_SOLUTIONS_DESIRED = 1
    solutions = []

    for _, B in enumerate(bases):
        if len(solutions) >= N_SOLUTIONS_DESIRED:
            break
        for __, B_ in enumerate(bases):
            if len(solutions) >= N_SOLUTIONS_DESIRED:
                break
            if _ == __:
                continue

            print(f"Checking bases {B} and {B_}")
            
            # Chequeo la intersección de B y B'
            interseccion = set(B).intersection(set(B_))

            for i in interseccion:
                sol = chequear_contraejemplo(B, B_, indice_seleccionado=i, modo='interseccion')
                if sol is not None:
                    found = True
                    p_sol, alpha_sol = sol
                    print(f"Bases {B} y {B_} satisfies the condition with (p, alpha) = {(p_sol, alpha_sol)}, for set {all_sets[i]}")
                    solutions.append((B, B_, p_sol, alpha_sol, all_sets[i]))
                    if len(solutions) >= N_SOLUTIONS_DESIRED:
                        break
            
            for sol in solutions:
                B_found, B__found, p_found, alpha_found, selected_set_found = sol
                print("----- Solution Found -----")
                print(f"B_1 = {B_found}")
                print(f"B_2 = {B__found}")
                print(f"p = {tuple(p_found.values())}")
                print(f"alpha = {tuple(alpha_found.values())}")
                print(f"S = {selected_set_found}")