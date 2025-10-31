# Imports
import itertools
import math
import numpy as np
import gurobipy as gp
from gurobipy import GRB

# Global parameters
N = 4
K = 2
M = math.comb(N, K)
tol = 1e-6

# Aux functions
def indicator_vector(S):
    """
    Returns the indicator vector of the set S.
    """
    vec = np.zeros(N)
    for i in S:
        vec[i] = 1
    return vec

def set_from_indicator(vec):
    """
    Returns the set corresponding to the indicator vector vec.
    """
    return {i for i in range(len(vec)) if vec[i] == 1}

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
    Busca un vector de probabilidades p(eps) y p'(eps, alpha) tal que p y p' sean factibles
    para las bases B y B' respectivamente, haciendo que las distribuciones óptimas no sea
    selection monotone.
    """
    
    # Model initialization
    model = gp.Model("contraejemplo_finder")
    model.Params.TimeLimit = 100            # (seconds) Give Gurobi a time budget
    model.Params.SolutionLimit = 2          # Stop after second feasible solution
    model.Params.OutputFlag = 1             # Show solver output
    model.Params.FeasibilityTol = 1e-09     # Tolerance for feasibility
    model.Params.PoolSolutions = 4  # Store up to 4 solutions
    model.Params.PoolSearchMode = 2 # Find the n best solutions, where n is given by PoolSolutions

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
    p = model.addVars(N, lb=0.0, ub=1, name="p")

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
    
    # p suma 2
    model.addConstr(gp.quicksum(p[i] for i in range(N)) == 2, name="p_suma_2")

    # p' suma 2 (alpha suma 0)
    model.addConstr(gp.quicksum(alpha[i] for i in selected_set) - 
                    gp.quicksum(alpha[i] for i in unselected_set) == 0, name="p_prime_suma_2")

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
            
    # Selection monotonicity violation
    if modo=='interseccion':
        model.addConstr(gp.quicksum(B_inv_[indice_seleccionado_en_B_, j] * (p[j] + alpha[j]) for j in selected_set) +
                        gp.quicksum(B_inv_[indice_seleccionado_en_B_, j] * (p[j] - alpha[j]) for j in unselected_set) <= 
                        gp.quicksum(B_inv[indice_seleccionado_en_B, j] * p[j] for j in selected_set) - tol,
                        name="selection_monotonicity_violation")
    else:
        model.addConstr(gp.quicksum(B_inv[indice_seleccionado_en_B, j] * p[j] for j in selected_set) >= tol,
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
    c = [6,5,4,3,2,1]
    all_sets = list(itertools.combinations(range(N), K))
    A = build_matrix(all_sets)
    bases = build_bases(A)

    for _, B in enumerate(bases):
        for __, B_ in enumerate(bases):
            if _ >= __:
                continue
            print(f"Checking bases {B} and {B_}")
            
            # Chequeo la intersección de B y B'
            interseccion = set(B).intersection(set(B_))

            for i in interseccion:
                sol = chequear_contraejemplo(B, B_, indice_seleccionado=i, modo='interseccion')
                if sol is not None:
                    p_sol, alpha_sol = sol
                    print(f"Bases {B} y {B_} satisfies the condition with (p, alpha) = {(p_sol, alpha_sol)}, for set {all_sets[i]}")

    print("-"*50)
    B_encontrado = A[:, (1, 3, 4, 5)]
    B_encontrado_ = A[:, (2, 3, 4, 5)]
    print("B_inv = ", np.linalg.inv(B_encontrado))
    print("B'_inv = ", np.linalg.inv(B_encontrado_))
    print("-"*50)