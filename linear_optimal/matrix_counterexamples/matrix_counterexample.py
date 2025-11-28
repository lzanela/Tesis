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

def build_bases(A: np.ndarray):
    """
    Builds the bases of the matrix A.
    """
    m, n = A.shape
    bases = []
    for cols in itertools.combinations(range(n), m):
        submatrix = A[:, cols]
        if np.linalg.matrix_rank(submatrix) == m:
            bases.append(cols)
    return bases

def complement(S):
    """
    Returns the complement of the set S with respect to the universal set {0, 1, ..., N-1}.
    """
    return {i for i in range(N) if i not in S}

def verify_condition(B: np.ndarray):
    """
    Checks if the inverse of the basis matrix B_inv satisfies the condition:
    It exists a pair of elements (i, j) such that there's a row l in B_inv such that i in set_l, j not in set_l and B_inv[l, i] - B_inv[l, j] < 0. 
    """
    B_inv = np.linalg.inv(B)
    m, n = B_inv.shape

    for l in range(m):
        set_l = set_from_indicator(B[:, l])
        row = B_inv[l, :]
        for i in set_l:
            for j in complement(set_l):
                if row[i] - row[j] < 0:
                    return True, (l, i, j)
    return False, None
    
def buscar_p(bases_not_satisfying):
    """
    Looks for a probability vector p such that the every basis satisfying the condition has a vector p that makes the optimal distribution non-selection monotone.
    """
    #assert (len(p) == N) and (len(c) == M), "Length of p must match N and length of c must match M"
    
    # Model initialization
    model = gp.Model("p_finder")
    model.Params.TimeLimit = 100            # (seconds) Give Gurobi a time budget
    model.Params.SolutionLimit = 2          # Stop after second feasible solution
    model.Params.OutputFlag = 1             # Show solver output
    model.Params.FeasibilityTol = 1e-09     # Tolerance for feasibility
    model.Params.PoolSolutions = 4  # Store up to 4 solutions
    model.Params.PoolSearchMode = 2 # Find the n best solutions, where n is given by PoolSolutions

    # Variables p determining the marginal probabilities vector.
    p = model.addVars(N, lb=0.0, ub=1.0, name="p")

    # Unfeasibility constraints for each basis not satisfying the condition
    for base in bases_not_satisfying:
        B = A[:, base]
        B_inv = np.linalg.inv(B)
        for i in range(N):
            model.addConstr(gp.quicksum(B_inv[i, j] * p[j] for j in range(N)) <= 0, 
                            name=f"unfeasibility_base_{base}_row_{i}")
            
    # Force p to be different than 0
    model.addConstr(gp.quicksum(p[i] for i in range(N)) >= 1e-3, name="p_nonzero")

    # Objective function
    model.setObjective(gp.quicksum(p[i] for i in range(N)), GRB.MAXIMIZE)
    
    # Solve the problem
    model.optimize()

    n_sols = model.SolCount
    print(f"Optimization ended with status {model.Status} and found {n_sols} solutions.")

    # Retrieve results
    if model.Status == GRB.OPTIMAL or model.Status == GRB.SUBOPTIMAL:
        p_sol = {i: p[i].X for i in range(N)}
        return p_sol
        
    else:
        print(f"Optimization ended with status {model.Status}")
        return None

if __name__ == "__main__":

    # Global setup
    all_sets = list(itertools.combinations(range(N), K))
    A = build_matrix(all_sets)
    bases = build_bases(A)
    bases_satisfying = []
    bases_not_satisfying = []

    for base in bases:
        B = A[:, base]
        condition, indices = verify_condition(B)
        if condition:
            bases_satisfying.append((base, indices))
            print(f"Base {base} satisfies the condition with (l, i, j) = {indices}")
        else:
            bases_not_satisfying.append(base)
            print(f"Base {base} does not satisfy the condition.")
    
    print(f"\nTotal bases: {len(bases)}")
    print(f"Bases not satisfying the condition: {len(bases_not_satisfying)}")
    
    p_sol = buscar_p(bases_not_satisfying)
    print(f"p_sol = {p_sol}")