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

def build_bases_with_A1(A: np.ndarray):
    """
    Builds the bases of the matrix A, guaranteeing that the first column of A is in the base.
    """
    m, n = A.shape
    bases = []
    for cols in itertools.combinations(range(1, n), m-1):
        submatrix = A[:, (0,) + cols]
        if np.linalg.matrix_rank(submatrix) == m:
            print(submatrix)
            bases.append((0,) + cols)
    return bases
    
def buscar_contraejemplo(B: np.array):
    """
    Looks for a probability vector p(eps) and p'(eps, alpha, beta, gamma, delta) such that the base B makes the optimal distribution non-selection monotone for p and p'.
    """
    
    # Model initialization
    model = gp.Model("contraejemplo_finder")
    model.Params.TimeLimit = 100            # (seconds) Give Gurobi a time budget
    model.Params.SolutionLimit = 2          # Stop after second feasible solution
    model.Params.OutputFlag = 1             # Show solver output
    model.Params.FeasibilityTol = 1e-09     # Tolerance for feasibility
    model.Params.PoolSolutions = 4  # Store up to 4 solutions
    model.Params.PoolSearchMode = 2 # Find the n best solutions, where n is given by PoolSolutions

    # Variables eps, alpha_1,...,alpha_4 determining the marginal probabilities vector.
    eps = model.addVar(lb=tol, ub=1/5, name="eps")
    alpha = model.addVars(N, lb=0.0, ub=1, name="alpha")

    # Definition of p(eps) and p'(eps, alpha)
    p = [1 - eps, (1+eps)/4, (1+eps)/4, (1+eps)/4, (1+eps)/4]
    p_prime = [1 - eps + alpha[0], (1+eps)/4 + alpha[1], (1+eps)/4 - alpha[2], (1+eps)/4 - alpha[3], (1+eps)/4 - alpha[4]]

    # Valid range for p'
    for i in range(N):
        model.addConstr(alpha[i] <= eps-tol, name=f"p_prime_in_range_{i}")
    
    # Feasibility of B for p and p'
    B_inv = np.linalg.inv(B)
    for i in range(N):
        model.addConstr(gp.quicksum(B_inv[i, j] * p[j] for j in range(N)) >= 0, 
                        name=f"feasibility_p_row_{i}")
        model.addConstr(gp.quicksum(B_inv[i, j] * p_prime[j] for j in range(N)) >= 0, 
                        name=f"feasibility_p_prime_row_{i}")
            
    # Selection monotonicity violation
    model.addConstr(gp.quicksum(B_inv[0,j] * alpha[j] for j in range(N)) <= -tol, name="selection_monotonicity_violation")

    # Objective function
    model.setObjective(gp.quicksum(alpha[i] for i in range(N)) + eps, GRB.MAXIMIZE)
    
    # Solve the problem
    model.write("contraejemplo_finder.lp")
    model.optimize()

    n_sols = model.SolCount
    print(f"Optimization ended with status {model.Status} and found {n_sols} solutions.")

    # Retrieve results
    if model.Status == GRB.OPTIMAL or model.Status == GRB.SUBOPTIMAL:
        eps_sol = eps.X
        alpha_sol = {i: alpha[i].X for i in range(N)}
        return eps_sol, alpha_sol
        
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
    A = build_matrix(all_sets)
    bases = build_bases_with_A1(A)
    bases_satisfying = []
    bases_not_satisfying = []

    print(f"\nTotal bases: {len(bases)}")

    for base in bases:
        B = A[:, base]
        eps_sol, alpha_sol = buscar_contraejemplo(B)
        if eps_sol is not None:
            bases_satisfying.append((base, (eps_sol, alpha_sol)))
            print(f"Base {base} satisfies the condition with (eps, alpha) = {(eps_sol, alpha_sol)}")
        else:
            bases_not_satisfying.append(base)
            print(f"Base {base} does not satisfy the condition.")