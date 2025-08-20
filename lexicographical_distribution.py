# Imports
import itertools
import gurobipy as gp
import math
import numpy as np
from gurobipy import GRB


# Aux functions
def kSubsets(arr, k):
    return list(itertools.combinations(arr, k))


# Global parameters
N = 6
K = 3
MAX_TRIALS = 1000000
eps = 1e-7
all_sets = kSubsets(range(N), K)
M = math.comb(N, K)
p = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
p = sorted(p, reverse=True)  # Sort probabilities in descending order
sets_containing = {i: [] for i in range(N)}

# Define, for each element, the sets it belongs to
for j, S in enumerate(all_sets):
    for i in S:
        sets_containing[i].append(j)


# Model initialization
model = gp.Model("lexicographical_distribution_finder")
model.Params.TimeLimit = 100          # (seconds) Give Gurobi a time budget
model.Params.SolutionLimit = 6       # Stop after second feasible solution
model.Params.OutputFlag = 1  # Show solver output
model.Params.FeasibilityTol = 1e-09

# Variables mu determining the distribution over all sets of K elements
mu = model.addVars(M, lb=0.0, ub=1.0, name="mu")

# Marginal constraints:
for i in range(N):
    model.addConstr(gp.quicksum(mu[j] for j in sets_containing[i]) == p[i], 
                    name=f"marginal_{i}")

# Probability distribution constraint
model.addConstr(gp.quicksum(mu) == 1.0, name="total_probability")

# Objective function
model.setObjective(gp.quicksum(mu[j] * 2**(M-j) for j in range(M)), GRB.MAXIMIZE)

if __name__ == "__main__":
    # Solve the problem
    model.optimize()

    # Retrieve results
    if model.Status == GRB.OPTIMAL or model.Status == GRB.SUBOPTIMAL:
        mu_sol = {i: mu[i].X for i in range(M)}
        
        print("-" * 40)
        print("Optimal μ:")    
        for j in range(M):
            print(f"μ({all_sets[j]}) = {mu_sol[j]:.6f}")

        print("-" * 40)
        print("p vs marginal probabilities:")
        for i in range(N):
            print(f"p({i}) = {p[i]:.6f} vs P({i} ∈ S) = {sum(mu_sol[j] for j in sets_containing[i]):.6f}")
    else:
        print(f"Optimization ended with status {model.Status}")

    n_sols = model.SolCount
    print(" ")
    print(f"Found {n_sols} feasible solutions")