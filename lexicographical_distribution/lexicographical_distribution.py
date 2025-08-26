# Imports
import itertools
import gurobipy as gp
import math
import numpy as np
from gurobipy import GRB


# Aux functions
def kSubsets(arr, k):
    return list(itertools.combinations(arr, k))


def setsContaining(all_sets):
    """
    Returns a dictionary where each key is an element and the value is a list of sets containing that element.
    """
    sets_containing = {i: [] for i in range(N)}
    for j, S in enumerate(all_sets):
        for i in S:
            sets_containing[i].append(j)
    return sets_containing


def lexicographical_distribution(p: list):
    """
    Generates a distribution over all k-subsets of [n] that gives priority to the subsets that appear before 
    in lexicographical order.
    """
    assert len(p) == N, "Length of p must match N"
    
    # Model initialization
    model = gp.Model("lexicographical_distribution_finder")
    model.Params.TimeLimit = 100            # (seconds) Give Gurobi a time budget
    model.Params.SolutionLimit = 6          # Stop after second feasible solution
    model.Params.OutputFlag = 1             # Show solver output
    model.Params.FeasibilityTol = 1e-09     # Tolerance for feasibility

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
    
    # Solve the problem
    model.optimize()

    n_sols = model.SolCount
    print(f"Optimization ended with status {model.Status} and found {n_sols} solutions.")

    # Retrieve results
    if model.Status == GRB.OPTIMAL or model.Status == GRB.SUBOPTIMAL:
        mu_sol = {i: mu[i].X for i in range(M)}
        return mu_sol
        
    else:
        print(f"Optimization ended with status {model.Status}")
        return None


def results_print(mu_sol):
    """
    Prints the results of the optimization.
    """
    print("-" * 40)
    if mu_sol is None:
        print("No solution found.")
        return

    # Print the optimal distribution
    print("Optimal μ:")    
    for j in range(M):
        print(f"μ({all_sets[j]}) = {mu_sol[j]:.6f}")

    # Print marginal probabilities
    print("-" * 40)
    print("p vs marginal probabilities:")
    for i in range(N):
        print(f"p({i}) = {p[i]:.6f} vs P({i} ∈ S) = {sum(mu_sol[j] for j in sets_containing[i]):.6f}")
    
    # Print the total probability
    print("-" * 40)
    total_prob = sum(mu_sol.values())
    print(f"Total Probability: {total_prob:.4f}")


# Global parameters
N = 5
K = 2
M = math.comb(N, K)
all_sets = kSubsets(range(N), K)
sets_containing = setsContaining(all_sets)

if __name__ == "__main__":
    
    # Define the marginal probabilities vector
    p = [0.5, 0.5, 0.4, 0.35, 0.25]
    p = sorted(p, reverse=True)  # Sort probabilities in descending order

    # Run the lexicographical distribution optimization
    mu_sol = lexicographical_distribution(p)
    
    # Print the results
    results_print(mu_sol)