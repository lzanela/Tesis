# Imports
import itertools
import gurobipy as gp
import math
import numpy as np
import random
from fractions import Fraction
from gurobipy import GRB



# Global parameters
N = 4
K = 2
M = math.comb(N, K)
c = [random.randint(0, 1000) for _ in range(M)]


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


def rand_vector_sum_k_fraction(dim=4, total=Fraction(2)):
    while True:
        # Generate dim - 1 random Fractions between 0 and total
        cuts = sorted(Fraction(random.uniform(0, float(total))) for _ in range(dim - 1))
        parts = [cuts[0]] + \
                [cuts[i] - cuts[i - 1] for i in range(1, dim - 1)] + \
                [total - cuts[-1]]
        if all(x < 1 for x in parts):
            return parts


def perturb_vector(p: list):
    """
    Perturbs the vector p by increasing the first component by eps and decreasing the last by eps, where eps is a 
    random value between 0 and min{1-p[0], p[n-1]}.
    Ensures that all components remain in [0, 1].
    """
    n = len(p)
    p_ = p.copy()
    delta_max = min(1-p[0], p[n-1])
    delta = random.uniform(0, delta_max)
    eps = min(delta, 0.05)  # Limit the perturbation to 0.05

    p_[0] = p_[0] + eps
    p_[n-1] = p_[n-1] - eps
    
    return p_


def optimal_distribution(p: list, c: list = c):
    """
    Generates a distribution over all k-subsets of [n] that maximizes sum(c[j]*mu[j] for j in range(M)).
    """
    assert (len(p) == N) and (len(c) == M), "Length of p must match N and length of c must match M"
    
    # Model initialization
    model = gp.Model("optimal_distribution_finder")
    model.Params.TimeLimit = 100            # (seconds) Give Gurobi a time budget
    model.Params.SolutionLimit = 2          # Stop after second feasible solution
    model.Params.OutputFlag = 1             # Show solver output
    model.Params.FeasibilityTol = 1e-09     # Tolerance for feasibility

    # Variables mu determining the distribution over all sets of K elements
    mu = model.addVars(M, lb=0.0, name="mu")

    # Marginal constraints:
    for i in range(N):
        model.addConstr(gp.quicksum(mu[j] for j in sets_containing[i]) == p[i], 
                        name=f"marginal_{i}")

    # Objective function
    model.setObjective(gp.quicksum(mu[j] * c[j] for j in range(M)), GRB.MAXIMIZE)
    
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
    

def dual_optimal(p: list, c: list = c):
    """
    Solves the dual of the other problem.
    """
    assert (len(p) == N) and (len(c) == M), "Length of p must match N and length of c must match M"
    
    # Model initialization
    model = gp.Model("optimal_dual_distribution_finder")
    model.Params.TimeLimit = 100            # (seconds) Give Gurobi a time budget
    model.Params.SolutionLimit = 2          # Stop after second feasible solution
    model.Params.OutputFlag = 1             # Show solver output
    model.Params.FeasibilityTol = 1e-09     # Tolerance for feasibility

    # Variables mu determining the distribution over all sets of K elements
    alpha = model.addVars(N, lb=-GRB.INFINITY, name="alpha")

    # Marginal constraints:
    for j in range(M):
        model.addConstr(gp.quicksum(alpha[i] for i in all_sets[j]) >= c[j], 
                        name=f"dual_set_{j}")

    # Objective function
    model.setObjective(gp.quicksum(alpha[i] * p[i] for i in range(N)), GRB.MINIMIZE)
    
    # Solve the problem
    model.optimize()

    n_sols = model.SolCount
    print(f"Optimization ended with status {model.Status} and found {n_sols} solutions.")

    # Retrieve results
    if model.Status == GRB.OPTIMAL or model.Status == GRB.SUBOPTIMAL:
        alpha_sol = {i: alpha[i].X for i in range(N)}
        return alpha_sol
        
    else:
        print(f"Optimization ended with status {model.Status}")
        return None


def results_summary(mu_sol, mu_sol_, p, p_, c):
    """
    Summarizes the results of the optimization.
    """
    if mu_sol is None or mu_sol_ is None:
        print("No solution found for one of the distributions.")
        return

    print("-" * 40)
    print("Summary of Results")
    print("-" * 40)
    print(f"eps = {abs(p_[0]-p[0]):.6f}")
    print(f"c = {c}")

    # Print marginal probabilities
    print("p and p':")
    for i in range(N):
        print(f"p({i}) = {p[i]:.6f}  |  p'({i}) = {p_[i]:.6f}")
    print("-" * 40)

    # Print the optimal distributions
    print("Optimal μ and μ':")    
    for j in range(M):
        print(f"μ({all_sets[j]}) = {mu_sol[j]:.6f}  |  μ'({all_sets[j]}) = {mu_sol_[j]:.6f}")
    print("-" * 40)

    # Print optimal distributions differences
    print("Differences in μ and μ':")    
    for j in range(M):
        print(f"μ({all_sets[j]}) - μ'({all_sets[j]}) = {mu_sol[j]-mu_sol_[j]}")
    print("-" * 40)

    # Print objective values difference
    print("Objective values:")
    print(f"opt = {sum(mu_sol[j] * c[j] for j in range(M))} | opt' {sum(mu_sol_[j] * c[j] for j in range(M))}")
    print("-" * 40)

    # Print selection monotonicity condition for {1,...,K}
    print("Selection monotonicity condition for {0,...,K-1}:")
    cond = mu_sol_[0] >= mu_sol[0]
    if cond:
        print(f"Condition holds: μ'({all_sets[0]}) >= μ({all_sets[0]})")
        print(f"{mu_sol_[0]:.6f} >= {mu_sol[0]:.6f}")
    else:
        print(f"Condition does not hold: μ'({all_sets[0]}) < μ({all_sets[0]})")
        print(f"{mu_sol_[0]} < {mu_sol[0]}")
    
    return (not cond) # Return True if a counterexample is found

def dual_results_summary(alpha_sol, alpha_sol_, p, p_, c):
    """
    Summarizes the results of the dual optimization.
    """
    if alpha_sol is None or alpha_sol_ is None:
        print("No solution found for one of the distributions.")
        return

    print("-" * 40)
    print("Summary of Dual Results")
    print("-" * 40)
    print(f"eps = {abs(p_[0]-p[0]):.6f}")
    print(f"c = {c}")

    # Print the optimal values
    print("Optimal α and α':")    
    for i in range(N):
        print(f"α({i}) = {alpha_sol[i]:.6f}  |  α'({i}) = {alpha_sol_[i]:.6f}")
    print("-" * 40)

    # Print optimal values differences
    print("Differences in α and α':")    
    for i in range(N):
        print(f"α({i}) - α'({i}) = {alpha_sol[i]-alpha_sol_[i]}")
    print("-" * 40)

    # Print objective values difference
    print("Objective values for the dual:")
    print(f"opt = {sum(alpha_sol[i] * p[i] for i in range(N))} | opt' {sum(alpha_sol_[i] * p_[i] for i in range(N))}")
    print("-" * 40)

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


if __name__ == "__main__":
    # Global setup
    random.seed(43)
    all_sets = kSubsets(range(N), K)
    sets_containing = setsContaining(all_sets)
    #c = [random.randint(0, 1000) for _ in range(M)]
    c = [(6-i) for i in range(M)]
    contraejemplo_found = False

    while not contraejemplo_found:
        # Define the marginal probabilities vectors
        p = rand_vector_sum_k_fraction(dim=N, total=Fraction(K))
        p_ = perturb_vector(p)

        # Run the linear optimal distribution optimization
        mu_sol = optimal_distribution(p, c=c)
        mu_sol_ = optimal_distribution(p_, c=c)
        
        # Run the dual optimization
        alpha_sol = dual_optimal(p, c=c)
        alpha_sol_ = dual_optimal(p_, c=c)

        # Print the results
        contraejemplo_found = results_summary(mu_sol, mu_sol_, p, p_, c)
        dual_results_summary(alpha_sol, alpha_sol_, p, p_, c)