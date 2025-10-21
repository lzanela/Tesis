# Imports
import matplotlib.pyplot as plt
import gurobipy as gp
from gurobipy import GRB

# Global parameters
K = 6
LEN = 1000
tol = 1e-6

# Aux functions
def plot_results(K, interval_length, p_sol, p_r_sol, p_b_sol, w_r_sol, w_b_sol):
    # Create the plot
    plt.figure(figsize=(10, 6))

    # Plot the line from 0 to LEN
    plt.plot([0, LEN], [0, 0], 'k-', lw=1)

    # Plot vertical lines for p_sol
    for i, p_val in p_sol.items():
        plt.axvline(x=p_val, color='black', linestyle='--', label=f'p_{i+1}={p_val:.2f}')
        plt.text(p_val, 0.02, f'p_{i+1}', rotation=90, verticalalignment='bottom')

    # Plot red points with weights
    for i, p_r_val in p_r_sol.items():
        plt.scatter(p_r_val, 0.1, color='red', label=f'p_r_{i}={p_r_val:.2f}, w_r={w_r_sol[i]:.2f}')
        plt.text(p_r_val, 0.12, f'w_r={w_r_sol[i]:.2f}', color='red', fontsize=8, ha='center')

    # Plot blue points with weights
    for i, p_b_val in p_b_sol.items():
        plt.scatter(p_b_val, -0.1, color='blue', label=f'p_b_{i}={p_b_val:.2f}, w_b={w_b_sol[i]:.2f}')
        plt.text(p_b_val, -0.12, f'w_b={w_b_sol[i]:.2f}', color='blue', fontsize=8, ha='center')

    # Add labels and legend
    plt.title("Optimal Distribution of Red and Blue Points")
    plt.xlabel("Interval")
    plt.ylabel("Weights")
    plt.ylim(-0.2, 0.2)
    plt.grid(True)
    plt.legend(loc='upper right', fontsize='small', frameon=False)

    # Show the plot
    plt.show()

def run_optimization(K: int = K, interval_length: int = LEN):
    """
    Looks for a distribution of red and blue points on the interval [0,interval_length] with the quantiles T_1,...T_k, and a bound v to optimize.
    """
    
    # Model initialization
    model = gp.Model("red_and_blue_points")
    model.Params.TimeLimit = 120            # (seconds) Give Gurobi a time budget
    #model.Params.SolutionLimit = 2         # Stop after second feasible solution
    model.Params.OutputFlag = 1             # Show solver output
    model.Params.FeasibilityTol = 1e-09     # Tolerance for feasibility

    # Variables p_1,...,p_k determining the thresholds
    p = model.addVars(K, lb=0.0, ub=interval_length, name="p")

    # Variables p_i^r and p_i^b determining the number of red and blue points in each interval
    p_r = model.addVars(K+1, lb=0.0, ub=interval_length, name="p_r")
    p_b = model.addVars(K+1, lb=0.0, ub=interval_length, name="p_b")

    # Variables w_i^r and w_i^b determining the weight of red and blue points in each interval
    w_r = model.addVars(K+1, lb=0.0, ub=1.0, name="w_r")
    w_b = model.addVars(K+1, lb=0.0, ub=1.0, name="w_b")

    # Variable v determining the bound to optimize
    v = model.addVar(lb=0.0, ub=1.0, name="v")

    # Thresholds order constraints
    model.addConstr(p[0] >= 0, name="p_1_lower_bound")
    model.addConstr(p[K-1] <= interval_length, name=f"p_{K}_upper_bound")
    for i in range(K-1):
        model.addConstr(p[i] <= p[i+1] - 1, name=f"p_{i+1}_less_equal_p_{i+2}")
    
    # Correct intervals for p_i^r and p_i^b
    model.addConstr(p_r[0] <= p[0], name="p_r_0_correct")
    model.addConstr(p_b[0] <= p[0], name="p_b_0_correct")
    model.addConstr(p[K-1] <= p_r[K], name="p_r_K_correct")
    model.addConstr(p[K-1] <= p_b[K], name="p_b_K_correct")
    for i in range(1, K):
        model.addConstr(p[i-1] <= p_r[i], name=f"p_r_{i}_left_correct")
        model.addConstr(p_r[i] <= p[i], name=f"p_r_{i}_right_correct")
        model.addConstr(p[i-1] <= p_b[i], name=f"p_b_{i}_left_correct")
        model.addConstr(p_b[i] <= p[i], name=f"p_b_{i}_right_correct")

    # Weight definition for red points
    for i in range(K+1):
        model.addConstr(w_r[i] + w_b[i] == 2/(K+1) , name=f"weights_{i}_correct_sum")
    model.addConstr(gp.quicksum(w_r[i] for i in range(K+1)) == 1 , name=f"correct_sum_of_r_weights")
    model.addConstr(gp.quicksum(w_b[i] for i in range(K+1)) == 1, name=f"correct_sum_of_b_weights")

    # Definition of C(p_i) and UB(p_i)
    C = list(range(K))
    UBv = list(range(K))
    y = list(range(K))
    for i in range(K):
        C[i] = gp.quicksum(gp.quicksum((p_b[l] - p_r[j])*(w_b[l]*w_r[j]) + (p_r[l] - p_b[j])*(w_r[l]*w_b[j]) for l in range(i+1, K+1)) for j in range(i+1))
        UBv[i] = gp.quicksum(gp.quicksum(((p[i] - p_r[j]) + (p[i] - p_b[l])) * (w_b[l]*w_r[j]) for l in range(i+1)) for j in range(i+1)) + \
                    gp.quicksum(gp.quicksum(((p_r[j] - p[i]) + (p_b[l] - p[i])) * (w_b[l]*w_r[j]) for l in range(i+1, K+1)) for j in range(i+1, K+1)) * v
        y[i] = model.addVar(name=f"y_{i}")
        model.addConstr(y[i] == UBv[i] - C[i], name=f"UBv_minus_C_{i}")

    # Bound constraints
    for i in range(K):
        model.addConstr(y[i] >= 0, name=f"bound_{i}")
        
    # Objective function
    model.setObjective(v, GRB.MINIMIZE)
    
    # Solve the problem
    model.write("red_and_blue_points.lp")
    model.optimize()

    n_sols = model.SolCount
    print(f"Optimization ended with status {model.Status} and found {n_sols} solutions.")

    # Retrieve results
    if model.Status == GRB.INFEASIBLE:
        print("Model is infeasible. Computing IIS...")
        model.computeIIS()
        print("Constraints in IIS:")
        for c in model.getConstrs():
            if c.IISConstr:
                print(f"  {c.ConstrName}: {model.getRow(c)} {c.Sense} {c.RHS}")

    else:
        v_sol = v.X
        p_sol = {i: p[i].X for i in range(K)}
        p_r_sol = {i: p_r[i].X for i in range(K+1)}
        p_b_sol = {i: p_b[i].X for i in range(K+1)}
        w_r_sol = {i: w_r[i].X for i in range(K+1)}
        w_b_sol = {i: w_b[i].X for i in range(K+1)}
        return v_sol, p_sol, p_r_sol, p_b_sol, w_r_sol, w_b_sol
        

if __name__ == "__main__":

    # Global setup
    v_sol, p_sol, p_r_sol, p_b_sol, w_r_sol, w_b_sol = run_optimization()

    #Plot results
    print(f"\nOptimal value v: {v_sol}")
    plot_results(K, LEN, p_sol, p_r_sol, p_b_sol, w_r_sol, w_b_sol)