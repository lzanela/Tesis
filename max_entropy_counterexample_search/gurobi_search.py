import gurobipy as gp
import math
from gurobipy import GRB
import itertools
from fractions import Fraction

# Auxiliary functions
def prod(s, q):
    return math.prod([q[i] for i in s])

def p(q, i):
    return sum([prod(s, q) for s in all_sets if i in s])

# Problem parameters
N = 10
K = 3
eps = 1e-15
all_sets = list(itertools.combinations(range(N), K))
S = tuple(range(K))
S_c = tuple(range(K, N))

# Jamie's starting point
Q = list(map(Fraction, [99620001435175085845613951348591, 33206667145059699577734936400435, 33206667145059699577734936400435, 23244667001544291253373835102276586, 23244667001544291253373835102276586, 1660333357252963458777541885429371]))
Q = Q + [Fraction(0)] * (N - len(Q))
Q_prod_normalization = (1/sum([prod(s, Q) for s in all_sets]))**(1/K)
Q = [x * Q_prod_normalization for x in Q]
Q_ = list(map(Fraction, [99620001435175193801835755646020, 33206667145059681577227243883092, 33206667145059681577227243883092, 23244667001544299141767505142336500, 23244667001544299141767505142336500, 1660333357252962147206216649823732]))
Q_ = Q_ + [Fraction(0)] * (N - len(Q_))
Q__prod_normalization = (1/sum([prod(s, Q_) for s in all_sets]))**(1/K)
Q_ = [x * Q__prod_normalization for x in Q_]

# Model initialization
model = gp.Model("counterexample_finder")
model.Params.TimeLimit = 100          # (seconds) Give Gurobi a time budget
model.Params.SolutionLimit = 6       # Stop after second feasible solution
model.Params.NonConvex = 2  # Enable non-convex QC constraints
model.Params.OutputFlag = 1  # Show solver output
model.Params.FeasibilityTol = 1e-09

# Variables q and q'
q = model.addVars(N, lb=0.0, ub=1.0, name="q")
q_ = model.addVars(N, lb=0.0, ub=1.0, name="q_")

# Constraints: 0 < q[i] for i in S
for i in S:
    model.addConstr(q[i] >= eps, name=f"q_positive_{i}")
    model.addConstr(q_[i] >= eps, name=f"q_prime_positive_{i}")

# Auxiliary variables for monomials (products of 3 variables)
aux = {}
aux_ = {}
for s in all_sets:
    aux[s] = model.addVar(lb=0.0, ub=1.0, name=f"aux_{s}")
    aux_[s] = model.addVar(lb=0.0, ub=1.0, name=f"aux_{s}_prime")

# Define bilinear product constraints: aux_s = q[i] * q[j] * q[k]
# We'll introduce intermediate variables to break down trilinear into bilinear forms
bilin_aux = {}
bilin_aux_ = {}
for s in all_sets:
    i, j, k = s
    # Intermediate variables for bilinear products
    m1 = model.addVar(lb=0.0, ub=1.0, name=f"m1_{s}")
    m1_ = model.addVar(lb=0.0, ub=1.0, name=f"m1_{s}_prime")
    bilin_aux[s] = m1
    bilin_aux_[s] = m1_

    # m1 = q[i] * q[j]
    model.addConstr(m1 == q[i] * q[j], name=f"bilin_{i}_{j}")
    model.addConstr(m1_ == q_[i] * q_[j], name=f"bilin_{i}_{j}_prime")

    # aux_s = m1 * q[k]
    model.addConstr(aux[s] == m1 * q[k], name=f"trilin_{i}_{j}_{k}")
    model.addConstr(aux_[s] == m1_ * q_[k], name=f"trilin_{i}_{j}_{k}_prime")

# Constraint: prod(S, q_) <= prod(S, q) - eps
# and also prod(S, q), prod(S, q_) >= eps 
model.addConstr(aux_[S] <= aux[S] - eps, name="main_prod_constraint")
model.addConstr(aux[S] >= eps, name="positive_selection_prob")
model.addConstr(aux_[S] >= eps, name="positive_selection_prob_prime")

# Constraints: p(q_, i) >= p(q, i) for i in S
for i in S:
    lhs = gp.quicksum(aux_[s] for s in all_sets if i in s)
    rhs = gp.quicksum(aux[s] for s in all_sets if i in s)
    model.addConstr(lhs >= rhs, name=f"ineq_S_{i}")

# Constraints: p(q_, i) <= p(q, i) for i in S_c
for i in S_c:
    lhs = gp.quicksum(aux_[s] for s in all_sets if i in s)
    rhs = gp.quicksum(aux[s] for s in all_sets if i in s)
    model.addConstr(lhs <= rhs, name=f"ineq_Sc_{i}")

# Normalization constraints: sum of all monomials == 1
model.addConstr(gp.quicksum(aux[s] for s in all_sets) == 1, name="norm_q")
model.addConstr(gp.quicksum(aux_[s] for s in all_sets) == 1, name="norm_q_prime")

# Objective: Minimize sum of squared differences (q[i] - q_[i])^2
#model.setObjective(gp.quicksum((q[i] - q_[i]) * (q[i] - q_[i]) for i in range(N)), GRB.MINIMIZE)
model.setObjective(0, GRB.MAXIMIZE)

# Warm start for all variables
for i in range(N):
    q[i].Start = Q[i]
    q_[i].Start = Q_[i]
for s in all_sets:
    i, j, k = s
    aux[s].Start = Q[i] * Q[j] * Q[k]
    aux_[s].Start = Q_[i] * Q_[j] * Q_[k]
    bilin_aux[s].Start = Q[i] * Q[j]
    bilin_aux_[s].Start = Q_[i] * Q_[j]

model.update()

print("Warm start for q:")
for i in range(N):
    print(f"q[{i}] Start =", q[i].Start)

print("\nWarm start for q_:")
for i in range(N):
    print(f"q_[{i}] Start =", q_[i].Start)

print("\nWarm start for aux:")
for s in all_sets:
    print(f"aux[{s}] Start =", aux[s].Start)

print("\nWarm start for aux_:")
for s in all_sets:
    print(f"aux_[{s}] Start =", aux_[s].Start)

print("\nWarm start for bilin_aux:")
for s in all_sets:
    print(f"bilin_aux[{s}] Start =", bilin_aux[s].Start)

print("\nWarm start for bilin_aux_:")
for s in all_sets:
    print(f"bilin_aux_[{s}] Start =", bilin_aux_[s].Start)

# Solve the problem
model.optimize()

# Retrieve results
if model.Status == GRB.OPTIMAL or model.Status == GRB.SUBOPTIMAL:
    q_sol = {i: q[i].X for i in range(N)}
    q__sol = {i: q_[i].X for i in range(N)}
    print("Optimal q:", q_sol)
    print("Optimal q_:", q__sol)
else:
    print(f"Optimization ended with status {model.Status}")

n_sols = model.SolCount
print(f"Found {n_sols} feasible solutions")

for k in range(n_sols):
    model.setParam(GRB.Param.SolutionNumber, k)
    q_sol_k = [q[i].Xn for i in range(N)]
    q__sol_k = [q_[i].Xn for i in range(N)]
    
    q_normalization = sum([prod(s, q_sol_k) for s in all_sets])
    q__normalization = sum([prod(s, q__sol_k) for s in all_sets])
    p_sol_k = [p(q_sol_k, i) / q_normalization for i in range(N)]
    p__sol_k = [p(q__sol_k, i) / q__normalization for i in range(N)]
    selection_prob_k = prod(S, q_sol_k) / q_normalization
    selection_prob__k = prod(S, q__sol_k) / q__normalization

    print("-" * 40)
    print(f"Solution {k+1}: q = {q_sol_k}")
    print(f"           q' = {q__sol_k} \n")
    print(f"           p(q) = {p_sol_k}")
    print(f"          p(q') = {p__sol_k} \n")
    print(f"Normalization (q): {q_normalization}")
    print(f"Normalization (q'): {q__normalization} \n")
    print(f"Selection prob (q):  {selection_prob_k}")
    print(f"Selection prob (q'): {selection_prob__k} \n")
    print(f"Sum of p(q): {sum(p_sol_k)}")
    print(f"Sum of p(q'): {sum(p__sol_k)} \n")
