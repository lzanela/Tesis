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

# Gurobi's solution
B_1 = (0, 1, 2, 3)
B_2 = (0, 1, 2, 4)
p = (0.10000199999625015, 0.8999990000012499, 0.8999990000012498, 0.1)
alpha = (0.8999970000050003, 0.8999970000037498, 0.0, 0.0)
S = (0, 3)
S_index = 2

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
    
def verificar_contraejemplo(B: tuple, B_:tuple, p:tuple, alpha:tuple, S: set, S_index: int, modo='interseccion') -> dict:
    """
    Se asume que los vectores p y p'(p, alpha) están en [0,1]^N, que B y B_ son bases de A representados como tuplas
    con los índices de las columnas seleccionadas,

    Verifica si los vectores de probabilidades p y p'(p, alpha) cumplen las siguientes condiciones: 
    - p y p' son factibles para las bases B y B' respectivamente,
    - p y p' suman K (dentro de una tolerancia),
    - Las distribuciones mu y mu' obtenidas a partir de B^(-1)*p y B_^(-1)*p' son factibles,
    - La selection monotonicity no se cumple entre mu y mu' para el conjunto S.
    """
    # Define the matrices and their inverses
    B_matrix = A[:, B]
    B_matrix_ = A[:, B_]
    B_inv = np.linalg.inv(B_matrix)
    B_inv_ = np.linalg.inv(B_matrix_)

    # Determine the index of the selected set in B and B'
    indice_seleccionado_en_B = B.index(S_index)
    indice_seleccionado_en_B_ = B_.index(S_index)

    # Define a dictionary with the conditions to verify
    conditions = {}
    
    # Definition of p'(p, alpha)
    p_prime = tuple(
        p[i] + alpha[i] if i in S else p[i] - alpha[i]
        for i in range(N)
    )
    
    # p suma 2
    conditions["p_suma_2"] = abs(sum(p[i] for i in range(N)) - 2) <= tol

    # alpha suma 0 (p' suma 2)
    conditions["alpha_suma_0"] = abs(sum(alpha[i] for i in S) - sum(alpha[i] for i in range(N) if i not in S)) <= tol

    # Feasibility of B and B' for p and p' resp
    for i in range(N):
        conditions[f"feasibility_mu_{i}_geq_0"] = sum(B_inv[i, j] * p[j] for j in range(N)) >= 0
        conditions[f"feasibility_mu_{i}_leq_1"] = sum(B_inv[i, j] * p[j] for j in range(N)) <= 1
        conditions[f"feasibility_mu_prime_{i}_geq_0"] = sum(B_inv_[i, j] * p_prime[j] for j in range(N)) >= 0
        conditions[f"feasibility_mu_prime_{i}_leq_1"] = sum(B_inv_[i, j] * p_prime[j] for j in range(N)) <= 1
    
    # Unfeasibility of all other bases for p and p'
    for i in range(len(bases)):
        base_aux = bases[i]
        if base_aux == B or base_aux == B_:
            continue
        base_aux_matrix = A[:, base_aux]
        base_aux_inv = np.linalg.inv(base_aux_matrix)
        
        conditions[f"infeasibility_p_base_{i}"] = False

        # At least one constraint violated for p
        for row in range(N):
            if (sum(base_aux_inv[row, j] * p[j] for j in range(N)) < 0):
                conditions[f"infeasibility_p_base_{i}"] = True
                break
        
        # At least one constraint violated for p'
        for row in range(N):
            if (sum(base_aux_inv[row, j] * p_prime[j] for j in range(N)) < 0):
                conditions[f"infeasibility_p_prime_base_{i}"] = True
                break

    # Unfeasibility of B for p' and of B' for p
    conditions[f"infeasibility_p_base_B_prime"] = False
    conditions[f"infeasibility_p_prime_base_B"] = False
    # At least one constraint violated for p
    for row in range(N):
            if (sum(B_inv_[row, j] * p[j] for j in range(N)) < 0):
                conditions[f"infeasibility_p_base_B_prime"] = True
                break
    # At least one constraint violated for p'
    for row in range(N):
        if (sum(B_inv[row, j] * p_prime[j] for j in range(N)) < 0):
            conditions[f"infeasibility_p_prime_base_{i}"] = True
            break
            
    # Selection monotonicity violation
    if modo == 'interseccion':
        conditions["selection_monotonicity_violation"] = sum(B_inv_[indice_seleccionado_en_B_, j] * p_prime[j] for j in range(N)) < \
                                                         sum(B_inv[indice_seleccionado_en_B, j] * p[j] for j in range(N))
    else:
        conditions["selection_monotonicity_violation"] = 0 < sum(B_inv[indice_seleccionado_en_B, j] * p[j] for j in range(N))

    return conditions

if __name__ == "__main__":
    # Global setup
    all_sets = list(itertools.combinations(range(N), K))
    
    A = build_matrix(all_sets)
    bases = build_bases(A)
    conditions = verificar_contraejemplo(B_1, B_2, p, alpha, S, S_index, modo='interseccion')
    all_conditions = all(conditions.values())

    print(f"All conditions met: {all_conditions}")
    for key, value in conditions.items():
        print(f"{key}: {value}")