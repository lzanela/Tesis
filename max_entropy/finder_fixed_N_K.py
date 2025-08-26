from Tesis.max_entropy.max_entropy_distribution import max_entropy, max_entropy_fraction
from fractions import Fraction
import numpy as np
import random
import itertools
import math


N = 6
K = 3
MAX_TRIALS = 1000000

def prod(s, q):
    return math.prod([q[i] for i in s])


def p(q, i):
    return sum([q[s] for s in all_sets if i in s])


def rand_vector_sum_k_fraction(dim=6, total=Fraction(3)):
    while True:
        # Generate dim - 1 random Fractions between 0 and total
        cuts = sorted(Fraction(random.uniform(0, float(total))) for _ in range(dim - 1))
        parts = [cuts[0]] + \
                [cuts[i] - cuts[i - 1] for i in range(1, dim - 1)] + \
                [total - cuts[-1]]
        if all(x < 1 for x in parts):
            return parts


def perturb_vector_fraction(p, delta=Fraction(1, 10), k=3):
    d = len(p)
    assert k < d, "k must be less than the length of the vector"
    admissible = False

    while not admissible:
        p_aux = p.copy()
        shift_total = delta * k
        shift_per_rest = shift_total / (d - k)

        increase_indices = list(range(k))
        decrease_indices = list(range(k, d))

        for i in increase_indices:
            p_aux[i] += delta
        for i in decrease_indices:
            p_aux[i] -= shift_per_rest

        if all(x >= 0 and x < 1 for x in p_aux) and sum(p_aux) == sum(p):
            return p_aux
        else:
            delta *= Fraction(95, 100)


if __name__ == "__main__":
    
    all_sets = list(itertools.combinations(range(N), K))
    contraejemplo_encontrado = False
    p1 = [0] * N
    p2 = [0] * N
    trials = 0

    while not contraejemplo_encontrado and trials < MAX_TRIALS:
        trials += 1
        print(f"Trial {trials} of {MAX_TRIALS}")

        # Generate a random vector p1 that sums to K and its perturbed version p2
        p1 = rand_vector_sum_k_fraction(N, Fraction(K))
        p2 = perturb_vector_fraction(p1, delta=Fraction(1, 10), k=K)

        # Calculate max entropy distributions for p1 and p2
        q1 = max_entropy({i: p1[i] for i in range(N)})
        q2 = max_entropy({i: p2[i] for i in range(N)})

        # Calculate marginal probabilities for p1 and p2 to use as the actual probabilities
        p1_marg = [p(q1, i) for i in range(N)]
        p2_marg = [p(q2, i) for i in range(N)]

        # Check if the generated p1 and p2 are valid
        if not all(0 <= x < 1 for x in p1_marg) or not all(0 <= x < 1 for x in p2_marg):
            print("Invalid p1 or p2, retrying...")
            continue

        # Check if the generated p1 and p2 are valid distributions
        # print(f"p1 sums: {sum(p1_marg)}, p2 sums: {sum(p2_marg)}")
        
        # Calculate selection probabilities for {1,...,K} (actually {0,...,K-1} in Python) with respect to q1 and q2
        q1_selection_prob = q1[(0, 1, 2)]
        q2_selection_prob = q2[(0, 1, 2)]
        
        # Check if the conditions for the counterexample are met
        if q2_selection_prob < q1_selection_prob \
            and p2[0] >= p1[0] and p2[1] >= p1[1] and p2[2] >= p1[2] \
            and p2[3] <= p1[3] and p2[4] <= p1[4] and p2[5] <= p1[5]:

            print(f"Found!")
            print("p1 as floats:")
            print(list(map(float, p1)))

            print("p2 as floats:")
            print(list(map(float, p2)), "\n\n")

            print("p1 as fracs:")
            print("\n".join(map(str, p1)))

            print("p2 as fracs:")
            print("\n".join(map(str, p2)), "\n\n")
            
            print(f"p1 sums: {sum(p1)}, p2 sums: {sum(p2)}", "\n\n")

            print("Selection probabilities:")
            print(f"q1_selection_prob: {float(q1_selection_prob)}")
            print(f"q2_selection_prob: {float(q2_selection_prob)}", "\n\n")