from lexicographical_distribution import lexicographical_distribution, results_print
from fractions import Fraction
import random

def rand_vector_sum_k_fraction(dim=6, total=Fraction(3)):
    while True:
        # Generate dim - 1 random Fractions between 0 and total
        cuts = sorted(Fraction(random.uniform(0, float(total))) for _ in range(dim - 1))
        parts = [cuts[0]] + \
                [cuts[i] - cuts[i - 1] for i in range(1, dim - 1)] + \
                [total - cuts[-1]]
        if all(x < 1 for x in parts):
            parts = sorted(parts, reverse=True)
            return parts


def perturb_vector_fraction(p, delta=Fraction(1, 10)):
    
    d = len(p)
    k = int(sum(p) + 0.1)

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
            p_aux = sorted(p_aux, reverse=True)
            return p_aux
        else:
            delta *= Fraction(95, 100)

if __name__ == "__main__":
    
    p = rand_vector_sum_k_fraction(dim=6, total=Fraction(3))
    p_ = perturb_vector_fraction(p, delta=Fraction(1, 10), k=3)
    
    mu = lexicographical_distribution(p)
    mu_ = lexicographical_distribution(p_)