from fractions import Fraction
import itertools
from math import prod

def perturb_vector_fraction(p, delta=Fraction(1), k=3):
    d = len(p)
    assert k < d, "k must be less than the length of the vector"
    admissible = False

    shift_total = delta * k
    shift_per_rest = shift_total / (d - k)

    increase_indices = list(range(k))
    decrease_indices = list(range(k, d))

    for i in increase_indices:
        p[i] += delta
    for i in decrease_indices:
        p[i] -= shift_per_rest

    return p

N = 3
K = 2
Q = [Fraction(52), Fraction(23), Fraction(62)]
all_sets = list(itertools.combinations(range(N), K))
z = sum(prod(Q[j] for j in s) for s in all_sets)
p = [sum(prod(Q[j] for j in s) for s in all_sets if i in s) / z for i in range(N)]

print("Initial Q =", Q)
print("All sets:", all_sets)
print("Initial Z =", z)
print("Initial p =", p)
print("Sum of probabilities:", sum(p))
print("Float sum of probabilities:", float(sum(p)))
print("Individual probabilities:")
for i in range(N):
    print(f"P({i}) = {p[i]}")

if __name__ == "__main__":
    Q_ = [Q[0] + 1] + Q[1:]
    z_ = sum(prod(Q_[j] for j in s) for s in all_sets)
    p_ = [sum(prod(Q_[j] for j in s) for s in all_sets if i in s) / z_ for i in range(N)]
    sum_p = sum(p_)
    
    print("Q  = ", Q)
    print("Q' = ", Q_)
    
    for i in range(N):
        print(f"p({i+1}) = {(p[i])}, p'({i+1}) = {(p_[i])}")
        print(f"p({i+1}) = {float(p[i])}, p'({i+1}) = {float(p_[i])}")
        print(p[i] < p_[i])