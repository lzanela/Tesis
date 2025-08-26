from fractions import Fraction
import random
import math
import itertools

# Variables globales
n = 6
k = 3
j = 2
all_sets = list(itertools.combinations(range(n), k))

def f(p: list[Fraction], A: set) -> Fraction:
    """Función f que aparece en la fórmula de probabilidad de Sampford."""
    return (k - sum((p[i]) for i in A)) * math.prod(p[i] for i in A) * math.prod((1 - p[i]) for i in range(n) if i not in A)

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
    
    not_contraejemplo = True
    A = set(range(k))
    
    while not_contraejemplo:
        # Generar un vector aleatorio de Fracciones que sume k
        p = rand_vector_sum_k_fraction(dim=n, total=Fraction(k))
        # Perturbar el vector para que las primeras j componentes aumenten y las restantes disminuyan
        p_ = perturb_vector_fraction(p, delta=Fraction(1, 10), k=j)

        # Imprimir los resultados
        print("p =  ", p)
        print("p' = ", p_)
        f_p = f(p, A)
        f_p_ = f(p_, A)
        
        print("f(p, [k]):", f_p)
        print("f(p', [k]):", f_p_)

        # Calcular f(p,A') y f(p', A') para todos los conjuntos A'
        for A_ in all_sets:
            if A_ == A:
                continue
            if f_p * f(p, A_) > f_p_ * f(p_, A_):
                not_contraejemplo = False
                print(f"Counterexample found for A = {A_}: f(p, A) * f(p', [k]) > f(p', A) * f(p, [k])")
                print(f"f(p, {A_}): {f(p, A_)}")
                print(f"f(p', {A_}): {f(p_, A_)}")
                break