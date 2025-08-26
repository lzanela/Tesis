from fractions import Fraction
import math
import random
import itertools

TRIALS = 1000000

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

    intentos = 0
    n = random.randint(5, 10)
    k = random.randint(min(n-1,3), n - 1)
    j = random.randint(min(k-1,2), k - 1)
    A = set(range(j))
    all_sets = list(itertools.combinations(range(n), k))
    not_contraejemplo = True
    
    while not_contraejemplo and intentos < TRIALS:
        # Generar un vector aleatorio de Fracciones que sume k
        p = rand_vector_sum_k_fraction(dim=n, total=Fraction(k))
        # Perturbar el vector para que las primeras j componentes aumenten y las restantes disminuyan
        p_ = perturb_vector_fraction(p, delta=Fraction(1, 10), k=j)

        order_condition = all(p[i] <= p_[i] for i in range(j)) and all(p_[i] <= p[i] for i in range(j, n))
        sum_condition = sum(p) == sum(p_) == k
        print(f"Intento {intentos}")
        intentos += 1

        if order_condition and sum_condition:
            suma_denom = 0
            suma_denom_ = 0
            suma_num = 0
            suma_num_ = 0

            for subset in all_sets:
                A_ = set(subset)
                prob_p = f(p, A_)
                prob_p_ = f(p_, A_)
                suma_denom += prob_p
                suma_denom_ += prob_p_

                if A.issubset(A_):
                    suma_num += prob_p
                    suma_num_ += prob_p_

            prob_p = suma_num / suma_denom
            prob_p_ = suma_num_ / suma_denom_

            if prob_p > prob_p_:
                
                print("Se encontró un contraejemplo:")
                print("N = ", n)
                print("k = ", k)
                print("j = ", j)
                print("A = ", A)

                # Imprimir los resultados
                print("p  = ", p)
                print("p' = ", p_)
                print("p as floats  = ", list(map(float, p)))
                print("p' as floats = ", list(map(float, p)))
                print("P_p([j] in S) = ", prob_p)
                print("P_p'([j] in S) = ", prob_p_)

                not_contraejemplo = False