from fractions import Fraction
import itertools
from math import prod
import random

if __name__ == "__main__":
    
    encontrado = False
    intentos = 0

    while not encontrado:
        
        N = random.randint(5, 13)
        K = random.randint(2, N-1)
        all_sets = all_sets = list(itertools.combinations(range(N), K))
        print(all_sets)

        intentos += 1

        Q = [Fraction(random.randint(500000, 1000000)) for _ in range(N)]
        z = sum(prod(Q[j] for j in s) for s in all_sets)
        p = [sum(prod(Q[j] for j in s) for s in all_sets if i in s) / z for i in range(N)]

        Q_ = [Q[i] + Fraction(random.randint(1, 5)) if i < K else Q[i] for i in range(N)]
        z_ = sum(prod(Q_[j] for j in s) for s in all_sets)
        p_ = [sum(prod(Q_[j] for j in s) for s in all_sets if i in s) / z_ for i in range(N)]

        while not ((p_[i] >= p[i] for i in range(K)) and (p_[i] <= p[i] for i in range(K, N))):
            Q_ = [Q[i] + Fraction(random.randint(1, 5)) if i < K else Q[i] for i in range(N)]
            z_ = sum(prod(Q_[j] for j in s) for s in all_sets)
            p_ = [sum(prod(Q_[j] for j in s) for s in all_sets if i in s) / z_ for i in range(N)]
        
        print(f"{intentos} intentos realizados.")
        
        if (prod(Q_[i] for i in range(K))/z < prod(Q[i] for i in range(K))/z_):
            
            print("Contraejemplo encontrado!")
            
            encontrado = True
            
            print("N  = ", N)
            print("K  = ", K)
            print("Q  = ", Q)
            print("Q' = ", Q_)
            print("Z  = ", z)
            print("Z' = ", z_)
            
            for i in range(N):
                print(f"p({i+1}) = {(p[i])}, p'({i+1}) = {(p_[i])}")
                print(f"p({i+1}) = {float(p[i])}, p'({i+1}) = {float(p_[i])}")
                print(p[i] < p_[i])
            
            with open("output.txt", "w") as f:
                f.write(f"N  = {N}\n")
                f.write(f"K  = {K}\n")
                f.write(f"Q  = {Q}\n")
                f.write(f"Q' = {Q_}\n")
                f.write(f"Z  = {z}\n")
                f.write(f"Z' = {z_}\n\n")
                
                for i in range(N):
                    f.write(f"p({i+1}) = {p[i]}, p'({i+1}) = {p_[i]}\n")
                    f.write(f"p({i+1}) = {float(p[i])}, p'({i+1}) = {float(p_[i])}\n")
                    f.write(f"{p[i] < p_[i]}\n")