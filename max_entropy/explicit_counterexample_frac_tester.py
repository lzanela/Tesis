from fractions import Fraction
import itertools
import math

# Auxiliary functions
def prod(s, q):
    return math.prod([q[i] for i in s])

def p(q, i):
    return sum([prod(s, q) for s in all_sets if i in s])


# Constants
N = 10
K = 3
S = tuple(range(K))
all_sets = list(itertools.combinations(range(N), K))


# Counterexample
q = [22.236543497581884, 22.236464321505565, 19.009094793387913, 20.420862051422475, 20.436600942729598, 20.473442188977877, 20.516409107104156, 20.549416517809327, 20.643922807410245, 16.355484636937295]
q_ = [22.236543497596703, 22.236464321520383, 19.009094793327144, 20.42086205142962, 20.43660094273796, 20.473442188985494, 20.51640910711013, 20.54941651781356, 20.64392280741342, 16.355484636939636]
max_denominator = 10**6

q = [Fraction(x).limit_denominator(max_denominator) for x in q]
q_ = [Fraction(x).limit_denominator(max_denominator) for x in q_]

# Jamie's counterexample
#q = list(map(Fraction, [99620001435175085845613951348591, 33206667145059699577734936400435, 33206667145059699577734936400435, 23244667001544291253373835102276586, 23244667001544291253373835102276586, 1660333357252963458777541885429371]))
#q = q + [Fraction(0)] * (N - len(q))
#q_prod_normalization = (1/sum([prod(s, q) for s in all_sets]))^{1/K}
#q = [x * q_prod_normalization for x in q]
#q_ = list(map(Fraction, [99620001435175193801835755646020, 33206667145059681577227243883092, 33206667145059681577227243883092, 23244667001544299141767505142336500, 23244667001544299141767505142336500, 1660333357252962147206216649823732]))
#q_ = q_ + [Fraction(0)] * (N - len(q_))
#q__prod_normalization = (1/sum([prod(s, q_) for s in all_sets]))^{1/K}
#q_ = [x * q__prod_normalization for x in q_]

if __name__ == "__main__":
    
    print("-" * 40)
    for i, f in enumerate(q):
        print(f"q[{i}]  = {f}  ≈ {float(f)}")
    
    q_normalization = sum([prod(s, q) for s in all_sets])
    probas = [p(q, i)/q_normalization for i in range(N)]
    selection_prob = prod(S, q)/q_normalization

    q__normalization = sum([prod(s, q_) for s in all_sets])
    probas_ = [p(q_, i)/q__normalization for i in range(N)]
    selection_prob_ = prod(S, q_)/q__normalization

    print("p and p' as floats:")
    print(list(map(float, probas)))
    print(list(map(float, probas_)))
    
    print("p and p' as fracs:")
    print("\n".join(map(str, probas)))
    print("\n".join(map(str, probas_)))
    
    print("Sum of p(q):")
    print(sum(probas))
    
    print("Selection prob (q):")
    print(selection_prob)
    print(float(selection_prob))
    
    print("Selection prob (q'):")
    print(selection_prob_)
    print(float(selection_prob_))
    
    print("")
    print("-" * 40)
    for i, f in enumerate(q_):
        print(f"q'[{i}]  = {f}  ≈ {float(f)}")

    # Define subconditions
    cond1 = selection_prob_ < selection_prob
    cond2 = probas_[0] >= probas[0]
    cond3 = probas_[1] >= probas[1]
    cond4 = probas_[2] >= probas[2]
    cond5 = probas_[3] <= probas[3]
    cond6 = probas_[4] <= probas[4]
    cond7 = probas_[5] <= probas[5]
    cond8 = probas_[6] <= probas[6]
    cond9 = probas_[7] <= probas[7]
    cond10 = probas_[8] <= probas[8]
    cond11 = probas_[9] <= probas[9]

    # Put them in a list with labels
    checks = [
        ("selection_prob_ < selection_prob", cond1),
        ("probas_[0] >= probas[0]", cond2),
        ("probas_[1] >= probas[1]", cond3),
        ("probas_[2] >= probas[2]", cond4),
        ("probas_[3] <= probas[3]", cond5),
        ("probas_[4] <= probas[4]", cond6),
        ("probas_[5] <= probas[5]", cond7),
        ("probas_[6] <= probas[6]", cond8),
        ("probas_[7] <= probas[7]", cond9),
        ("probas_[8] <= probas[8]", cond10),
        ("probas_[9] <= probas[9]", cond11),
    ]

    if all(cond for _, cond in checks):
        print(f"Found!")
        #print('q normalization', q__normalization)
        #print("p as floats :", list(map(float, probas)))
        #print("p' as floats:", list(map(float, probas_)))
        #print("p as fracs:")
        #print("\n".join(map(str, probas)))
        #print("p' as fracs:")
        #print("\n".join(map(str, probas_)))
        print('Sum of probabilities', 'p:', sum(probas), ", p'",sum(probas_))
        #print("Selection prob p :", selection_prob)
        #print("Selection prob p':", selection_prob_)
        print("Selection prob diff:", selection_prob_ - selection_prob)
        print("\n\n\n")
        p_diffs = [probas_[i] - probas[i] for i in range(N)]
        set_diffs = {s: prod(s, q_)/q__normalization - prod(s, q)/q_normalization for s in all_sets}
        print("p_diffs:")
        for i in range(N):
            print(f"{i}: {p_diffs[i]}")
        print("\nset_diffs:")
        for k, v in set_diffs.items():
            print(f"{k}: {v}")
    else:
        print("Not a counterexample.")
        for label, cond in checks:
            if not cond:
                print(f"   - {label} is FALSE")
