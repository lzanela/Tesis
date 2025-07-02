from fractions import Fraction
import itertools
import math


all_sets = list(itertools.combinations(range(6), 3))


def prod(s, q):
    return math.prod([q[i] for i in s])


def p(q, i):
    return sum([prod(s, q) for s in all_sets if i in s])


q1 = list(map(Fraction, [99620001435175085845613951348591, 33206667145059699577734936400435, 33206667145059699577734936400435, 23244667001544291253373835102276586, 23244667001544291253373835102276586, 1660333357252963458777541885429371]))
q1_normalization = sum([prod(s, q1) for s in all_sets])
p1 = [p(q1, i)/q1_normalization for i in range(6)]

if __name__ == "__main__":
    print(q1_normalization)
    print("p1 as floats:")
    print(list(map(float, p1)))
    print("p1 as fracs:")
    print("\n".join(map(str, p1)))
    print(sum(p1))
    selection_prob = prod((0, 1, 2), q1)/q1_normalization
    print(selection_prob)
    print(float(selection_prob))
    print("")
    size = 1000000000000000000000000
    for q2 in [list(map(Fraction, [99620001435175193801835755646020, 33206667145059681577227243883092, 33206667145059681577227243883092, 23244667001544299141767505142336500, 23244667001544299141767505142336500, 1660333357252962147206216649823732]))]:
        q2_normalization = sum([prod(s, q2) for s in all_sets])
        p2 = [p(q2, i)/q2_normalization for i in range(6)]
        new_selection_prob = prod((0, 1, 2), q2)/q2_normalization
        if new_selection_prob < selection_prob \
            and p2[0] >= p1[0] and p2[1] >= p1[1] and p2[2] >= p1[2] \
            and p2[3] <= p1[3] and p2[4] <= p1[4] and p2[5] <= p1[5]:
            print(f"Found!")
            print(q2_normalization)
            print("p2 as floats:")
            print(list(map(float, p2)))
            print("p2 as fracs:")
            print("\n".join(map(str, p2)))
            print(sum(p2))
            print(new_selection_prob)
            print(float(new_selection_prob))
            print("\n\n\n")
            p_diffs = [p2[i] - p1[i] for i in range(6)]
            set_diffs = {s: prod(s, q2)/q2_normalization - prod(s, q1)/q1_normalization for s in all_sets}
            print("p_diffs:")
            for i in range(6):
                print(f"{i}: {p_diffs[i]}")
            print("\nset_diffs:")
            for k, v in set_diffs.items():
                print(f"{k}: {v}")

                
