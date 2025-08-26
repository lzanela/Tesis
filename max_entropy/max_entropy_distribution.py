import itertools
from cvxpy import Variable, Maximize, Problem, entr, log
from p_tqdm import p_map
from fractions import Fraction

if (1/2)*2 != 1:
    raise Exception("Need to run with Python 3.")


def max_entropy(residue):
    n = len(residue)
    assignments = list(itertools.combinations(range(n), int(sum(residue.values()) + 0.0001)))
    p = {ass : Variable() for ass in assignments}
    constraints = [p[ass] >= 0 for ass in assignments] \
        + [sum(p[ass] for ass in assignments) == 1] \
        + [sum(p[ass] for ass in assignments if i in ass) == residue[i] for i in range(n)]
    obj = Maximize(sum(entr(p[ass]) for ass in assignments))
    prob = Problem(obj, constraints)
    prob.solve()
    return {ass : float(p[ass].value) for ass in assignments}


def max_entropy_fraction(residue, max_iter=200, tol=1e-12):
    """
    Args:
        residue (dict): A dictionary mapping index `i` to its target marginal.
        max_iter (int): The maximum number of iterations.
        tol (float): The tolerance for checking convergence.

    Returns:
        dict: A dictionary mapping each assignment (tuple) to its probability (Fraction).
    
    Raises:
        ValueError: If sum of residues is not an integer.
        RuntimeError: If the algorithm fails to converge.
    """
    print(residue)
    n = len(residue)
    residue_frac = {k: Fraction(v).limit_denominator() for k, v in residue.items()}

    # --- Robust check for integer sum (from previous discussion) ---
    k_frac = sum(residue_frac.values())
    k_rounded = round(k_frac)
    if abs(k_frac - k_rounded) > Fraction(1, 10**9):
        raise ValueError(f"Sum of residue probabilities ({float(k_frac)}) is not an integer.")
    k = k_rounded
    
    if k < 0 or k > n:
        return {}
    if k == 0:
        return {(): Fraction(1)}

    # --- Iterative Scaling with Dynamic Programming ---
    mus = [Fraction(1) for _ in range(n)]

    for iteration in range(max_iter):
        # --- DP to calculate Elementary Symmetric Polynomials ---
        # dp[i][j] will be the sum of products of choosing j items from the first i mus.
        # This is e_j({mu_0, ..., mu_{i-1}}).
        dp = [[Fraction(0)] * (k + 1) for _ in range(n + 1)]
        for i in range(n + 1):
            dp[i][0] = Fraction(1) # Base case: an empty product is 1

        for i in range(1, n + 1):
            mu_prev = mus[i - 1]
            for j in range(1, k + 1):
                # We can either:
                # 1. Not include mu_{i-1}: sum is dp[i-1][j]
                # 2. Include mu_{i-1}: sum is mu_{i-1} * (sum of j-1 items from first i-1)
                dp[i][j] = dp[i-1][j] + mu_prev * dp[i-1][j-1]
        
        # The partition function Z is the sum of products of k items from all n mus.
        Z = dp[n][k]

        if Z == 0:
            break

        # --- Calculate Model Marginals using the DP table ---
        model_marginals = [Fraction(0)] * n
        max_error = 0
        
        # To calculate the marginal for mu_i, we need e_{k-1}(S \ {mu_i}).
        # We can compute this by reversing the DP recurrence:
        # e_k(S) = e_k(S\{mu_i}) + mu_i * e_{k-1}(S\{mu_i})
        # Therefore, e_k(S\{mu_i}) = e_k(S) - mu_i * e_{k-1}(S\{mu_i}), which isn't helpful.
        # Instead, we calculate it directly. The numerator for marginal i is
        # mu_i * e_{k-1}(S \ {mu_i}).
        
        # A simple way is to re-calculate DP for each exclusion, but that is O(n^2*k).
        # A more advanced way uses forward/backward tables.
        # For clarity and a massive improvement already, we do a quick calculation.
        
        # dp_without_i[j] will be e_j(S \ {mu_i})
        for i in range(n):
            if mus[i] == 0:
                # If mu_i is 0, its marginal contribution is 0, unless target is also 0.
                model_marginals[i] = Fraction(0)
                continue

            # Calculate e_{k-1}(S \ {mu_i})
            # This is (e_k(S) - e_k(S\{mu_i})) / mu_i but requires e_k(S\{mu_i}).
            # Let's find e_{k-1} of the set without mu_i.
            # Using the identity: e_j(S) = e_j(S\{x}) + x * e_{j-1}(S\{x})
            # => e_{j-1}(S\{x}) = (e_j(S) - e_j(S\{x})) / x
            # This allows calculating one DP row based on another.
            
            # dp_row_full is e_j(S) for j=0..k
            dp_row_full = dp[n] 
            # dp_row_without_i will be e_j(S\{mu_i})
            dp_row_without_i = [Fraction(0)] * (k + 1)
            dp_row_without_i[0] = Fraction(1)
            
            for j in range(1, k + 1):
                # Inverting the recurrence: e_j(S\{x}) = e_j(S) - x*e_{j-1}(S\{x})
                dp_row_without_i[j] = dp_row_full[j] - mus[i] * dp_row_without_i[j-1]
            
            numerator = mus[i] * dp_row_without_i[k-1]
            model_marginals[i] = numerator / Z
        
        # --- Check for convergence and update mus ---
        max_error = max(abs(float(model_marginals[i] - residue_frac[i])) for i in range(n))
        if max_error < tol:
            # print(f"Converged in {iteration + 1} iterations.")
            break

        for i in range(n):
            target_marginal = residue_frac.get(i, Fraction(0))
            if model_marginals[i] > 0 and target_marginal > 0:
                mus[i] *= (target_marginal / model_marginals[i])
            elif target_marginal == 0:
                mus[i] = Fraction(0)
            # If model_marginal is 0 but target is >0, mu should increase,
            # but it is likely 0 and needs to be handled by a small nudge if stuck.
            # Current update rule handles this implicitly in later iterations.
    else:
        raise RuntimeError(f"Failed to converge within {max_iter} iterations. Max error: {max_error:.2e}")

    # --- Final Calculation: Compute explicit probabilities ---
    # This step is only done ONCE. Its cost is unavoidable if the output format
    # requires one probability per assignment.
    final_p = {}
    
    # Recalculate final Z with final mus
    dp = [[Fraction(0)] * (k + 1) for _ in range(n + 1)]
    for i in range(n + 1): dp[i][0] = Fraction(1)
    for i in range(1, n + 1):
        for j in range(1, k + 1):
            dp[i][j] = dp[i - 1][j] + mus[i - 1] * dp[i - 1][j - 1]
    final_Z = dp[n][k]

    if final_Z == 0:
        return {}

    for ass in itertools.combinations(range(n), k):
        prod = Fraction(1)
        for idx in ass:
            prod *= mus[idx]
        final_p[ass] = prod / final_Z
        
    return final_p

if __name__=="__main__":
    
    test_num = 1

    if test_num == 1:
        q = max_entropy_fraction({0: 1/3, 1: 1/3, 2: 2/3, 3: 2/3})
        print(q)
        print(float(q[(0,1)]))
        
    
    #if test_num == 2:
    #    print(max_entropy({0: 0.0618342562928861, 1: 0.0207176796116814, 2: 0.0207176796116814, 3: 0.9933997806603289, 4: 0.9933997806603289, 5: 0.9099308231630933}))