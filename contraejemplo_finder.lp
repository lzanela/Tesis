\ Model contraejemplo_finder
\ LP format - for model browsing. Use MPS format to capture full model detail.
\ Signature: 0x7b92baeeccadbd50
Maximize
  eps + alpha[0] + alpha[1] + alpha[2] + alpha[3] + alpha[4]
Subject To
 p_prime_in_range_0: - eps + alpha[0] <= -1e-06
 p_prime_in_range_1: - eps + alpha[1] <= -1e-06
 p_prime_in_range_2: - eps + alpha[2] <= -1e-06
 p_prime_in_range_3: - eps + alpha[3] <= -1e-06
 p_prime_in_range_4: - eps + alpha[4] <= -1e-06
 feasibility_p_row_0: 0.25 eps >= -0.25
 feasibility_p_prime_row_0: 0.25 eps + alpha[1] >= -0.25
 feasibility_p_row_1: - 0.75 eps >= -0.25
 feasibility_p_prime_row_1: - 0.75 eps + 0.5 alpha[0] - 0.5 alpha[1]
   - 0.5 alpha[2] + 0.5 alpha[3] + 0.5 alpha[4] >= -0.25
 feasibility_p_row_2: - 0.75 eps >= -0.25
 feasibility_p_prime_row_2: - 0.75 eps + 0.5 alpha[0] - 0.5 alpha[1]
   + 0.5 alpha[2] - 0.5 alpha[3] + 0.5 alpha[4] >= -0.25
 feasibility_p_row_3: 0.25 eps >= -0.25
 feasibility_p_prime_row_3: 0.25 eps - alpha[4] >= -0.25
 feasibility_p_row_4: eps >= 0
 feasibility_p_prime_row_4: eps - 0.5 alpha[0] + 0.5 alpha[1]
   - 0.5 alpha[2] - 0.5 alpha[3] - 0.5 alpha[4] >= 0
 selection_monotonicity_violation: alpha[1] <= -1e-06
Bounds
 1e-06 <= eps <= 0.2
 alpha[0] <= 1
 alpha[1] <= 1
 alpha[2] <= 1
 alpha[3] <= 1
 alpha[4] <= 1
End
