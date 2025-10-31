\ Model contraejemplo_finder
\ LP format - for model browsing. Use MPS format to capture full model detail.
\ Signature: 0x54dd26993f63a4f0
Maximize
  alpha[0] + alpha[1]
Subject To
 p_prime_0_leq_1: alpha[0] + p[0] <= 1
 p_prime_0_geq_0: alpha[0] + p[0] >= 0
 p_prime_1_leq_1: alpha[1] + p[1] <= 1
 p_prime_1_geq_0: alpha[1] + p[1] >= 0
 p_prime_2_leq_1: - alpha[2] + p[2] <= 1
 p_prime_2_geq_0: - alpha[2] + p[2] >= 0
 p_prime_3_leq_1: - alpha[3] + p[3] <= 1
 p_prime_3_geq_0: - alpha[3] + p[3] >= 0
 p_suma_2: p[0] + p[1] + p[2] + p[3] = 2
 alpha_suma_0: alpha[0] + alpha[1] - alpha[2] - alpha[3] = 0
 feasibility_p_geq_0_row_0: p[1] >= 0
 feasibility_p_leq_1_row_0: p[1] <= 1
 feasibility_p_prime_geq_0_row_0: 0.5 alpha[0] + 0.5 alpha[1]
   + 0.5 alpha[2] + 0.5 alpha[3] + 0.5 p[0] + 0.5 p[1] - 0.5 p[2]
   - 0.5 p[3] >= 0
 feasibility_p_prime_leq_1_row_0: 0.5 alpha[0] + 0.5 alpha[1]
   + 0.5 alpha[2] + 0.5 alpha[3] + 0.5 p[0] + 0.5 p[1] - 0.5 p[2]
   - 0.5 p[3] <= 1
 feasibility_p_geq_0_row_1: 0.5 p[0] - 0.5 p[1] + 0.5 p[2] - 0.5 p[3] >= 0
 feasibility_p_leq_1_row_1: 0.5 p[0] - 0.5 p[1] + 0.5 p[2] - 0.5 p[3] <= 1
 feasibility_p_prime_geq_0_row_1: 0.5 alpha[0] - 0.5 alpha[1]
   - 0.5 alpha[2] + 0.5 alpha[3] + 0.5 p[0] - 0.5 p[1] + 0.5 p[2]
   - 0.5 p[3] >= 0
 feasibility_p_prime_leq_1_row_1: 0.5 alpha[0] - 0.5 alpha[1]
   - 0.5 alpha[2] + 0.5 alpha[3] + 0.5 p[0] - 0.5 p[1] + 0.5 p[2]
   - 0.5 p[3] <= 1
 feasibility_p_geq_0_row_2: 0.5 p[0] - 0.5 p[1] - 0.5 p[2] + 0.5 p[3] >= 0
 feasibility_p_leq_1_row_2: 0.5 p[0] - 0.5 p[1] - 0.5 p[2] + 0.5 p[3] <= 1
 feasibility_p_prime_geq_0_row_2: - alpha[3] + p[3] >= 0
 feasibility_p_prime_leq_1_row_2: - alpha[3] + p[3] <= 1
 feasibility_p_geq_0_row_3: - 0.5 p[0] + 0.5 p[1] + 0.5 p[2] + 0.5 p[3]
   >= 0
 feasibility_p_leq_1_row_3: - 0.5 p[0] + 0.5 p[1] + 0.5 p[2] + 0.5 p[3]
   <= 1
 feasibility_p_prime_geq_0_row_3: - 0.5 alpha[0] + 0.5 alpha[1]
   - 0.5 alpha[2] - 0.5 alpha[3] - 0.5 p[0] + 0.5 p[1] + 0.5 p[2]
   + 0.5 p[3] >= 0
 feasibility_p_prime_leq_1_row_3: - 0.5 alpha[0] + 0.5 alpha[1]
   - 0.5 alpha[2] - 0.5 alpha[3] - 0.5 p[0] + 0.5 p[1] + 0.5 p[2]
   + 0.5 p[3] <= 1
 infeasibility_p_base_1_at_least_one_violated: z_p_base_1[0]
   + z_p_base_1[1] + z_p_base_1[2] + z_p_base_1[3] = 1
 infeasibility_p_prime_base_1_at_least_one_violated: z_p_prime_base_1[0]
   + z_p_prime_base_1[1] + z_p_prime_base_1[2] + z_p_prime_base_1[3] = 1
 infeasibility_p_base_3_at_least_one_violated: z_p_base_3[0]
   + z_p_base_3[1] + z_p_base_3[2] + z_p_base_3[3] = 1
 infeasibility_p_prime_base_3_at_least_one_violated: z_p_prime_base_3[0]
   + z_p_prime_base_3[1] + z_p_prime_base_3[2] + z_p_prime_base_3[3] = 1
 infeasibility_p_base_4_at_least_one_violated: z_p_base_4[0]
   + z_p_base_4[1] + z_p_base_4[2] + z_p_base_4[3] = 1
 infeasibility_p_prime_base_4_at_least_one_violated: z_p_prime_base_4[0]
   + z_p_prime_base_4[1] + z_p_prime_base_4[2] + z_p_prime_base_4[3] = 1
 infeasibility_p_base_5_at_least_one_violated: z_p_base_5[0]
   + z_p_base_5[1] + z_p_base_5[2] + z_p_base_5[3] = 1
 infeasibility_p_prime_base_5_at_least_one_violated: z_p_prime_base_5[0]
   + z_p_prime_base_5[1] + z_p_prime_base_5[2] + z_p_prime_base_5[3] = 1
 infeasibility_p_base_6_at_least_one_violated: z_p_base_6[0]
   + z_p_base_6[1] + z_p_base_6[2] + z_p_base_6[3] = 1
 infeasibility_p_prime_base_6_at_least_one_violated: z_p_prime_base_6[0]
   + z_p_prime_base_6[1] + z_p_prime_base_6[2] + z_p_prime_base_6[3] = 1
 infeasibility_p_base_7_at_least_one_violated: z_p_base_7[0]
   + z_p_base_7[1] + z_p_base_7[2] + z_p_base_7[3] = 1
 infeasibility_p_prime_base_7_at_least_one_violated: z_p_prime_base_7[0]
   + z_p_prime_base_7[1] + z_p_prime_base_7[2] + z_p_prime_base_7[3] = 1
 infeasibility_p_base_8_at_least_one_violated: z_p_base_8[0]
   + z_p_base_8[1] + z_p_base_8[2] + z_p_base_8[3] = 1
 infeasibility_p_prime_base_8_at_least_one_violated: z_p_prime_base_8[0]
   + z_p_prime_base_8[1] + z_p_prime_base_8[2] + z_p_prime_base_8[3] = 1
 infeasibility_p_base_9_at_least_one_violated: z_p_base_9[0]
   + z_p_base_9[1] + z_p_base_9[2] + z_p_base_9[3] = 1
 infeasibility_p_prime_base_9_at_least_one_violated: z_p_prime_base_9[0]
   + z_p_prime_base_9[1] + z_p_prime_base_9[2] + z_p_prime_base_9[3] = 1
 infeasibility_p_base_10_at_least_one_violated: z_p_base_10[0]
   + z_p_base_10[1] + z_p_base_10[2] + z_p_base_10[3] = 1
 infeasibility_p_prime_base_10_at_least_one_violated: z_p_prime_base_10[0]
   + z_p_prime_base_10[1] + z_p_prime_base_10[2] + z_p_prime_base_10[3]
   = 1
 infeasibility_p_base_11_at_least_one_violated: z_p_base_11[0]
   + z_p_base_11[1] + z_p_base_11[2] + z_p_base_11[3] = 1
 infeasibility_p_prime_base_11_at_least_one_violated: z_p_prime_base_11[0]
   + z_p_prime_base_11[1] + z_p_prime_base_11[2] + z_p_prime_base_11[3]
   = 1
 infeasibility_p_base_B_prime_at_least_one_violated: z_p_base_B_prime[0]
   + z_p_base_B_prime[1] + z_p_base_B_prime[2] + z_p_base_B_prime[3] = 1
 infeasibility_p_prime_base_B_at_least_one_violated: z_p_prime_base_B[0]
   + z_p_prime_base_B[1] + z_p_prime_base_B[2] + z_p_prime_base_B[3] = 1
 selection_monotonicity_violation: 0.5 alpha[0] + 0.5 alpha[1]
   + 0.5 alpha[2] + 0.5 alpha[3] + 0.5 p[0] - 0.5 p[1] - 0.5 p[2]
   - 0.5 p[3] <= -1e-07
 infeasibility_p_base_1_row_0_leq_0: 1e-07 z_p_base_1[0] + [
   0.5 p[0] * z_p_base_1[0] + 0.5 p[1] * z_p_base_1[0]
   - 0.5 p[2] * z_p_base_1[0] - 0.5 p[3] * z_p_base_1[0] ] <= 0
 infeasibility_p_base_1_row_1_leq_0: 1e-07 z_p_base_1[1] + [
   p[2] * z_p_base_1[1] ] <= 0
 infeasibility_p_base_1_row_2_leq_0: 1e-07 z_p_base_1[2] + [
   0.5 p[0] * z_p_base_1[2] - 0.5 p[1] * z_p_base_1[2]
   - 0.5 p[2] * z_p_base_1[2] + 0.5 p[3] * z_p_base_1[2] ] <= 0
 infeasibility_p_base_1_row_3_leq_0: 1e-07 z_p_base_1[3] + [
   - 0.5 p[0] * z_p_base_1[3] + 0.5 p[1] * z_p_base_1[3]
   + 0.5 p[2] * z_p_base_1[3] + 0.5 p[3] * z_p_base_1[3] ] <= 0
 infeasibility_p_prime_base_1_row_0_leq_0: 1e-07 z_p_prime_base_1[0] + [
   0.5 alpha[0] * z_p_prime_base_1[0] + 0.5 alpha[1] * z_p_prime_base_1[0]
   + 0.5 alpha[2] * z_p_prime_base_1[0]
   + 0.5 alpha[3] * z_p_prime_base_1[0] + 0.5 p[0] * z_p_prime_base_1[0]
   + 0.5 p[1] * z_p_prime_base_1[0] - 0.5 p[2] * z_p_prime_base_1[0]
   - 0.5 p[3] * z_p_prime_base_1[0] ] <= 0
 infeasibility_p_prime_base_1_row_1_leq_0: 1e-07 z_p_prime_base_1[1] + [
   - alpha[2] * z_p_prime_base_1[1] + p[2] * z_p_prime_base_1[1] ] <= 0
 infeasibility_p_prime_base_1_row_2_leq_0: 1e-07 z_p_prime_base_1[2] + [
   0.5 alpha[0] * z_p_prime_base_1[2] - 0.5 alpha[1] * z_p_prime_base_1[2]
   + 0.5 alpha[2] * z_p_prime_base_1[2]
   - 0.5 alpha[3] * z_p_prime_base_1[2] + 0.5 p[0] * z_p_prime_base_1[2]
   - 0.5 p[1] * z_p_prime_base_1[2] - 0.5 p[2] * z_p_prime_base_1[2]
   + 0.5 p[3] * z_p_prime_base_1[2] ] <= 0
 infeasibility_p_prime_base_1_row_3_leq_0: 1e-07 z_p_prime_base_1[3] + [
   - 0.5 alpha[0] * z_p_prime_base_1[3]
   + 0.5 alpha[1] * z_p_prime_base_1[3]
   - 0.5 alpha[2] * z_p_prime_base_1[3]
   - 0.5 alpha[3] * z_p_prime_base_1[3] - 0.5 p[0] * z_p_prime_base_1[3]
   + 0.5 p[1] * z_p_prime_base_1[3] + 0.5 p[2] * z_p_prime_base_1[3]
   + 0.5 p[3] * z_p_prime_base_1[3] ] <= 0
 infeasibility_p_base_3_row_0_leq_0: 1e-07 z_p_base_3[0] + [
   0.5 p[0] * z_p_base_3[0] + 0.5 p[1] * z_p_base_3[0]
   - 0.5 p[2] * z_p_base_3[0] - 0.5 p[3] * z_p_base_3[0] ] <= 0
 infeasibility_p_base_3_row_1_leq_0: 1e-07 z_p_base_3[1] + [
   0.5 p[0] * z_p_base_3[1] - 0.5 p[1] * z_p_base_3[1]
   + 0.5 p[2] * z_p_base_3[1] + 0.5 p[3] * z_p_base_3[1] ] <= 0
 infeasibility_p_base_3_row_2_leq_0: 1e-07 z_p_base_3[2] + [
   - 0.5 p[0] * z_p_base_3[2] + 0.5 p[1] * z_p_base_3[2]
   + 0.5 p[2] * z_p_base_3[2] - 0.5 p[3] * z_p_base_3[2] ] <= 0
 infeasibility_p_base_3_row_3_leq_0: 1e-07 z_p_base_3[3] + [
   p[3] * z_p_base_3[3] ] <= 0
 infeasibility_p_prime_base_3_row_0_leq_0: 1e-07 z_p_prime_base_3[0] + [
   0.5 alpha[0] * z_p_prime_base_3[0] + 0.5 alpha[1] * z_p_prime_base_3[0]
   + 0.5 alpha[2] * z_p_prime_base_3[0]
   + 0.5 alpha[3] * z_p_prime_base_3[0] + 0.5 p[0] * z_p_prime_base_3[0]
   + 0.5 p[1] * z_p_prime_base_3[0] - 0.5 p[2] * z_p_prime_base_3[0]
   - 0.5 p[3] * z_p_prime_base_3[0] ] <= 0
 infeasibility_p_prime_base_3_row_1_leq_0: 1e-07 z_p_prime_base_3[1] + [
   0.5 alpha[0] * z_p_prime_base_3[1] - 0.5 alpha[1] * z_p_prime_base_3[1]
   - 0.5 alpha[2] * z_p_prime_base_3[1]
   - 0.5 alpha[3] * z_p_prime_base_3[1] + 0.5 p[0] * z_p_prime_base_3[1]
   - 0.5 p[1] * z_p_prime_base_3[1] + 0.5 p[2] * z_p_prime_base_3[1]
   + 0.5 p[3] * z_p_prime_base_3[1] ] <= 0
 infeasibility_p_prime_base_3_row_2_leq_0: 1e-07 z_p_prime_base_3[2] + [
   - 0.5 alpha[0] * z_p_prime_base_3[2]
   + 0.5 alpha[1] * z_p_prime_base_3[2]
   - 0.5 alpha[2] * z_p_prime_base_3[2]
   + 0.5 alpha[3] * z_p_prime_base_3[2] - 0.5 p[0] * z_p_prime_base_3[2]
   + 0.5 p[1] * z_p_prime_base_3[2] + 0.5 p[2] * z_p_prime_base_3[2]
   - 0.5 p[3] * z_p_prime_base_3[2] ] <= 0
 infeasibility_p_prime_base_3_row_3_leq_0: 1e-07 z_p_prime_base_3[3] + [
   - alpha[3] * z_p_prime_base_3[3] + p[3] * z_p_prime_base_3[3] ] <= 0
 infeasibility_p_base_4_row_0_leq_0: 1e-07 z_p_base_4[0] + [
   0.5 p[0] * z_p_base_4[0] + 0.5 p[1] * z_p_base_4[0]
   - 0.5 p[2] * z_p_base_4[0] + 0.5 p[3] * z_p_base_4[0] ] <= 0
 infeasibility_p_base_4_row_1_leq_0: 1e-07 z_p_base_4[1] + [
   0.5 p[0] * z_p_base_4[1] - 0.5 p[1] * z_p_base_4[1]
   + 0.5 p[2] * z_p_base_4[1] - 0.5 p[3] * z_p_base_4[1] ] <= 0
 infeasibility_p_base_4_row_2_leq_0: 1e-07 z_p_base_4[2] + [
   - 0.5 p[0] * z_p_base_4[2] + 0.5 p[1] * z_p_base_4[2]
   + 0.5 p[2] * z_p_base_4[2] - 0.5 p[3] * z_p_base_4[2] ] <= 0
 infeasibility_p_base_4_row_3_leq_0: 1e-07 z_p_base_4[3] + [
   p[3] * z_p_base_4[3] ] <= 0
 infeasibility_p_prime_base_4_row_0_leq_0: 1e-07 z_p_prime_base_4[0] + [
   0.5 alpha[0] * z_p_prime_base_4[0] + 0.5 alpha[1] * z_p_prime_base_4[0]
   + 0.5 alpha[2] * z_p_prime_base_4[0]
   - 0.5 alpha[3] * z_p_prime_base_4[0] + 0.5 p[0] * z_p_prime_base_4[0]
   + 0.5 p[1] * z_p_prime_base_4[0] - 0.5 p[2] * z_p_prime_base_4[0]
   + 0.5 p[3] * z_p_prime_base_4[0] ] <= 0
 infeasibility_p_prime_base_4_row_1_leq_0: 1e-07 z_p_prime_base_4[1] + [
   0.5 alpha[0] * z_p_prime_base_4[1] - 0.5 alpha[1] * z_p_prime_base_4[1]
   - 0.5 alpha[2] * z_p_prime_base_4[1]
   + 0.5 alpha[3] * z_p_prime_base_4[1] + 0.5 p[0] * z_p_prime_base_4[1]
   - 0.5 p[1] * z_p_prime_base_4[1] + 0.5 p[2] * z_p_prime_base_4[1]
   - 0.5 p[3] * z_p_prime_base_4[1] ] <= 0
 infeasibility_p_prime_base_4_row_2_leq_0: 1e-07 z_p_prime_base_4[2] + [
   - 0.5 alpha[0] * z_p_prime_base_4[2]
   + 0.5 alpha[1] * z_p_prime_base_4[2]
   - 0.5 alpha[2] * z_p_prime_base_4[2]
   + 0.5 alpha[3] * z_p_prime_base_4[2] - 0.5 p[0] * z_p_prime_base_4[2]
   + 0.5 p[1] * z_p_prime_base_4[2] + 0.5 p[2] * z_p_prime_base_4[2]
   - 0.5 p[3] * z_p_prime_base_4[2] ] <= 0
 infeasibility_p_prime_base_4_row_3_leq_0: 1e-07 z_p_prime_base_4[3] + [
   - alpha[3] * z_p_prime_base_4[3] + p[3] * z_p_prime_base_4[3] ] <= 0
 infeasibility_p_base_5_row_0_leq_0: 1e-07 z_p_base_5[0] + [
   0.5 p[0] * z_p_base_5[0] + 0.5 p[1] * z_p_base_5[0]
   - 0.5 p[2] * z_p_base_5[0] - 0.5 p[3] * z_p_base_5[0] ] <= 0
 infeasibility_p_base_5_row_1_leq_0: 1e-07 z_p_base_5[1] + [
   0.5 p[0] * z_p_base_5[1] - 0.5 p[1] * z_p_base_5[1]
   + 0.5 p[2] * z_p_base_5[1] + 0.5 p[3] * z_p_base_5[1] ] <= 0
 infeasibility_p_base_5_row_2_leq_0: 1e-07 z_p_base_5[2] + [
   p[2] * z_p_base_5[2] ] <= 0
 infeasibility_p_base_5_row_3_leq_0: 1e-07 z_p_base_5[3] + [
   - 0.5 p[0] * z_p_base_5[3] + 0.5 p[1] * z_p_base_5[3]
   - 0.5 p[2] * z_p_base_5[3] + 0.5 p[3] * z_p_base_5[3] ] <= 0
 infeasibility_p_prime_base_5_row_0_leq_0: 1e-07 z_p_prime_base_5[0] + [
   0.5 alpha[0] * z_p_prime_base_5[0] + 0.5 alpha[1] * z_p_prime_base_5[0]
   + 0.5 alpha[2] * z_p_prime_base_5[0]
   + 0.5 alpha[3] * z_p_prime_base_5[0] + 0.5 p[0] * z_p_prime_base_5[0]
   + 0.5 p[1] * z_p_prime_base_5[0] - 0.5 p[2] * z_p_prime_base_5[0]
   - 0.5 p[3] * z_p_prime_base_5[0] ] <= 0
 infeasibility_p_prime_base_5_row_1_leq_0: 1e-07 z_p_prime_base_5[1] + [
   0.5 alpha[0] * z_p_prime_base_5[1] - 0.5 alpha[1] * z_p_prime_base_5[1]
   - 0.5 alpha[2] * z_p_prime_base_5[1]
   - 0.5 alpha[3] * z_p_prime_base_5[1] + 0.5 p[0] * z_p_prime_base_5[1]
   - 0.5 p[1] * z_p_prime_base_5[1] + 0.5 p[2] * z_p_prime_base_5[1]
   + 0.5 p[3] * z_p_prime_base_5[1] ] <= 0
 infeasibility_p_prime_base_5_row_2_leq_0: 1e-07 z_p_prime_base_5[2] + [
   - alpha[2] * z_p_prime_base_5[2] + p[2] * z_p_prime_base_5[2] ] <= 0
 infeasibility_p_prime_base_5_row_3_leq_0: 1e-07 z_p_prime_base_5[3] + [
   - 0.5 alpha[0] * z_p_prime_base_5[3]
   + 0.5 alpha[1] * z_p_prime_base_5[3]
   + 0.5 alpha[2] * z_p_prime_base_5[3]
   - 0.5 alpha[3] * z_p_prime_base_5[3] - 0.5 p[0] * z_p_prime_base_5[3]
   + 0.5 p[1] * z_p_prime_base_5[3] - 0.5 p[2] * z_p_prime_base_5[3]
   + 0.5 p[3] * z_p_prime_base_5[3] ] <= 0
 infeasibility_p_base_6_row_0_leq_0: 1e-07 z_p_base_6[0] + [
   0.5 p[0] * z_p_base_6[0] + 0.5 p[1] * z_p_base_6[0]
   + 0.5 p[2] * z_p_base_6[0] - 0.5 p[3] * z_p_base_6[0] ] <= 0
 infeasibility_p_base_6_row_1_leq_0: 1e-07 z_p_base_6[1] + [
   0.5 p[0] * z_p_base_6[1] - 0.5 p[1] * z_p_base_6[1]
   - 0.5 p[2] * z_p_base_6[1] + 0.5 p[3] * z_p_base_6[1] ] <= 0
 infeasibility_p_base_6_row_2_leq_0: 1e-07 z_p_base_6[2] + [
   - 0.5 p[0] * z_p_base_6[2] + 0.5 p[1] * z_p_base_6[2]
   - 0.5 p[2] * z_p_base_6[2] + 0.5 p[3] * z_p_base_6[2] ] <= 0
 infeasibility_p_base_6_row_3_leq_0: 1e-07 z_p_base_6[3] + [
   p[2] * z_p_base_6[3] ] <= 0
 infeasibility_p_prime_base_6_row_0_leq_0: 1e-07 z_p_prime_base_6[0] + [
   0.5 alpha[0] * z_p_prime_base_6[0] + 0.5 alpha[1] * z_p_prime_base_6[0]
   - 0.5 alpha[2] * z_p_prime_base_6[0]
   + 0.5 alpha[3] * z_p_prime_base_6[0] + 0.5 p[0] * z_p_prime_base_6[0]
   + 0.5 p[1] * z_p_prime_base_6[0] + 0.5 p[2] * z_p_prime_base_6[0]
   - 0.5 p[3] * z_p_prime_base_6[0] ] <= 0
 infeasibility_p_prime_base_6_row_1_leq_0: 1e-07 z_p_prime_base_6[1] + [
   0.5 alpha[0] * z_p_prime_base_6[1] - 0.5 alpha[1] * z_p_prime_base_6[1]
   + 0.5 alpha[2] * z_p_prime_base_6[1]
   - 0.5 alpha[3] * z_p_prime_base_6[1] + 0.5 p[0] * z_p_prime_base_6[1]
   - 0.5 p[1] * z_p_prime_base_6[1] - 0.5 p[2] * z_p_prime_base_6[1]
   + 0.5 p[3] * z_p_prime_base_6[1] ] <= 0
 infeasibility_p_prime_base_6_row_2_leq_0: 1e-07 z_p_prime_base_6[2] + [
   - 0.5 alpha[0] * z_p_prime_base_6[2]
   + 0.5 alpha[1] * z_p_prime_base_6[2]
   + 0.5 alpha[2] * z_p_prime_base_6[2]
   - 0.5 alpha[3] * z_p_prime_base_6[2] - 0.5 p[0] * z_p_prime_base_6[2]
   + 0.5 p[1] * z_p_prime_base_6[2] - 0.5 p[2] * z_p_prime_base_6[2]
   + 0.5 p[3] * z_p_prime_base_6[2] ] <= 0
 infeasibility_p_prime_base_6_row_3_leq_0: 1e-07 z_p_prime_base_6[3] + [
   - alpha[2] * z_p_prime_base_6[3] + p[2] * z_p_prime_base_6[3] ] <= 0
 infeasibility_p_base_7_row_0_leq_0: 1e-07 z_p_base_7[0] + [
   p[0] * z_p_base_7[0] ] <= 0
 infeasibility_p_base_7_row_1_leq_0: 1e-07 z_p_base_7[1] + [
   - 0.5 p[0] * z_p_base_7[1] + 0.5 p[1] * z_p_base_7[1]
   + 0.5 p[2] * z_p_base_7[1] - 0.5 p[3] * z_p_base_7[1] ] <= 0
 infeasibility_p_base_7_row_2_leq_0: 1e-07 z_p_base_7[2] + [
   - 0.5 p[0] * z_p_base_7[2] + 0.5 p[1] * z_p_base_7[2]
   - 0.5 p[2] * z_p_base_7[2] + 0.5 p[3] * z_p_base_7[2] ] <= 0
 infeasibility_p_base_7_row_3_leq_0: 1e-07 z_p_base_7[3] + [
   0.5 p[0] * z_p_base_7[3] - 0.5 p[1] * z_p_base_7[3]
   + 0.5 p[2] * z_p_base_7[3] + 0.5 p[3] * z_p_base_7[3] ] <= 0
 infeasibility_p_prime_base_7_row_0_leq_0: 1e-07 z_p_prime_base_7[0] + [
   alpha[0] * z_p_prime_base_7[0] + p[0] * z_p_prime_base_7[0] ] <= 0
 infeasibility_p_prime_base_7_row_1_leq_0: 1e-07 z_p_prime_base_7[1] + [
   - 0.5 alpha[0] * z_p_prime_base_7[1]
   + 0.5 alpha[1] * z_p_prime_base_7[1]
   - 0.5 alpha[2] * z_p_prime_base_7[1]
   + 0.5 alpha[3] * z_p_prime_base_7[1] - 0.5 p[0] * z_p_prime_base_7[1]
   + 0.5 p[1] * z_p_prime_base_7[1] + 0.5 p[2] * z_p_prime_base_7[1]
   - 0.5 p[3] * z_p_prime_base_7[1] ] <= 0
 infeasibility_p_prime_base_7_row_2_leq_0: 1e-07 z_p_prime_base_7[2] + [
   - 0.5 alpha[0] * z_p_prime_base_7[2]
   + 0.5 alpha[1] * z_p_prime_base_7[2]
   + 0.5 alpha[2] * z_p_prime_base_7[2]
   - 0.5 alpha[3] * z_p_prime_base_7[2] - 0.5 p[0] * z_p_prime_base_7[2]
   + 0.5 p[1] * z_p_prime_base_7[2] - 0.5 p[2] * z_p_prime_base_7[2]
   + 0.5 p[3] * z_p_prime_base_7[2] ] <= 0
 infeasibility_p_prime_base_7_row_3_leq_0: 1e-07 z_p_prime_base_7[3] + [
   0.5 alpha[0] * z_p_prime_base_7[3] - 0.5 alpha[1] * z_p_prime_base_7[3]
   - 0.5 alpha[2] * z_p_prime_base_7[3]
   - 0.5 alpha[3] * z_p_prime_base_7[3] + 0.5 p[0] * z_p_prime_base_7[3]
   - 0.5 p[1] * z_p_prime_base_7[3] + 0.5 p[2] * z_p_prime_base_7[3]
   + 0.5 p[3] * z_p_prime_base_7[3] ] <= 0
 infeasibility_p_base_8_row_0_leq_0: 1e-07 z_p_base_8[0] + [
   0.5 p[0] * z_p_base_8[0] - 0.5 p[1] * z_p_base_8[0]
   + 0.5 p[2] * z_p_base_8[0] - 0.5 p[3] * z_p_base_8[0] ] <= 0
 infeasibility_p_base_8_row_1_leq_0: 1e-07 z_p_base_8[1] + [
   0.5 p[0] * z_p_base_8[1] + 0.5 p[1] * z_p_base_8[1]
   - 0.5 p[2] * z_p_base_8[1] + 0.5 p[3] * z_p_base_8[1] ] <= 0
 infeasibility_p_base_8_row_2_leq_0: 1e-07 z_p_base_8[2] + [
   p[1] * z_p_base_8[2] ] <= 0
 infeasibility_p_base_8_row_3_leq_0: 1e-07 z_p_base_8[3] + [
   - 0.5 p[0] * z_p_base_8[3] - 0.5 p[1] * z_p_base_8[3]
   + 0.5 p[2] * z_p_base_8[3] + 0.5 p[3] * z_p_base_8[3] ] <= 0
 infeasibility_p_prime_base_8_row_0_leq_0: 1e-07 z_p_prime_base_8[0] + [
   0.5 alpha[0] * z_p_prime_base_8[0] - 0.5 alpha[1] * z_p_prime_base_8[0]
   - 0.5 alpha[2] * z_p_prime_base_8[0]
   + 0.5 alpha[3] * z_p_prime_base_8[0] + 0.5 p[0] * z_p_prime_base_8[0]
   - 0.5 p[1] * z_p_prime_base_8[0] + 0.5 p[2] * z_p_prime_base_8[0]
   - 0.5 p[3] * z_p_prime_base_8[0] ] <= 0
 infeasibility_p_prime_base_8_row_1_leq_0: 1e-07 z_p_prime_base_8[1] + [
   0.5 alpha[0] * z_p_prime_base_8[1] + 0.5 alpha[1] * z_p_prime_base_8[1]
   + 0.5 alpha[2] * z_p_prime_base_8[1]
   - 0.5 alpha[3] * z_p_prime_base_8[1] + 0.5 p[0] * z_p_prime_base_8[1]
   + 0.5 p[1] * z_p_prime_base_8[1] - 0.5 p[2] * z_p_prime_base_8[1]
   + 0.5 p[3] * z_p_prime_base_8[1] ] <= 0
 infeasibility_p_prime_base_8_row_2_leq_0: 1e-07 z_p_prime_base_8[2] + [
   alpha[1] * z_p_prime_base_8[2] + p[1] * z_p_prime_base_8[2] ] <= 0
 infeasibility_p_prime_base_8_row_3_leq_0: 1e-07 z_p_prime_base_8[3] + [
   - 0.5 alpha[0] * z_p_prime_base_8[3]
   - 0.5 alpha[1] * z_p_prime_base_8[3]
   - 0.5 alpha[2] * z_p_prime_base_8[3]
   - 0.5 alpha[3] * z_p_prime_base_8[3] - 0.5 p[0] * z_p_prime_base_8[3]
   - 0.5 p[1] * z_p_prime_base_8[3] + 0.5 p[2] * z_p_prime_base_8[3]
   + 0.5 p[3] * z_p_prime_base_8[3] ] <= 0
 infeasibility_p_base_9_row_0_leq_0: 1e-07 z_p_base_9[0] + [
   0.5 p[0] * z_p_base_9[0] + 0.5 p[1] * z_p_base_9[0]
   + 0.5 p[2] * z_p_base_9[0] - 0.5 p[3] * z_p_base_9[0] ] <= 0
 infeasibility_p_base_9_row_1_leq_0: 1e-07 z_p_base_9[1] + [
   0.5 p[0] * z_p_base_9[1] - 0.5 p[1] * z_p_base_9[1]
   - 0.5 p[2] * z_p_base_9[1] + 0.5 p[3] * z_p_base_9[1] ] <= 0
 infeasibility_p_base_9_row_2_leq_0: 1e-07 z_p_base_9[2] + [
   p[1] * z_p_base_9[2] ] <= 0
 infeasibility_p_base_9_row_3_leq_0: 1e-07 z_p_base_9[3] + [
   - 0.5 p[0] * z_p_base_9[3] - 0.5 p[1] * z_p_base_9[3]
   + 0.5 p[2] * z_p_base_9[3] + 0.5 p[3] * z_p_base_9[3] ] <= 0
 infeasibility_p_prime_base_9_row_0_leq_0: 1e-07 z_p_prime_base_9[0] + [
   0.5 alpha[0] * z_p_prime_base_9[0] + 0.5 alpha[1] * z_p_prime_base_9[0]
   - 0.5 alpha[2] * z_p_prime_base_9[0]
   + 0.5 alpha[3] * z_p_prime_base_9[0] + 0.5 p[0] * z_p_prime_base_9[0]
   + 0.5 p[1] * z_p_prime_base_9[0] + 0.5 p[2] * z_p_prime_base_9[0]
   - 0.5 p[3] * z_p_prime_base_9[0] ] <= 0
 infeasibility_p_prime_base_9_row_1_leq_0: 1e-07 z_p_prime_base_9[1] + [
   0.5 alpha[0] * z_p_prime_base_9[1] - 0.5 alpha[1] * z_p_prime_base_9[1]
   + 0.5 alpha[2] * z_p_prime_base_9[1]
   - 0.5 alpha[3] * z_p_prime_base_9[1] + 0.5 p[0] * z_p_prime_base_9[1]
   - 0.5 p[1] * z_p_prime_base_9[1] - 0.5 p[2] * z_p_prime_base_9[1]
   + 0.5 p[3] * z_p_prime_base_9[1] ] <= 0
 infeasibility_p_prime_base_9_row_2_leq_0: 1e-07 z_p_prime_base_9[2] + [
   alpha[1] * z_p_prime_base_9[2] + p[1] * z_p_prime_base_9[2] ] <= 0
 infeasibility_p_prime_base_9_row_3_leq_0: 1e-07 z_p_prime_base_9[3] + [
   - 0.5 alpha[0] * z_p_prime_base_9[3]
   - 0.5 alpha[1] * z_p_prime_base_9[3]
   - 0.5 alpha[2] * z_p_prime_base_9[3]
   - 0.5 alpha[3] * z_p_prime_base_9[3] - 0.5 p[0] * z_p_prime_base_9[3]
   - 0.5 p[1] * z_p_prime_base_9[3] + 0.5 p[2] * z_p_prime_base_9[3]
   + 0.5 p[3] * z_p_prime_base_9[3] ] <= 0
 infeasibility_p_base_10_row_0_leq_0: 1e-07 z_p_base_10[0] + [
   p[0] * z_p_base_10[0] ] <= 0
 infeasibility_p_base_10_row_1_leq_0: 1e-07 z_p_base_10[1] + [
   - 0.5 p[0] * z_p_base_10[1] + 0.5 p[1] * z_p_base_10[1]
   + 0.5 p[2] * z_p_base_10[1] - 0.5 p[3] * z_p_base_10[1] ] <= 0
 infeasibility_p_base_10_row_2_leq_0: 1e-07 z_p_base_10[2] + [
   0.5 p[0] * z_p_base_10[2] + 0.5 p[1] * z_p_base_10[2]
   - 0.5 p[2] * z_p_base_10[2] + 0.5 p[3] * z_p_base_10[2] ] <= 0
 infeasibility_p_base_10_row_3_leq_0: 1e-07 z_p_base_10[3] + [
   - 0.5 p[0] * z_p_base_10[3] - 0.5 p[1] * z_p_base_10[3]
   + 0.5 p[2] * z_p_base_10[3] + 0.5 p[3] * z_p_base_10[3] ] <= 0
 infeasibility_p_prime_base_10_row_0_leq_0: 1e-07 z_p_prime_base_10[0] + [
   alpha[0] * z_p_prime_base_10[0] + p[0] * z_p_prime_base_10[0] ] <= 0
 infeasibility_p_prime_base_10_row_1_leq_0: 1e-07 z_p_prime_base_10[1] + [
   - 0.5 alpha[0] * z_p_prime_base_10[1]
   + 0.5 alpha[1] * z_p_prime_base_10[1]
   - 0.5 alpha[2] * z_p_prime_base_10[1]
   + 0.5 alpha[3] * z_p_prime_base_10[1] - 0.5 p[0] * z_p_prime_base_10[1]
   + 0.5 p[1] * z_p_prime_base_10[1] + 0.5 p[2] * z_p_prime_base_10[1]
   - 0.5 p[3] * z_p_prime_base_10[1] ] <= 0
 infeasibility_p_prime_base_10_row_2_leq_0: 1e-07 z_p_prime_base_10[2] + [
   0.5 alpha[0] * z_p_prime_base_10[2]
   + 0.5 alpha[1] * z_p_prime_base_10[2]
   + 0.5 alpha[2] * z_p_prime_base_10[2]
   - 0.5 alpha[3] * z_p_prime_base_10[2] + 0.5 p[0] * z_p_prime_base_10[2]
   + 0.5 p[1] * z_p_prime_base_10[2] - 0.5 p[2] * z_p_prime_base_10[2]
   + 0.5 p[3] * z_p_prime_base_10[2] ] <= 0
 infeasibility_p_prime_base_10_row_3_leq_0: 1e-07 z_p_prime_base_10[3] + [
   - 0.5 alpha[0] * z_p_prime_base_10[3]
   - 0.5 alpha[1] * z_p_prime_base_10[3]
   - 0.5 alpha[2] * z_p_prime_base_10[3]
   - 0.5 alpha[3] * z_p_prime_base_10[3] - 0.5 p[0] * z_p_prime_base_10[3]
   - 0.5 p[1] * z_p_prime_base_10[3] + 0.5 p[2] * z_p_prime_base_10[3]
   + 0.5 p[3] * z_p_prime_base_10[3] ] <= 0
 infeasibility_p_base_11_row_0_leq_0: 1e-07 z_p_base_11[0] + [
   p[0] * z_p_base_11[0] ] <= 0
 infeasibility_p_base_11_row_1_leq_0: 1e-07 z_p_base_11[1] + [
   0.5 p[0] * z_p_base_11[1] + 0.5 p[1] * z_p_base_11[1]
   + 0.5 p[2] * z_p_base_11[1] - 0.5 p[3] * z_p_base_11[1] ] <= 0
 infeasibility_p_base_11_row_2_leq_0: 1e-07 z_p_base_11[2] + [
   - 0.5 p[0] * z_p_base_11[2] + 0.5 p[1] * z_p_base_11[2]
   - 0.5 p[2] * z_p_base_11[2] + 0.5 p[3] * z_p_base_11[2] ] <= 0
 infeasibility_p_base_11_row_3_leq_0: 1e-07 z_p_base_11[3] + [
   - 0.5 p[0] * z_p_base_11[3] - 0.5 p[1] * z_p_base_11[3]
   + 0.5 p[2] * z_p_base_11[3] + 0.5 p[3] * z_p_base_11[3] ] <= 0
 infeasibility_p_prime_base_11_row_0_leq_0: 1e-07 z_p_prime_base_11[0] + [
   alpha[0] * z_p_prime_base_11[0] + p[0] * z_p_prime_base_11[0] ] <= 0
 infeasibility_p_prime_base_11_row_1_leq_0: 1e-07 z_p_prime_base_11[1] + [
   0.5 alpha[0] * z_p_prime_base_11[1]
   + 0.5 alpha[1] * z_p_prime_base_11[1]
   - 0.5 alpha[2] * z_p_prime_base_11[1]
   + 0.5 alpha[3] * z_p_prime_base_11[1] + 0.5 p[0] * z_p_prime_base_11[1]
   + 0.5 p[1] * z_p_prime_base_11[1] + 0.5 p[2] * z_p_prime_base_11[1]
   - 0.5 p[3] * z_p_prime_base_11[1] ] <= 0
 infeasibility_p_prime_base_11_row_2_leq_0: 1e-07 z_p_prime_base_11[2] + [
   - 0.5 alpha[0] * z_p_prime_base_11[2]
   + 0.5 alpha[1] * z_p_prime_base_11[2]
   + 0.5 alpha[2] * z_p_prime_base_11[2]
   - 0.5 alpha[3] * z_p_prime_base_11[2] - 0.5 p[0] * z_p_prime_base_11[2]
   + 0.5 p[1] * z_p_prime_base_11[2] - 0.5 p[2] * z_p_prime_base_11[2]
   + 0.5 p[3] * z_p_prime_base_11[2] ] <= 0
 infeasibility_p_prime_base_11_row_3_leq_0: 1e-07 z_p_prime_base_11[3] + [
   - 0.5 alpha[0] * z_p_prime_base_11[3]
   - 0.5 alpha[1] * z_p_prime_base_11[3]
   - 0.5 alpha[2] * z_p_prime_base_11[3]
   - 0.5 alpha[3] * z_p_prime_base_11[3] - 0.5 p[0] * z_p_prime_base_11[3]
   - 0.5 p[1] * z_p_prime_base_11[3] + 0.5 p[2] * z_p_prime_base_11[3]
   + 0.5 p[3] * z_p_prime_base_11[3] ] <= 0
 infeasibility_p_base_B_prime_row_0_leq_0: 1e-07 z_p_base_B_prime[0] + [
   0.5 p[0] * z_p_base_B_prime[0] + 0.5 p[1] * z_p_base_B_prime[0]
   - 0.5 p[2] * z_p_base_B_prime[0] - 0.5 p[3] * z_p_base_B_prime[0] ]
   <= 0
 infeasibility_p_base_B_prime_row_1_leq_0: 1e-07 z_p_base_B_prime[1] + [
   0.5 p[0] * z_p_base_B_prime[1] - 0.5 p[1] * z_p_base_B_prime[1]
   + 0.5 p[2] * z_p_base_B_prime[1] - 0.5 p[3] * z_p_base_B_prime[1] ]
   <= 0
 infeasibility_p_base_B_prime_row_2_leq_0: 1e-07 z_p_base_B_prime[2] + [
   p[3] * z_p_base_B_prime[2] ] <= 0
 infeasibility_p_base_B_prime_row_3_leq_0: 1e-07 z_p_base_B_prime[3] + [
   - 0.5 p[0] * z_p_base_B_prime[3] + 0.5 p[1] * z_p_base_B_prime[3]
   + 0.5 p[2] * z_p_base_B_prime[3] + 0.5 p[3] * z_p_base_B_prime[3] ]
   <= 0
 infeasibility_p_prime_base_B_row_0_leq_0: 1e-07 z_p_prime_base_B[0] + [
   alpha[1] * z_p_prime_base_B[0] + p[1] * z_p_prime_base_B[0] ] <= 0
 infeasibility_p_prime_base_B_row_1_leq_0: 1e-07 z_p_prime_base_B[1] + [
   0.5 alpha[0] * z_p_prime_base_B[1] - 0.5 alpha[1] * z_p_prime_base_B[1]
   - 0.5 alpha[2] * z_p_prime_base_B[1]
   + 0.5 alpha[3] * z_p_prime_base_B[1] + 0.5 p[0] * z_p_prime_base_B[1]
   - 0.5 p[1] * z_p_prime_base_B[1] + 0.5 p[2] * z_p_prime_base_B[1]
   - 0.5 p[3] * z_p_prime_base_B[1] ] <= 0
 infeasibility_p_prime_base_B_row_2_leq_0: 1e-07 z_p_prime_base_B[2] + [
   0.5 alpha[0] * z_p_prime_base_B[2] - 0.5 alpha[1] * z_p_prime_base_B[2]
   + 0.5 alpha[2] * z_p_prime_base_B[2]
   - 0.5 alpha[3] * z_p_prime_base_B[2] + 0.5 p[0] * z_p_prime_base_B[2]
   - 0.5 p[1] * z_p_prime_base_B[2] - 0.5 p[2] * z_p_prime_base_B[2]
   + 0.5 p[3] * z_p_prime_base_B[2] ] <= 0
 infeasibility_p_prime_base_B_row_3_leq_0: 1e-07 z_p_prime_base_B[3] + [
   - 0.5 alpha[0] * z_p_prime_base_B[3]
   + 0.5 alpha[1] * z_p_prime_base_B[3]
   - 0.5 alpha[2] * z_p_prime_base_B[3]
   - 0.5 alpha[3] * z_p_prime_base_B[3] - 0.5 p[0] * z_p_prime_base_B[3]
   + 0.5 p[1] * z_p_prime_base_B[3] + 0.5 p[2] * z_p_prime_base_B[3]
   + 0.5 p[3] * z_p_prime_base_B[3] ] <= 0
Bounds
 alpha[0] <= 1
 alpha[1] <= 1
 alpha[2] <= 1
 alpha[3] <= 1
 1e-07 <= p[0] <= 0.9999999000000001
 1e-07 <= p[1] <= 0.9999999000000001
 1e-07 <= p[2] <= 0.9999999000000001
 1e-07 <= p[3] <= 0.9999999000000001
Binaries
 z_p_base_1[0] z_p_base_1[1] z_p_base_1[2] z_p_base_1[3]
 z_p_prime_base_1[0] z_p_prime_base_1[1] z_p_prime_base_1[2]
 z_p_prime_base_1[3] z_p_base_3[0] z_p_base_3[1] z_p_base_3[2]
 z_p_base_3[3] z_p_prime_base_3[0] z_p_prime_base_3[1] z_p_prime_base_3[2]
 z_p_prime_base_3[3] z_p_base_4[0] z_p_base_4[1] z_p_base_4[2]
 z_p_base_4[3] z_p_prime_base_4[0] z_p_prime_base_4[1] z_p_prime_base_4[2]
 z_p_prime_base_4[3] z_p_base_5[0] z_p_base_5[1] z_p_base_5[2]
 z_p_base_5[3] z_p_prime_base_5[0] z_p_prime_base_5[1] z_p_prime_base_5[2]
 z_p_prime_base_5[3] z_p_base_6[0] z_p_base_6[1] z_p_base_6[2]
 z_p_base_6[3] z_p_prime_base_6[0] z_p_prime_base_6[1] z_p_prime_base_6[2]
 z_p_prime_base_6[3] z_p_base_7[0] z_p_base_7[1] z_p_base_7[2]
 z_p_base_7[3] z_p_prime_base_7[0] z_p_prime_base_7[1] z_p_prime_base_7[2]
 z_p_prime_base_7[3] z_p_base_8[0] z_p_base_8[1] z_p_base_8[2]
 z_p_base_8[3] z_p_prime_base_8[0] z_p_prime_base_8[1] z_p_prime_base_8[2]
 z_p_prime_base_8[3] z_p_base_9[0] z_p_base_9[1] z_p_base_9[2]
 z_p_base_9[3] z_p_prime_base_9[0] z_p_prime_base_9[1] z_p_prime_base_9[2]
 z_p_prime_base_9[3] z_p_base_10[0] z_p_base_10[1] z_p_base_10[2]
 z_p_base_10[3] z_p_prime_base_10[0] z_p_prime_base_10[1]
 z_p_prime_base_10[2] z_p_prime_base_10[3] z_p_base_11[0] z_p_base_11[1]
 z_p_base_11[2] z_p_base_11[3] z_p_prime_base_11[0] z_p_prime_base_11[1]
 z_p_prime_base_11[2] z_p_prime_base_11[3] z_p_base_B_prime[0]
 z_p_base_B_prime[1] z_p_base_B_prime[2] z_p_base_B_prime[3]
 z_p_prime_base_B[0] z_p_prime_base_B[1] z_p_prime_base_B[2]
 z_p_prime_base_B[3]
End
