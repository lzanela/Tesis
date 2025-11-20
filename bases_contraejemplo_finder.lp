\ Model contraejemplo_finder
\ LP format - for model browsing. Use MPS format to capture full model detail.
\ Signature: 0x8bd5eb3276e07055
Maximize
  alpha[0] + alpha[4]
Subject To
 p_prime_0_leq_1: alpha[0] + p[0] <= 1
 p_prime_0_geq_0: alpha[0] + p[0] >= 0
 p_prime_4_leq_1: alpha[4] + p[4] <= 1
 p_prime_4_geq_0: alpha[4] + p[4] >= 0
 p_prime_1_leq_1: - alpha[1] + p[1] <= 1
 p_prime_1_geq_0: - alpha[1] + p[1] >= 0
 p_prime_2_leq_1: - alpha[2] + p[2] <= 1
 p_prime_2_geq_0: - alpha[2] + p[2] >= 0
 p_prime_3_leq_1: - alpha[3] + p[3] <= 1
 p_prime_3_geq_0: - alpha[3] + p[3] >= 0
 p_suma_K: p[0] + p[1] + p[2] + p[3] + p[4] = 2
 alpha_suma_0: alpha[0] - alpha[1] - alpha[2] - alpha[3] + alpha[4] = 0
 feasibility_p_geq_0_row_0: 0.5 p[0] + 0.5 p[1] - 0.5 p[2] - 0.5 p[3]
   - 0.5 p[4] >= 0
 feasibility_p_leq_1_row_0: 0.5 p[0] + 0.5 p[1] - 0.5 p[2] - 0.5 p[3]
   - 0.5 p[4] <= 1
 feasibility_p_prime_geq_0_row_0: 0.5 alpha[0] - 0.5 alpha[1]
   - 0.5 alpha[2] + 0.5 alpha[3] - 0.5 alpha[4] + 0.5 p[0] + 0.5 p[1]
   + 0.5 p[2] - 0.5 p[3] - 0.5 p[4] >= 0
 feasibility_p_prime_leq_1_row_0: 0.5 alpha[0] - 0.5 alpha[1]
   - 0.5 alpha[2] + 0.5 alpha[3] - 0.5 alpha[4] + 0.5 p[0] + 0.5 p[1]
   + 0.5 p[2] - 0.5 p[3] - 0.5 p[4] <= 1
 feasibility_p_geq_0_row_1: 0.5 p[0] - 0.5 p[1] + 0.5 p[2] - 0.5 p[3]
   - 0.5 p[4] >= 0
 feasibility_p_leq_1_row_1: 0.5 p[0] - 0.5 p[1] + 0.5 p[2] - 0.5 p[3]
   - 0.5 p[4] <= 1
 feasibility_p_prime_geq_0_row_1: 0.5 alpha[0] + 0.5 alpha[1]
   + 0.5 alpha[2] - 0.5 alpha[3] + 0.5 alpha[4] + 0.5 p[0] - 0.5 p[1]
   - 0.5 p[2] + 0.5 p[3] + 0.5 p[4] >= 0
 feasibility_p_prime_leq_1_row_1: 0.5 alpha[0] + 0.5 alpha[1]
   + 0.5 alpha[2] - 0.5 alpha[3] + 0.5 alpha[4] + 0.5 p[0] - 0.5 p[1]
   - 0.5 p[2] + 0.5 p[3] + 0.5 p[4] <= 1
 feasibility_p_geq_0_row_2: p[3] >= 0
 feasibility_p_leq_1_row_2: p[3] <= 1
 feasibility_p_prime_geq_0_row_2: - 0.5 alpha[0] - 0.5 alpha[1]
   + 0.5 alpha[2] - 0.5 alpha[3] + 0.5 alpha[4] - 0.5 p[0] + 0.5 p[1]
   - 0.5 p[2] + 0.5 p[3] + 0.5 p[4] >= 0
 feasibility_p_prime_leq_1_row_2: - 0.5 alpha[0] - 0.5 alpha[1]
   + 0.5 alpha[2] - 0.5 alpha[3] + 0.5 alpha[4] - 0.5 p[0] + 0.5 p[1]
   - 0.5 p[2] + 0.5 p[3] + 0.5 p[4] <= 1
 feasibility_p_geq_0_row_3: p[4] >= 0
 feasibility_p_leq_1_row_3: p[4] <= 1
 feasibility_p_prime_geq_0_row_3: - alpha[3] + p[3] >= 0
 feasibility_p_prime_leq_1_row_3: - alpha[3] + p[3] <= 1
 feasibility_p_geq_0_row_4: - 0.5 p[0] + 0.5 p[1] + 0.5 p[2] + 0.5 p[3]
   + 0.5 p[4] >= 0
 feasibility_p_leq_1_row_4: - 0.5 p[0] + 0.5 p[1] + 0.5 p[2] + 0.5 p[3]
   + 0.5 p[4] <= 1
 feasibility_p_prime_geq_0_row_4: - alpha[2] + alpha[3] + p[2] - p[3] >= 0
 feasibility_p_prime_leq_1_row_4: - alpha[2] + alpha[3] + p[2] - p[3] <= 1
 infeasibility_p_base_1_at_least_one_violated: z_p_base_1[0]
   + z_p_base_1[1] + z_p_base_1[2] + z_p_base_1[3] + z_p_base_1[4] = 1
 infeasibility_p_prime_base_1_at_least_one_violated: z_p_prime_base_1[0]
   + z_p_prime_base_1[1] + z_p_prime_base_1[2] + z_p_prime_base_1[3]
   + z_p_prime_base_1[4] = 1
 infeasibility_p_base_2_at_least_one_violated: z_p_base_2[0]
   + z_p_base_2[1] + z_p_base_2[2] + z_p_base_2[3] + z_p_base_2[4] = 1
 infeasibility_p_prime_base_2_at_least_one_violated: z_p_prime_base_2[0]
   + z_p_prime_base_2[1] + z_p_prime_base_2[2] + z_p_prime_base_2[3]
   + z_p_prime_base_2[4] = 1
 infeasibility_p_base_3_at_least_one_violated: z_p_base_3[0]
   + z_p_base_3[1] + z_p_base_3[2] + z_p_base_3[3] + z_p_base_3[4] = 1
 infeasibility_p_prime_base_3_at_least_one_violated: z_p_prime_base_3[0]
   + z_p_prime_base_3[1] + z_p_prime_base_3[2] + z_p_prime_base_3[3]
   + z_p_prime_base_3[4] = 1
 infeasibility_p_base_4_at_least_one_violated: z_p_base_4[0]
   + z_p_base_4[1] + z_p_base_4[2] + z_p_base_4[3] + z_p_base_4[4] = 1
 infeasibility_p_prime_base_4_at_least_one_violated: z_p_prime_base_4[0]
   + z_p_prime_base_4[1] + z_p_prime_base_4[2] + z_p_prime_base_4[3]
   + z_p_prime_base_4[4] = 1
 infeasibility_p_base_5_at_least_one_violated: z_p_base_5[0]
   + z_p_base_5[1] + z_p_base_5[2] + z_p_base_5[3] + z_p_base_5[4] = 1
 infeasibility_p_prime_base_5_at_least_one_violated: z_p_prime_base_5[0]
   + z_p_prime_base_5[1] + z_p_prime_base_5[2] + z_p_prime_base_5[3]
   + z_p_prime_base_5[4] = 1
 infeasibility_p_base_6_at_least_one_violated: z_p_base_6[0]
   + z_p_base_6[1] + z_p_base_6[2] + z_p_base_6[3] + z_p_base_6[4] = 1
 infeasibility_p_prime_base_6_at_least_one_violated: z_p_prime_base_6[0]
   + z_p_prime_base_6[1] + z_p_prime_base_6[2] + z_p_prime_base_6[3]
   + z_p_prime_base_6[4] = 1
 infeasibility_p_base_7_at_least_one_violated: z_p_base_7[0]
   + z_p_base_7[1] + z_p_base_7[2] + z_p_base_7[3] + z_p_base_7[4] = 1
 infeasibility_p_prime_base_7_at_least_one_violated: z_p_prime_base_7[0]
   + z_p_prime_base_7[1] + z_p_prime_base_7[2] + z_p_prime_base_7[3]
   + z_p_prime_base_7[4] = 1
 infeasibility_p_base_8_at_least_one_violated: z_p_base_8[0]
   + z_p_base_8[1] + z_p_base_8[2] + z_p_base_8[3] + z_p_base_8[4] = 1
 infeasibility_p_prime_base_8_at_least_one_violated: z_p_prime_base_8[0]
   + z_p_prime_base_8[1] + z_p_prime_base_8[2] + z_p_prime_base_8[3]
   + z_p_prime_base_8[4] = 1
 infeasibility_p_base_9_at_least_one_violated: z_p_base_9[0]
   + z_p_base_9[1] + z_p_base_9[2] + z_p_base_9[3] + z_p_base_9[4] = 1
 infeasibility_p_prime_base_9_at_least_one_violated: z_p_prime_base_9[0]
   + z_p_prime_base_9[1] + z_p_prime_base_9[2] + z_p_prime_base_9[3]
   + z_p_prime_base_9[4] = 1
 infeasibility_p_base_10_at_least_one_violated: z_p_base_10[0]
   + z_p_base_10[1] + z_p_base_10[2] + z_p_base_10[3] + z_p_base_10[4] = 1
 infeasibility_p_prime_base_10_at_least_one_violated: z_p_prime_base_10[0]
   + z_p_prime_base_10[1] + z_p_prime_base_10[2] + z_p_prime_base_10[3]
   + z_p_prime_base_10[4] = 1
 infeasibility_p_base_11_at_least_one_violated: z_p_base_11[0]
   + z_p_base_11[1] + z_p_base_11[2] + z_p_base_11[3] + z_p_base_11[4] = 1
 infeasibility_p_prime_base_11_at_least_one_violated: z_p_prime_base_11[0]
   + z_p_prime_base_11[1] + z_p_prime_base_11[2] + z_p_prime_base_11[3]
   + z_p_prime_base_11[4] = 1
 infeasibility_p_base_12_at_least_one_violated: z_p_base_12[0]
   + z_p_base_12[1] + z_p_base_12[2] + z_p_base_12[3] + z_p_base_12[4] = 1
 infeasibility_p_prime_base_12_at_least_one_violated: z_p_prime_base_12[0]
   + z_p_prime_base_12[1] + z_p_prime_base_12[2] + z_p_prime_base_12[3]
   + z_p_prime_base_12[4] = 1
 infeasibility_p_base_13_at_least_one_violated: z_p_base_13[0]
   + z_p_base_13[1] + z_p_base_13[2] + z_p_base_13[3] + z_p_base_13[4] = 1
 infeasibility_p_prime_base_13_at_least_one_violated: z_p_prime_base_13[0]
   + z_p_prime_base_13[1] + z_p_prime_base_13[2] + z_p_prime_base_13[3]
   + z_p_prime_base_13[4] = 1
 infeasibility_p_base_14_at_least_one_violated: z_p_base_14[0]
   + z_p_base_14[1] + z_p_base_14[2] + z_p_base_14[3] + z_p_base_14[4] = 1
 infeasibility_p_prime_base_14_at_least_one_violated: z_p_prime_base_14[0]
   + z_p_prime_base_14[1] + z_p_prime_base_14[2] + z_p_prime_base_14[3]
   + z_p_prime_base_14[4] = 1
 infeasibility_p_base_15_at_least_one_violated: z_p_base_15[0]
   + z_p_base_15[1] + z_p_base_15[2] + z_p_base_15[3] + z_p_base_15[4] = 1
 infeasibility_p_prime_base_15_at_least_one_violated: z_p_prime_base_15[0]
   + z_p_prime_base_15[1] + z_p_prime_base_15[2] + z_p_prime_base_15[3]
   + z_p_prime_base_15[4] = 1
 infeasibility_p_base_16_at_least_one_violated: z_p_base_16[0]
   + z_p_base_16[1] + z_p_base_16[2] + z_p_base_16[3] + z_p_base_16[4] = 1
 infeasibility_p_prime_base_16_at_least_one_violated: z_p_prime_base_16[0]
   + z_p_prime_base_16[1] + z_p_prime_base_16[2] + z_p_prime_base_16[3]
   + z_p_prime_base_16[4] = 1
 infeasibility_p_base_17_at_least_one_violated: z_p_base_17[0]
   + z_p_base_17[1] + z_p_base_17[2] + z_p_base_17[3] + z_p_base_17[4] = 1
 infeasibility_p_prime_base_17_at_least_one_violated: z_p_prime_base_17[0]
   + z_p_prime_base_17[1] + z_p_prime_base_17[2] + z_p_prime_base_17[3]
   + z_p_prime_base_17[4] = 1
 infeasibility_p_base_18_at_least_one_violated: z_p_base_18[0]
   + z_p_base_18[1] + z_p_base_18[2] + z_p_base_18[3] + z_p_base_18[4] = 1
 infeasibility_p_prime_base_18_at_least_one_violated: z_p_prime_base_18[0]
   + z_p_prime_base_18[1] + z_p_prime_base_18[2] + z_p_prime_base_18[3]
   + z_p_prime_base_18[4] = 1
 infeasibility_p_base_19_at_least_one_violated: z_p_base_19[0]
   + z_p_base_19[1] + z_p_base_19[2] + z_p_base_19[3] + z_p_base_19[4] = 1
 infeasibility_p_prime_base_19_at_least_one_violated: z_p_prime_base_19[0]
   + z_p_prime_base_19[1] + z_p_prime_base_19[2] + z_p_prime_base_19[3]
   + z_p_prime_base_19[4] = 1
 infeasibility_p_base_20_at_least_one_violated: z_p_base_20[0]
   + z_p_base_20[1] + z_p_base_20[2] + z_p_base_20[3] + z_p_base_20[4] = 1
 infeasibility_p_prime_base_20_at_least_one_violated: z_p_prime_base_20[0]
   + z_p_prime_base_20[1] + z_p_prime_base_20[2] + z_p_prime_base_20[3]
   + z_p_prime_base_20[4] = 1
 infeasibility_p_base_21_at_least_one_violated: z_p_base_21[0]
   + z_p_base_21[1] + z_p_base_21[2] + z_p_base_21[3] + z_p_base_21[4] = 1
 infeasibility_p_prime_base_21_at_least_one_violated: z_p_prime_base_21[0]
   + z_p_prime_base_21[1] + z_p_prime_base_21[2] + z_p_prime_base_21[3]
   + z_p_prime_base_21[4] = 1
 infeasibility_p_base_22_at_least_one_violated: z_p_base_22[0]
   + z_p_base_22[1] + z_p_base_22[2] + z_p_base_22[3] + z_p_base_22[4] = 1
 infeasibility_p_prime_base_22_at_least_one_violated: z_p_prime_base_22[0]
   + z_p_prime_base_22[1] + z_p_prime_base_22[2] + z_p_prime_base_22[3]
   + z_p_prime_base_22[4] = 1
 infeasibility_p_base_23_at_least_one_violated: z_p_base_23[0]
   + z_p_base_23[1] + z_p_base_23[2] + z_p_base_23[3] + z_p_base_23[4] = 1
 infeasibility_p_prime_base_23_at_least_one_violated: z_p_prime_base_23[0]
   + z_p_prime_base_23[1] + z_p_prime_base_23[2] + z_p_prime_base_23[3]
   + z_p_prime_base_23[4] = 1
 infeasibility_p_base_24_at_least_one_violated: z_p_base_24[0]
   + z_p_base_24[1] + z_p_base_24[2] + z_p_base_24[3] + z_p_base_24[4] = 1
 infeasibility_p_prime_base_24_at_least_one_violated: z_p_prime_base_24[0]
   + z_p_prime_base_24[1] + z_p_prime_base_24[2] + z_p_prime_base_24[3]
   + z_p_prime_base_24[4] = 1
 infeasibility_p_base_25_at_least_one_violated: z_p_base_25[0]
   + z_p_base_25[1] + z_p_base_25[2] + z_p_base_25[3] + z_p_base_25[4] = 1
 infeasibility_p_prime_base_25_at_least_one_violated: z_p_prime_base_25[0]
   + z_p_prime_base_25[1] + z_p_prime_base_25[2] + z_p_prime_base_25[3]
   + z_p_prime_base_25[4] = 1
 infeasibility_p_base_26_at_least_one_violated: z_p_base_26[0]
   + z_p_base_26[1] + z_p_base_26[2] + z_p_base_26[3] + z_p_base_26[4] = 1
 infeasibility_p_prime_base_26_at_least_one_violated: z_p_prime_base_26[0]
   + z_p_prime_base_26[1] + z_p_prime_base_26[2] + z_p_prime_base_26[3]
   + z_p_prime_base_26[4] = 1
 infeasibility_p_base_27_at_least_one_violated: z_p_base_27[0]
   + z_p_base_27[1] + z_p_base_27[2] + z_p_base_27[3] + z_p_base_27[4] = 1
 infeasibility_p_prime_base_27_at_least_one_violated: z_p_prime_base_27[0]
   + z_p_prime_base_27[1] + z_p_prime_base_27[2] + z_p_prime_base_27[3]
   + z_p_prime_base_27[4] = 1
 infeasibility_p_base_28_at_least_one_violated: z_p_base_28[0]
   + z_p_base_28[1] + z_p_base_28[2] + z_p_base_28[3] + z_p_base_28[4] = 1
 infeasibility_p_prime_base_28_at_least_one_violated: z_p_prime_base_28[0]
   + z_p_prime_base_28[1] + z_p_prime_base_28[2] + z_p_prime_base_28[3]
   + z_p_prime_base_28[4] = 1
 infeasibility_p_base_29_at_least_one_violated: z_p_base_29[0]
   + z_p_base_29[1] + z_p_base_29[2] + z_p_base_29[3] + z_p_base_29[4] = 1
 infeasibility_p_prime_base_29_at_least_one_violated: z_p_prime_base_29[0]
   + z_p_prime_base_29[1] + z_p_prime_base_29[2] + z_p_prime_base_29[3]
   + z_p_prime_base_29[4] = 1
 infeasibility_p_base_30_at_least_one_violated: z_p_base_30[0]
   + z_p_base_30[1] + z_p_base_30[2] + z_p_base_30[3] + z_p_base_30[4] = 1
 infeasibility_p_prime_base_30_at_least_one_violated: z_p_prime_base_30[0]
   + z_p_prime_base_30[1] + z_p_prime_base_30[2] + z_p_prime_base_30[3]
   + z_p_prime_base_30[4] = 1
 infeasibility_p_base_31_at_least_one_violated: z_p_base_31[0]
   + z_p_base_31[1] + z_p_base_31[2] + z_p_base_31[3] + z_p_base_31[4] = 1
 infeasibility_p_prime_base_31_at_least_one_violated: z_p_prime_base_31[0]
   + z_p_prime_base_31[1] + z_p_prime_base_31[2] + z_p_prime_base_31[3]
   + z_p_prime_base_31[4] = 1
 infeasibility_p_base_32_at_least_one_violated: z_p_base_32[0]
   + z_p_base_32[1] + z_p_base_32[2] + z_p_base_32[3] + z_p_base_32[4] = 1
 infeasibility_p_prime_base_32_at_least_one_violated: z_p_prime_base_32[0]
   + z_p_prime_base_32[1] + z_p_prime_base_32[2] + z_p_prime_base_32[3]
   + z_p_prime_base_32[4] = 1
 infeasibility_p_base_33_at_least_one_violated: z_p_base_33[0]
   + z_p_base_33[1] + z_p_base_33[2] + z_p_base_33[3] + z_p_base_33[4] = 1
 infeasibility_p_prime_base_33_at_least_one_violated: z_p_prime_base_33[0]
   + z_p_prime_base_33[1] + z_p_prime_base_33[2] + z_p_prime_base_33[3]
   + z_p_prime_base_33[4] = 1
 infeasibility_p_base_34_at_least_one_violated: z_p_base_34[0]
   + z_p_base_34[1] + z_p_base_34[2] + z_p_base_34[3] + z_p_base_34[4] = 1
 infeasibility_p_prime_base_34_at_least_one_violated: z_p_prime_base_34[0]
   + z_p_prime_base_34[1] + z_p_prime_base_34[2] + z_p_prime_base_34[3]
   + z_p_prime_base_34[4] = 1
 infeasibility_p_base_35_at_least_one_violated: z_p_base_35[0]
   + z_p_base_35[1] + z_p_base_35[2] + z_p_base_35[3] + z_p_base_35[4] = 1
 infeasibility_p_prime_base_35_at_least_one_violated: z_p_prime_base_35[0]
   + z_p_prime_base_35[1] + z_p_prime_base_35[2] + z_p_prime_base_35[3]
   + z_p_prime_base_35[4] = 1
 infeasibility_p_base_36_at_least_one_violated: z_p_base_36[0]
   + z_p_base_36[1] + z_p_base_36[2] + z_p_base_36[3] + z_p_base_36[4] = 1
 infeasibility_p_prime_base_36_at_least_one_violated: z_p_prime_base_36[0]
   + z_p_prime_base_36[1] + z_p_prime_base_36[2] + z_p_prime_base_36[3]
   + z_p_prime_base_36[4] = 1
 infeasibility_p_base_37_at_least_one_violated: z_p_base_37[0]
   + z_p_base_37[1] + z_p_base_37[2] + z_p_base_37[3] + z_p_base_37[4] = 1
 infeasibility_p_prime_base_37_at_least_one_violated: z_p_prime_base_37[0]
   + z_p_prime_base_37[1] + z_p_prime_base_37[2] + z_p_prime_base_37[3]
   + z_p_prime_base_37[4] = 1
 infeasibility_p_base_38_at_least_one_violated: z_p_base_38[0]
   + z_p_base_38[1] + z_p_base_38[2] + z_p_base_38[3] + z_p_base_38[4] = 1
 infeasibility_p_prime_base_38_at_least_one_violated: z_p_prime_base_38[0]
   + z_p_prime_base_38[1] + z_p_prime_base_38[2] + z_p_prime_base_38[3]
   + z_p_prime_base_38[4] = 1
 infeasibility_p_base_39_at_least_one_violated: z_p_base_39[0]
   + z_p_base_39[1] + z_p_base_39[2] + z_p_base_39[3] + z_p_base_39[4] = 1
 infeasibility_p_prime_base_39_at_least_one_violated: z_p_prime_base_39[0]
   + z_p_prime_base_39[1] + z_p_prime_base_39[2] + z_p_prime_base_39[3]
   + z_p_prime_base_39[4] = 1
 infeasibility_p_base_40_at_least_one_violated: z_p_base_40[0]
   + z_p_base_40[1] + z_p_base_40[2] + z_p_base_40[3] + z_p_base_40[4] = 1
 infeasibility_p_prime_base_40_at_least_one_violated: z_p_prime_base_40[0]
   + z_p_prime_base_40[1] + z_p_prime_base_40[2] + z_p_prime_base_40[3]
   + z_p_prime_base_40[4] = 1
 infeasibility_p_base_41_at_least_one_violated: z_p_base_41[0]
   + z_p_base_41[1] + z_p_base_41[2] + z_p_base_41[3] + z_p_base_41[4] = 1
 infeasibility_p_prime_base_41_at_least_one_violated: z_p_prime_base_41[0]
   + z_p_prime_base_41[1] + z_p_prime_base_41[2] + z_p_prime_base_41[3]
   + z_p_prime_base_41[4] = 1
 infeasibility_p_base_42_at_least_one_violated: z_p_base_42[0]
   + z_p_base_42[1] + z_p_base_42[2] + z_p_base_42[3] + z_p_base_42[4] = 1
 infeasibility_p_prime_base_42_at_least_one_violated: z_p_prime_base_42[0]
   + z_p_prime_base_42[1] + z_p_prime_base_42[2] + z_p_prime_base_42[3]
   + z_p_prime_base_42[4] = 1
 infeasibility_p_base_43_at_least_one_violated: z_p_base_43[0]
   + z_p_base_43[1] + z_p_base_43[2] + z_p_base_43[3] + z_p_base_43[4] = 1
 infeasibility_p_prime_base_43_at_least_one_violated: z_p_prime_base_43[0]
   + z_p_prime_base_43[1] + z_p_prime_base_43[2] + z_p_prime_base_43[3]
   + z_p_prime_base_43[4] = 1
 infeasibility_p_base_44_at_least_one_violated: z_p_base_44[0]
   + z_p_base_44[1] + z_p_base_44[2] + z_p_base_44[3] + z_p_base_44[4] = 1
 infeasibility_p_prime_base_44_at_least_one_violated: z_p_prime_base_44[0]
   + z_p_prime_base_44[1] + z_p_prime_base_44[2] + z_p_prime_base_44[3]
   + z_p_prime_base_44[4] = 1
 infeasibility_p_base_45_at_least_one_violated: z_p_base_45[0]
   + z_p_base_45[1] + z_p_base_45[2] + z_p_base_45[3] + z_p_base_45[4] = 1
 infeasibility_p_prime_base_45_at_least_one_violated: z_p_prime_base_45[0]
   + z_p_prime_base_45[1] + z_p_prime_base_45[2] + z_p_prime_base_45[3]
   + z_p_prime_base_45[4] = 1
 infeasibility_p_base_46_at_least_one_violated: z_p_base_46[0]
   + z_p_base_46[1] + z_p_base_46[2] + z_p_base_46[3] + z_p_base_46[4] = 1
 infeasibility_p_prime_base_46_at_least_one_violated: z_p_prime_base_46[0]
   + z_p_prime_base_46[1] + z_p_prime_base_46[2] + z_p_prime_base_46[3]
   + z_p_prime_base_46[4] = 1
 infeasibility_p_base_47_at_least_one_violated: z_p_base_47[0]
   + z_p_base_47[1] + z_p_base_47[2] + z_p_base_47[3] + z_p_base_47[4] = 1
 infeasibility_p_prime_base_47_at_least_one_violated: z_p_prime_base_47[0]
   + z_p_prime_base_47[1] + z_p_prime_base_47[2] + z_p_prime_base_47[3]
   + z_p_prime_base_47[4] = 1
 infeasibility_p_base_48_at_least_one_violated: z_p_base_48[0]
   + z_p_base_48[1] + z_p_base_48[2] + z_p_base_48[3] + z_p_base_48[4] = 1
 infeasibility_p_prime_base_48_at_least_one_violated: z_p_prime_base_48[0]
   + z_p_prime_base_48[1] + z_p_prime_base_48[2] + z_p_prime_base_48[3]
   + z_p_prime_base_48[4] = 1
 infeasibility_p_base_49_at_least_one_violated: z_p_base_49[0]
   + z_p_base_49[1] + z_p_base_49[2] + z_p_base_49[3] + z_p_base_49[4] = 1
 infeasibility_p_prime_base_49_at_least_one_violated: z_p_prime_base_49[0]
   + z_p_prime_base_49[1] + z_p_prime_base_49[2] + z_p_prime_base_49[3]
   + z_p_prime_base_49[4] = 1
 infeasibility_p_base_50_at_least_one_violated: z_p_base_50[0]
   + z_p_base_50[1] + z_p_base_50[2] + z_p_base_50[3] + z_p_base_50[4] = 1
 infeasibility_p_prime_base_50_at_least_one_violated: z_p_prime_base_50[0]
   + z_p_prime_base_50[1] + z_p_prime_base_50[2] + z_p_prime_base_50[3]
   + z_p_prime_base_50[4] = 1
 infeasibility_p_base_51_at_least_one_violated: z_p_base_51[0]
   + z_p_base_51[1] + z_p_base_51[2] + z_p_base_51[3] + z_p_base_51[4] = 1
 infeasibility_p_prime_base_51_at_least_one_violated: z_p_prime_base_51[0]
   + z_p_prime_base_51[1] + z_p_prime_base_51[2] + z_p_prime_base_51[3]
   + z_p_prime_base_51[4] = 1
 infeasibility_p_base_52_at_least_one_violated: z_p_base_52[0]
   + z_p_base_52[1] + z_p_base_52[2] + z_p_base_52[3] + z_p_base_52[4] = 1
 infeasibility_p_prime_base_52_at_least_one_violated: z_p_prime_base_52[0]
   + z_p_prime_base_52[1] + z_p_prime_base_52[2] + z_p_prime_base_52[3]
   + z_p_prime_base_52[4] = 1
 infeasibility_p_base_53_at_least_one_violated: z_p_base_53[0]
   + z_p_base_53[1] + z_p_base_53[2] + z_p_base_53[3] + z_p_base_53[4] = 1
 infeasibility_p_prime_base_53_at_least_one_violated: z_p_prime_base_53[0]
   + z_p_prime_base_53[1] + z_p_prime_base_53[2] + z_p_prime_base_53[3]
   + z_p_prime_base_53[4] = 1
 infeasibility_p_base_54_at_least_one_violated: z_p_base_54[0]
   + z_p_base_54[1] + z_p_base_54[2] + z_p_base_54[3] + z_p_base_54[4] = 1
 infeasibility_p_prime_base_54_at_least_one_violated: z_p_prime_base_54[0]
   + z_p_prime_base_54[1] + z_p_prime_base_54[2] + z_p_prime_base_54[3]
   + z_p_prime_base_54[4] = 1
 infeasibility_p_base_55_at_least_one_violated: z_p_base_55[0]
   + z_p_base_55[1] + z_p_base_55[2] + z_p_base_55[3] + z_p_base_55[4] = 1
 infeasibility_p_prime_base_55_at_least_one_violated: z_p_prime_base_55[0]
   + z_p_prime_base_55[1] + z_p_prime_base_55[2] + z_p_prime_base_55[3]
   + z_p_prime_base_55[4] = 1
 infeasibility_p_base_56_at_least_one_violated: z_p_base_56[0]
   + z_p_base_56[1] + z_p_base_56[2] + z_p_base_56[3] + z_p_base_56[4] = 1
 infeasibility_p_prime_base_56_at_least_one_violated: z_p_prime_base_56[0]
   + z_p_prime_base_56[1] + z_p_prime_base_56[2] + z_p_prime_base_56[3]
   + z_p_prime_base_56[4] = 1
 infeasibility_p_base_57_at_least_one_violated: z_p_base_57[0]
   + z_p_base_57[1] + z_p_base_57[2] + z_p_base_57[3] + z_p_base_57[4] = 1
 infeasibility_p_prime_base_57_at_least_one_violated: z_p_prime_base_57[0]
   + z_p_prime_base_57[1] + z_p_prime_base_57[2] + z_p_prime_base_57[3]
   + z_p_prime_base_57[4] = 1
 infeasibility_p_base_58_at_least_one_violated: z_p_base_58[0]
   + z_p_base_58[1] + z_p_base_58[2] + z_p_base_58[3] + z_p_base_58[4] = 1
 infeasibility_p_prime_base_58_at_least_one_violated: z_p_prime_base_58[0]
   + z_p_prime_base_58[1] + z_p_prime_base_58[2] + z_p_prime_base_58[3]
   + z_p_prime_base_58[4] = 1
 infeasibility_p_base_59_at_least_one_violated: z_p_base_59[0]
   + z_p_base_59[1] + z_p_base_59[2] + z_p_base_59[3] + z_p_base_59[4] = 1
 infeasibility_p_prime_base_59_at_least_one_violated: z_p_prime_base_59[0]
   + z_p_prime_base_59[1] + z_p_prime_base_59[2] + z_p_prime_base_59[3]
   + z_p_prime_base_59[4] = 1
 infeasibility_p_base_60_at_least_one_violated: z_p_base_60[0]
   + z_p_base_60[1] + z_p_base_60[2] + z_p_base_60[3] + z_p_base_60[4] = 1
 infeasibility_p_prime_base_60_at_least_one_violated: z_p_prime_base_60[0]
   + z_p_prime_base_60[1] + z_p_prime_base_60[2] + z_p_prime_base_60[3]
   + z_p_prime_base_60[4] = 1
 infeasibility_p_base_61_at_least_one_violated: z_p_base_61[0]
   + z_p_base_61[1] + z_p_base_61[2] + z_p_base_61[3] + z_p_base_61[4] = 1
 infeasibility_p_prime_base_61_at_least_one_violated: z_p_prime_base_61[0]
   + z_p_prime_base_61[1] + z_p_prime_base_61[2] + z_p_prime_base_61[3]
   + z_p_prime_base_61[4] = 1
 infeasibility_p_base_62_at_least_one_violated: z_p_base_62[0]
   + z_p_base_62[1] + z_p_base_62[2] + z_p_base_62[3] + z_p_base_62[4] = 1
 infeasibility_p_prime_base_62_at_least_one_violated: z_p_prime_base_62[0]
   + z_p_prime_base_62[1] + z_p_prime_base_62[2] + z_p_prime_base_62[3]
   + z_p_prime_base_62[4] = 1
 infeasibility_p_base_63_at_least_one_violated: z_p_base_63[0]
   + z_p_base_63[1] + z_p_base_63[2] + z_p_base_63[3] + z_p_base_63[4] = 1
 infeasibility_p_prime_base_63_at_least_one_violated: z_p_prime_base_63[0]
   + z_p_prime_base_63[1] + z_p_prime_base_63[2] + z_p_prime_base_63[3]
   + z_p_prime_base_63[4] = 1
 infeasibility_p_base_64_at_least_one_violated: z_p_base_64[0]
   + z_p_base_64[1] + z_p_base_64[2] + z_p_base_64[3] + z_p_base_64[4] = 1
 infeasibility_p_prime_base_64_at_least_one_violated: z_p_prime_base_64[0]
   + z_p_prime_base_64[1] + z_p_prime_base_64[2] + z_p_prime_base_64[3]
   + z_p_prime_base_64[4] = 1
 infeasibility_p_base_66_at_least_one_violated: z_p_base_66[0]
   + z_p_base_66[1] + z_p_base_66[2] + z_p_base_66[3] + z_p_base_66[4] = 1
 infeasibility_p_prime_base_66_at_least_one_violated: z_p_prime_base_66[0]
   + z_p_prime_base_66[1] + z_p_prime_base_66[2] + z_p_prime_base_66[3]
   + z_p_prime_base_66[4] = 1
 infeasibility_p_base_67_at_least_one_violated: z_p_base_67[0]
   + z_p_base_67[1] + z_p_base_67[2] + z_p_base_67[3] + z_p_base_67[4] = 1
 infeasibility_p_prime_base_67_at_least_one_violated: z_p_prime_base_67[0]
   + z_p_prime_base_67[1] + z_p_prime_base_67[2] + z_p_prime_base_67[3]
   + z_p_prime_base_67[4] = 1
 infeasibility_p_base_68_at_least_one_violated: z_p_base_68[0]
   + z_p_base_68[1] + z_p_base_68[2] + z_p_base_68[3] + z_p_base_68[4] = 1
 infeasibility_p_prime_base_68_at_least_one_violated: z_p_prime_base_68[0]
   + z_p_prime_base_68[1] + z_p_prime_base_68[2] + z_p_prime_base_68[3]
   + z_p_prime_base_68[4] = 1
 infeasibility_p_base_69_at_least_one_violated: z_p_base_69[0]
   + z_p_base_69[1] + z_p_base_69[2] + z_p_base_69[3] + z_p_base_69[4] = 1
 infeasibility_p_prime_base_69_at_least_one_violated: z_p_prime_base_69[0]
   + z_p_prime_base_69[1] + z_p_prime_base_69[2] + z_p_prime_base_69[3]
   + z_p_prime_base_69[4] = 1
 infeasibility_p_base_70_at_least_one_violated: z_p_base_70[0]
   + z_p_base_70[1] + z_p_base_70[2] + z_p_base_70[3] + z_p_base_70[4] = 1
 infeasibility_p_prime_base_70_at_least_one_violated: z_p_prime_base_70[0]
   + z_p_prime_base_70[1] + z_p_prime_base_70[2] + z_p_prime_base_70[3]
   + z_p_prime_base_70[4] = 1
 infeasibility_p_base_71_at_least_one_violated: z_p_base_71[0]
   + z_p_base_71[1] + z_p_base_71[2] + z_p_base_71[3] + z_p_base_71[4] = 1
 infeasibility_p_prime_base_71_at_least_one_violated: z_p_prime_base_71[0]
   + z_p_prime_base_71[1] + z_p_prime_base_71[2] + z_p_prime_base_71[3]
   + z_p_prime_base_71[4] = 1
 infeasibility_p_base_72_at_least_one_violated: z_p_base_72[0]
   + z_p_base_72[1] + z_p_base_72[2] + z_p_base_72[3] + z_p_base_72[4] = 1
 infeasibility_p_prime_base_72_at_least_one_violated: z_p_prime_base_72[0]
   + z_p_prime_base_72[1] + z_p_prime_base_72[2] + z_p_prime_base_72[3]
   + z_p_prime_base_72[4] = 1
 infeasibility_p_base_73_at_least_one_violated: z_p_base_73[0]
   + z_p_base_73[1] + z_p_base_73[2] + z_p_base_73[3] + z_p_base_73[4] = 1
 infeasibility_p_prime_base_73_at_least_one_violated: z_p_prime_base_73[0]
   + z_p_prime_base_73[1] + z_p_prime_base_73[2] + z_p_prime_base_73[3]
   + z_p_prime_base_73[4] = 1
 infeasibility_p_base_74_at_least_one_violated: z_p_base_74[0]
   + z_p_base_74[1] + z_p_base_74[2] + z_p_base_74[3] + z_p_base_74[4] = 1
 infeasibility_p_prime_base_74_at_least_one_violated: z_p_prime_base_74[0]
   + z_p_prime_base_74[1] + z_p_prime_base_74[2] + z_p_prime_base_74[3]
   + z_p_prime_base_74[4] = 1
 infeasibility_p_base_75_at_least_one_violated: z_p_base_75[0]
   + z_p_base_75[1] + z_p_base_75[2] + z_p_base_75[3] + z_p_base_75[4] = 1
 infeasibility_p_prime_base_75_at_least_one_violated: z_p_prime_base_75[0]
   + z_p_prime_base_75[1] + z_p_prime_base_75[2] + z_p_prime_base_75[3]
   + z_p_prime_base_75[4] = 1
 infeasibility_p_base_76_at_least_one_violated: z_p_base_76[0]
   + z_p_base_76[1] + z_p_base_76[2] + z_p_base_76[3] + z_p_base_76[4] = 1
 infeasibility_p_prime_base_76_at_least_one_violated: z_p_prime_base_76[0]
   + z_p_prime_base_76[1] + z_p_prime_base_76[2] + z_p_prime_base_76[3]
   + z_p_prime_base_76[4] = 1
 infeasibility_p_base_77_at_least_one_violated: z_p_base_77[0]
   + z_p_base_77[1] + z_p_base_77[2] + z_p_base_77[3] + z_p_base_77[4] = 1
 infeasibility_p_prime_base_77_at_least_one_violated: z_p_prime_base_77[0]
   + z_p_prime_base_77[1] + z_p_prime_base_77[2] + z_p_prime_base_77[3]
   + z_p_prime_base_77[4] = 1
 infeasibility_p_base_78_at_least_one_violated: z_p_base_78[0]
   + z_p_base_78[1] + z_p_base_78[2] + z_p_base_78[3] + z_p_base_78[4] = 1
 infeasibility_p_prime_base_78_at_least_one_violated: z_p_prime_base_78[0]
   + z_p_prime_base_78[1] + z_p_prime_base_78[2] + z_p_prime_base_78[3]
   + z_p_prime_base_78[4] = 1
 infeasibility_p_base_79_at_least_one_violated: z_p_base_79[0]
   + z_p_base_79[1] + z_p_base_79[2] + z_p_base_79[3] + z_p_base_79[4] = 1
 infeasibility_p_prime_base_79_at_least_one_violated: z_p_prime_base_79[0]
   + z_p_prime_base_79[1] + z_p_prime_base_79[2] + z_p_prime_base_79[3]
   + z_p_prime_base_79[4] = 1
 infeasibility_p_base_80_at_least_one_violated: z_p_base_80[0]
   + z_p_base_80[1] + z_p_base_80[2] + z_p_base_80[3] + z_p_base_80[4] = 1
 infeasibility_p_prime_base_80_at_least_one_violated: z_p_prime_base_80[0]
   + z_p_prime_base_80[1] + z_p_prime_base_80[2] + z_p_prime_base_80[3]
   + z_p_prime_base_80[4] = 1
 infeasibility_p_base_81_at_least_one_violated: z_p_base_81[0]
   + z_p_base_81[1] + z_p_base_81[2] + z_p_base_81[3] + z_p_base_81[4] = 1
 infeasibility_p_prime_base_81_at_least_one_violated: z_p_prime_base_81[0]
   + z_p_prime_base_81[1] + z_p_prime_base_81[2] + z_p_prime_base_81[3]
   + z_p_prime_base_81[4] = 1
 infeasibility_p_base_82_at_least_one_violated: z_p_base_82[0]
   + z_p_base_82[1] + z_p_base_82[2] + z_p_base_82[3] + z_p_base_82[4] = 1
 infeasibility_p_prime_base_82_at_least_one_violated: z_p_prime_base_82[0]
   + z_p_prime_base_82[1] + z_p_prime_base_82[2] + z_p_prime_base_82[3]
   + z_p_prime_base_82[4] = 1
 infeasibility_p_base_83_at_least_one_violated: z_p_base_83[0]
   + z_p_base_83[1] + z_p_base_83[2] + z_p_base_83[3] + z_p_base_83[4] = 1
 infeasibility_p_prime_base_83_at_least_one_violated: z_p_prime_base_83[0]
   + z_p_prime_base_83[1] + z_p_prime_base_83[2] + z_p_prime_base_83[3]
   + z_p_prime_base_83[4] = 1
 infeasibility_p_base_84_at_least_one_violated: z_p_base_84[0]
   + z_p_base_84[1] + z_p_base_84[2] + z_p_base_84[3] + z_p_base_84[4] = 1
 infeasibility_p_prime_base_84_at_least_one_violated: z_p_prime_base_84[0]
   + z_p_prime_base_84[1] + z_p_prime_base_84[2] + z_p_prime_base_84[3]
   + z_p_prime_base_84[4] = 1
 infeasibility_p_base_85_at_least_one_violated: z_p_base_85[0]
   + z_p_base_85[1] + z_p_base_85[2] + z_p_base_85[3] + z_p_base_85[4] = 1
 infeasibility_p_prime_base_85_at_least_one_violated: z_p_prime_base_85[0]
   + z_p_prime_base_85[1] + z_p_prime_base_85[2] + z_p_prime_base_85[3]
   + z_p_prime_base_85[4] = 1
 infeasibility_p_base_86_at_least_one_violated: z_p_base_86[0]
   + z_p_base_86[1] + z_p_base_86[2] + z_p_base_86[3] + z_p_base_86[4] = 1
 infeasibility_p_prime_base_86_at_least_one_violated: z_p_prime_base_86[0]
   + z_p_prime_base_86[1] + z_p_prime_base_86[2] + z_p_prime_base_86[3]
   + z_p_prime_base_86[4] = 1
 infeasibility_p_base_87_at_least_one_violated: z_p_base_87[0]
   + z_p_base_87[1] + z_p_base_87[2] + z_p_base_87[3] + z_p_base_87[4] = 1
 infeasibility_p_prime_base_87_at_least_one_violated: z_p_prime_base_87[0]
   + z_p_prime_base_87[1] + z_p_prime_base_87[2] + z_p_prime_base_87[3]
   + z_p_prime_base_87[4] = 1
 infeasibility_p_base_88_at_least_one_violated: z_p_base_88[0]
   + z_p_base_88[1] + z_p_base_88[2] + z_p_base_88[3] + z_p_base_88[4] = 1
 infeasibility_p_prime_base_88_at_least_one_violated: z_p_prime_base_88[0]
   + z_p_prime_base_88[1] + z_p_prime_base_88[2] + z_p_prime_base_88[3]
   + z_p_prime_base_88[4] = 1
 infeasibility_p_base_89_at_least_one_violated: z_p_base_89[0]
   + z_p_base_89[1] + z_p_base_89[2] + z_p_base_89[3] + z_p_base_89[4] = 1
 infeasibility_p_prime_base_89_at_least_one_violated: z_p_prime_base_89[0]
   + z_p_prime_base_89[1] + z_p_prime_base_89[2] + z_p_prime_base_89[3]
   + z_p_prime_base_89[4] = 1
 infeasibility_p_base_90_at_least_one_violated: z_p_base_90[0]
   + z_p_base_90[1] + z_p_base_90[2] + z_p_base_90[3] + z_p_base_90[4] = 1
 infeasibility_p_prime_base_90_at_least_one_violated: z_p_prime_base_90[0]
   + z_p_prime_base_90[1] + z_p_prime_base_90[2] + z_p_prime_base_90[3]
   + z_p_prime_base_90[4] = 1
 infeasibility_p_base_91_at_least_one_violated: z_p_base_91[0]
   + z_p_base_91[1] + z_p_base_91[2] + z_p_base_91[3] + z_p_base_91[4] = 1
 infeasibility_p_prime_base_91_at_least_one_violated: z_p_prime_base_91[0]
   + z_p_prime_base_91[1] + z_p_prime_base_91[2] + z_p_prime_base_91[3]
   + z_p_prime_base_91[4] = 1
 infeasibility_p_base_92_at_least_one_violated: z_p_base_92[0]
   + z_p_base_92[1] + z_p_base_92[2] + z_p_base_92[3] + z_p_base_92[4] = 1
 infeasibility_p_prime_base_92_at_least_one_violated: z_p_prime_base_92[0]
   + z_p_prime_base_92[1] + z_p_prime_base_92[2] + z_p_prime_base_92[3]
   + z_p_prime_base_92[4] = 1
 infeasibility_p_base_93_at_least_one_violated: z_p_base_93[0]
   + z_p_base_93[1] + z_p_base_93[2] + z_p_base_93[3] + z_p_base_93[4] = 1
 infeasibility_p_prime_base_93_at_least_one_violated: z_p_prime_base_93[0]
   + z_p_prime_base_93[1] + z_p_prime_base_93[2] + z_p_prime_base_93[3]
   + z_p_prime_base_93[4] = 1
 infeasibility_p_base_94_at_least_one_violated: z_p_base_94[0]
   + z_p_base_94[1] + z_p_base_94[2] + z_p_base_94[3] + z_p_base_94[4] = 1
 infeasibility_p_prime_base_94_at_least_one_violated: z_p_prime_base_94[0]
   + z_p_prime_base_94[1] + z_p_prime_base_94[2] + z_p_prime_base_94[3]
   + z_p_prime_base_94[4] = 1
 infeasibility_p_base_95_at_least_one_violated: z_p_base_95[0]
   + z_p_base_95[1] + z_p_base_95[2] + z_p_base_95[3] + z_p_base_95[4] = 1
 infeasibility_p_prime_base_95_at_least_one_violated: z_p_prime_base_95[0]
   + z_p_prime_base_95[1] + z_p_prime_base_95[2] + z_p_prime_base_95[3]
   + z_p_prime_base_95[4] = 1
 infeasibility_p_base_96_at_least_one_violated: z_p_base_96[0]
   + z_p_base_96[1] + z_p_base_96[2] + z_p_base_96[3] + z_p_base_96[4] = 1
 infeasibility_p_prime_base_96_at_least_one_violated: z_p_prime_base_96[0]
   + z_p_prime_base_96[1] + z_p_prime_base_96[2] + z_p_prime_base_96[3]
   + z_p_prime_base_96[4] = 1
 infeasibility_p_base_97_at_least_one_violated: z_p_base_97[0]
   + z_p_base_97[1] + z_p_base_97[2] + z_p_base_97[3] + z_p_base_97[4] = 1
 infeasibility_p_prime_base_97_at_least_one_violated: z_p_prime_base_97[0]
   + z_p_prime_base_97[1] + z_p_prime_base_97[2] + z_p_prime_base_97[3]
   + z_p_prime_base_97[4] = 1
 infeasibility_p_base_98_at_least_one_violated: z_p_base_98[0]
   + z_p_base_98[1] + z_p_base_98[2] + z_p_base_98[3] + z_p_base_98[4] = 1
 infeasibility_p_prime_base_98_at_least_one_violated: z_p_prime_base_98[0]
   + z_p_prime_base_98[1] + z_p_prime_base_98[2] + z_p_prime_base_98[3]
   + z_p_prime_base_98[4] = 1
 infeasibility_p_base_99_at_least_one_violated: z_p_base_99[0]
   + z_p_base_99[1] + z_p_base_99[2] + z_p_base_99[3] + z_p_base_99[4] = 1
 infeasibility_p_prime_base_99_at_least_one_violated: z_p_prime_base_99[0]
   + z_p_prime_base_99[1] + z_p_prime_base_99[2] + z_p_prime_base_99[3]
   + z_p_prime_base_99[4] = 1
 infeasibility_p_base_100_at_least_one_violated: z_p_base_100[0]
   + z_p_base_100[1] + z_p_base_100[2] + z_p_base_100[3] + z_p_base_100[4]
   = 1
 infeasibility_p_prime_base_100_at_least_one_violated:
   z_p_prime_base_100[0] + z_p_prime_base_100[1] + z_p_prime_base_100[2]
   + z_p_prime_base_100[3] + z_p_prime_base_100[4] = 1
 infeasibility_p_base_101_at_least_one_violated: z_p_base_101[0]
   + z_p_base_101[1] + z_p_base_101[2] + z_p_base_101[3] + z_p_base_101[4]
   = 1
 infeasibility_p_prime_base_101_at_least_one_violated:
   z_p_prime_base_101[0] + z_p_prime_base_101[1] + z_p_prime_base_101[2]
   + z_p_prime_base_101[3] + z_p_prime_base_101[4] = 1
 infeasibility_p_base_102_at_least_one_violated: z_p_base_102[0]
   + z_p_base_102[1] + z_p_base_102[2] + z_p_base_102[3] + z_p_base_102[4]
   = 1
 infeasibility_p_prime_base_102_at_least_one_violated:
   z_p_prime_base_102[0] + z_p_prime_base_102[1] + z_p_prime_base_102[2]
   + z_p_prime_base_102[3] + z_p_prime_base_102[4] = 1
 infeasibility_p_base_103_at_least_one_violated: z_p_base_103[0]
   + z_p_base_103[1] + z_p_base_103[2] + z_p_base_103[3] + z_p_base_103[4]
   = 1
 infeasibility_p_prime_base_103_at_least_one_violated:
   z_p_prime_base_103[0] + z_p_prime_base_103[1] + z_p_prime_base_103[2]
   + z_p_prime_base_103[3] + z_p_prime_base_103[4] = 1
 infeasibility_p_base_104_at_least_one_violated: z_p_base_104[0]
   + z_p_base_104[1] + z_p_base_104[2] + z_p_base_104[3] + z_p_base_104[4]
   = 1
 infeasibility_p_prime_base_104_at_least_one_violated:
   z_p_prime_base_104[0] + z_p_prime_base_104[1] + z_p_prime_base_104[2]
   + z_p_prime_base_104[3] + z_p_prime_base_104[4] = 1
 infeasibility_p_base_105_at_least_one_violated: z_p_base_105[0]
   + z_p_base_105[1] + z_p_base_105[2] + z_p_base_105[3] + z_p_base_105[4]
   = 1
 infeasibility_p_prime_base_105_at_least_one_violated:
   z_p_prime_base_105[0] + z_p_prime_base_105[1] + z_p_prime_base_105[2]
   + z_p_prime_base_105[3] + z_p_prime_base_105[4] = 1
 infeasibility_p_base_106_at_least_one_violated: z_p_base_106[0]
   + z_p_base_106[1] + z_p_base_106[2] + z_p_base_106[3] + z_p_base_106[4]
   = 1
 infeasibility_p_prime_base_106_at_least_one_violated:
   z_p_prime_base_106[0] + z_p_prime_base_106[1] + z_p_prime_base_106[2]
   + z_p_prime_base_106[3] + z_p_prime_base_106[4] = 1
 infeasibility_p_base_107_at_least_one_violated: z_p_base_107[0]
   + z_p_base_107[1] + z_p_base_107[2] + z_p_base_107[3] + z_p_base_107[4]
   = 1
 infeasibility_p_prime_base_107_at_least_one_violated:
   z_p_prime_base_107[0] + z_p_prime_base_107[1] + z_p_prime_base_107[2]
   + z_p_prime_base_107[3] + z_p_prime_base_107[4] = 1
 infeasibility_p_base_108_at_least_one_violated: z_p_base_108[0]
   + z_p_base_108[1] + z_p_base_108[2] + z_p_base_108[3] + z_p_base_108[4]
   = 1
 infeasibility_p_prime_base_108_at_least_one_violated:
   z_p_prime_base_108[0] + z_p_prime_base_108[1] + z_p_prime_base_108[2]
   + z_p_prime_base_108[3] + z_p_prime_base_108[4] = 1
 infeasibility_p_base_109_at_least_one_violated: z_p_base_109[0]
   + z_p_base_109[1] + z_p_base_109[2] + z_p_base_109[3] + z_p_base_109[4]
   = 1
 infeasibility_p_prime_base_109_at_least_one_violated:
   z_p_prime_base_109[0] + z_p_prime_base_109[1] + z_p_prime_base_109[2]
   + z_p_prime_base_109[3] + z_p_prime_base_109[4] = 1
 infeasibility_p_base_110_at_least_one_violated: z_p_base_110[0]
   + z_p_base_110[1] + z_p_base_110[2] + z_p_base_110[3] + z_p_base_110[4]
   = 1
 infeasibility_p_prime_base_110_at_least_one_violated:
   z_p_prime_base_110[0] + z_p_prime_base_110[1] + z_p_prime_base_110[2]
   + z_p_prime_base_110[3] + z_p_prime_base_110[4] = 1
 infeasibility_p_base_111_at_least_one_violated: z_p_base_111[0]
   + z_p_base_111[1] + z_p_base_111[2] + z_p_base_111[3] + z_p_base_111[4]
   = 1
 infeasibility_p_prime_base_111_at_least_one_violated:
   z_p_prime_base_111[0] + z_p_prime_base_111[1] + z_p_prime_base_111[2]
   + z_p_prime_base_111[3] + z_p_prime_base_111[4] = 1
 infeasibility_p_base_112_at_least_one_violated: z_p_base_112[0]
   + z_p_base_112[1] + z_p_base_112[2] + z_p_base_112[3] + z_p_base_112[4]
   = 1
 infeasibility_p_prime_base_112_at_least_one_violated:
   z_p_prime_base_112[0] + z_p_prime_base_112[1] + z_p_prime_base_112[2]
   + z_p_prime_base_112[3] + z_p_prime_base_112[4] = 1
 infeasibility_p_base_113_at_least_one_violated: z_p_base_113[0]
   + z_p_base_113[1] + z_p_base_113[2] + z_p_base_113[3] + z_p_base_113[4]
   = 1
 infeasibility_p_prime_base_113_at_least_one_violated:
   z_p_prime_base_113[0] + z_p_prime_base_113[1] + z_p_prime_base_113[2]
   + z_p_prime_base_113[3] + z_p_prime_base_113[4] = 1
 infeasibility_p_base_114_at_least_one_violated: z_p_base_114[0]
   + z_p_base_114[1] + z_p_base_114[2] + z_p_base_114[3] + z_p_base_114[4]
   = 1
 infeasibility_p_prime_base_114_at_least_one_violated:
   z_p_prime_base_114[0] + z_p_prime_base_114[1] + z_p_prime_base_114[2]
   + z_p_prime_base_114[3] + z_p_prime_base_114[4] = 1
 infeasibility_p_base_115_at_least_one_violated: z_p_base_115[0]
   + z_p_base_115[1] + z_p_base_115[2] + z_p_base_115[3] + z_p_base_115[4]
   = 1
 infeasibility_p_prime_base_115_at_least_one_violated:
   z_p_prime_base_115[0] + z_p_prime_base_115[1] + z_p_prime_base_115[2]
   + z_p_prime_base_115[3] + z_p_prime_base_115[4] = 1
 infeasibility_p_base_116_at_least_one_violated: z_p_base_116[0]
   + z_p_base_116[1] + z_p_base_116[2] + z_p_base_116[3] + z_p_base_116[4]
   = 1
 infeasibility_p_prime_base_116_at_least_one_violated:
   z_p_prime_base_116[0] + z_p_prime_base_116[1] + z_p_prime_base_116[2]
   + z_p_prime_base_116[3] + z_p_prime_base_116[4] = 1
 infeasibility_p_base_117_at_least_one_violated: z_p_base_117[0]
   + z_p_base_117[1] + z_p_base_117[2] + z_p_base_117[3] + z_p_base_117[4]
   = 1
 infeasibility_p_prime_base_117_at_least_one_violated:
   z_p_prime_base_117[0] + z_p_prime_base_117[1] + z_p_prime_base_117[2]
   + z_p_prime_base_117[3] + z_p_prime_base_117[4] = 1
 infeasibility_p_base_118_at_least_one_violated: z_p_base_118[0]
   + z_p_base_118[1] + z_p_base_118[2] + z_p_base_118[3] + z_p_base_118[4]
   = 1
 infeasibility_p_prime_base_118_at_least_one_violated:
   z_p_prime_base_118[0] + z_p_prime_base_118[1] + z_p_prime_base_118[2]
   + z_p_prime_base_118[3] + z_p_prime_base_118[4] = 1
 infeasibility_p_base_119_at_least_one_violated: z_p_base_119[0]
   + z_p_base_119[1] + z_p_base_119[2] + z_p_base_119[3] + z_p_base_119[4]
   = 1
 infeasibility_p_prime_base_119_at_least_one_violated:
   z_p_prime_base_119[0] + z_p_prime_base_119[1] + z_p_prime_base_119[2]
   + z_p_prime_base_119[3] + z_p_prime_base_119[4] = 1
 infeasibility_p_base_120_at_least_one_violated: z_p_base_120[0]
   + z_p_base_120[1] + z_p_base_120[2] + z_p_base_120[3] + z_p_base_120[4]
   = 1
 infeasibility_p_prime_base_120_at_least_one_violated:
   z_p_prime_base_120[0] + z_p_prime_base_120[1] + z_p_prime_base_120[2]
   + z_p_prime_base_120[3] + z_p_prime_base_120[4] = 1
 infeasibility_p_base_121_at_least_one_violated: z_p_base_121[0]
   + z_p_base_121[1] + z_p_base_121[2] + z_p_base_121[3] + z_p_base_121[4]
   = 1
 infeasibility_p_prime_base_121_at_least_one_violated:
   z_p_prime_base_121[0] + z_p_prime_base_121[1] + z_p_prime_base_121[2]
   + z_p_prime_base_121[3] + z_p_prime_base_121[4] = 1
 infeasibility_p_base_122_at_least_one_violated: z_p_base_122[0]
   + z_p_base_122[1] + z_p_base_122[2] + z_p_base_122[3] + z_p_base_122[4]
   = 1
 infeasibility_p_prime_base_122_at_least_one_violated:
   z_p_prime_base_122[0] + z_p_prime_base_122[1] + z_p_prime_base_122[2]
   + z_p_prime_base_122[3] + z_p_prime_base_122[4] = 1
 infeasibility_p_base_123_at_least_one_violated: z_p_base_123[0]
   + z_p_base_123[1] + z_p_base_123[2] + z_p_base_123[3] + z_p_base_123[4]
   = 1
 infeasibility_p_prime_base_123_at_least_one_violated:
   z_p_prime_base_123[0] + z_p_prime_base_123[1] + z_p_prime_base_123[2]
   + z_p_prime_base_123[3] + z_p_prime_base_123[4] = 1
 infeasibility_p_base_124_at_least_one_violated: z_p_base_124[0]
   + z_p_base_124[1] + z_p_base_124[2] + z_p_base_124[3] + z_p_base_124[4]
   = 1
 infeasibility_p_prime_base_124_at_least_one_violated:
   z_p_prime_base_124[0] + z_p_prime_base_124[1] + z_p_prime_base_124[2]
   + z_p_prime_base_124[3] + z_p_prime_base_124[4] = 1
 infeasibility_p_base_125_at_least_one_violated: z_p_base_125[0]
   + z_p_base_125[1] + z_p_base_125[2] + z_p_base_125[3] + z_p_base_125[4]
   = 1
 infeasibility_p_prime_base_125_at_least_one_violated:
   z_p_prime_base_125[0] + z_p_prime_base_125[1] + z_p_prime_base_125[2]
   + z_p_prime_base_125[3] + z_p_prime_base_125[4] = 1
 infeasibility_p_base_126_at_least_one_violated: z_p_base_126[0]
   + z_p_base_126[1] + z_p_base_126[2] + z_p_base_126[3] + z_p_base_126[4]
   = 1
 infeasibility_p_prime_base_126_at_least_one_violated:
   z_p_prime_base_126[0] + z_p_prime_base_126[1] + z_p_prime_base_126[2]
   + z_p_prime_base_126[3] + z_p_prime_base_126[4] = 1
 infeasibility_p_base_127_at_least_one_violated: z_p_base_127[0]
   + z_p_base_127[1] + z_p_base_127[2] + z_p_base_127[3] + z_p_base_127[4]
   = 1
 infeasibility_p_prime_base_127_at_least_one_violated:
   z_p_prime_base_127[0] + z_p_prime_base_127[1] + z_p_prime_base_127[2]
   + z_p_prime_base_127[3] + z_p_prime_base_127[4] = 1
 infeasibility_p_base_128_at_least_one_violated: z_p_base_128[0]
   + z_p_base_128[1] + z_p_base_128[2] + z_p_base_128[3] + z_p_base_128[4]
   = 1
 infeasibility_p_prime_base_128_at_least_one_violated:
   z_p_prime_base_128[0] + z_p_prime_base_128[1] + z_p_prime_base_128[2]
   + z_p_prime_base_128[3] + z_p_prime_base_128[4] = 1
 infeasibility_p_base_129_at_least_one_violated: z_p_base_129[0]
   + z_p_base_129[1] + z_p_base_129[2] + z_p_base_129[3] + z_p_base_129[4]
   = 1
 infeasibility_p_prime_base_129_at_least_one_violated:
   z_p_prime_base_129[0] + z_p_prime_base_129[1] + z_p_prime_base_129[2]
   + z_p_prime_base_129[3] + z_p_prime_base_129[4] = 1
 infeasibility_p_base_130_at_least_one_violated: z_p_base_130[0]
   + z_p_base_130[1] + z_p_base_130[2] + z_p_base_130[3] + z_p_base_130[4]
   = 1
 infeasibility_p_prime_base_130_at_least_one_violated:
   z_p_prime_base_130[0] + z_p_prime_base_130[1] + z_p_prime_base_130[2]
   + z_p_prime_base_130[3] + z_p_prime_base_130[4] = 1
 infeasibility_p_base_131_at_least_one_violated: z_p_base_131[0]
   + z_p_base_131[1] + z_p_base_131[2] + z_p_base_131[3] + z_p_base_131[4]
   = 1
 infeasibility_p_prime_base_131_at_least_one_violated:
   z_p_prime_base_131[0] + z_p_prime_base_131[1] + z_p_prime_base_131[2]
   + z_p_prime_base_131[3] + z_p_prime_base_131[4] = 1
 infeasibility_p_base_132_at_least_one_violated: z_p_base_132[0]
   + z_p_base_132[1] + z_p_base_132[2] + z_p_base_132[3] + z_p_base_132[4]
   = 1
 infeasibility_p_prime_base_132_at_least_one_violated:
   z_p_prime_base_132[0] + z_p_prime_base_132[1] + z_p_prime_base_132[2]
   + z_p_prime_base_132[3] + z_p_prime_base_132[4] = 1
 infeasibility_p_base_133_at_least_one_violated: z_p_base_133[0]
   + z_p_base_133[1] + z_p_base_133[2] + z_p_base_133[3] + z_p_base_133[4]
   = 1
 infeasibility_p_prime_base_133_at_least_one_violated:
   z_p_prime_base_133[0] + z_p_prime_base_133[1] + z_p_prime_base_133[2]
   + z_p_prime_base_133[3] + z_p_prime_base_133[4] = 1
 infeasibility_p_base_134_at_least_one_violated: z_p_base_134[0]
   + z_p_base_134[1] + z_p_base_134[2] + z_p_base_134[3] + z_p_base_134[4]
   = 1
 infeasibility_p_prime_base_134_at_least_one_violated:
   z_p_prime_base_134[0] + z_p_prime_base_134[1] + z_p_prime_base_134[2]
   + z_p_prime_base_134[3] + z_p_prime_base_134[4] = 1
 infeasibility_p_base_135_at_least_one_violated: z_p_base_135[0]
   + z_p_base_135[1] + z_p_base_135[2] + z_p_base_135[3] + z_p_base_135[4]
   = 1
 infeasibility_p_prime_base_135_at_least_one_violated:
   z_p_prime_base_135[0] + z_p_prime_base_135[1] + z_p_prime_base_135[2]
   + z_p_prime_base_135[3] + z_p_prime_base_135[4] = 1
 infeasibility_p_base_136_at_least_one_violated: z_p_base_136[0]
   + z_p_base_136[1] + z_p_base_136[2] + z_p_base_136[3] + z_p_base_136[4]
   = 1
 infeasibility_p_prime_base_136_at_least_one_violated:
   z_p_prime_base_136[0] + z_p_prime_base_136[1] + z_p_prime_base_136[2]
   + z_p_prime_base_136[3] + z_p_prime_base_136[4] = 1
 infeasibility_p_base_137_at_least_one_violated: z_p_base_137[0]
   + z_p_base_137[1] + z_p_base_137[2] + z_p_base_137[3] + z_p_base_137[4]
   = 1
 infeasibility_p_prime_base_137_at_least_one_violated:
   z_p_prime_base_137[0] + z_p_prime_base_137[1] + z_p_prime_base_137[2]
   + z_p_prime_base_137[3] + z_p_prime_base_137[4] = 1
 infeasibility_p_base_138_at_least_one_violated: z_p_base_138[0]
   + z_p_base_138[1] + z_p_base_138[2] + z_p_base_138[3] + z_p_base_138[4]
   = 1
 infeasibility_p_prime_base_138_at_least_one_violated:
   z_p_prime_base_138[0] + z_p_prime_base_138[1] + z_p_prime_base_138[2]
   + z_p_prime_base_138[3] + z_p_prime_base_138[4] = 1
 infeasibility_p_base_139_at_least_one_violated: z_p_base_139[0]
   + z_p_base_139[1] + z_p_base_139[2] + z_p_base_139[3] + z_p_base_139[4]
   = 1
 infeasibility_p_prime_base_139_at_least_one_violated:
   z_p_prime_base_139[0] + z_p_prime_base_139[1] + z_p_prime_base_139[2]
   + z_p_prime_base_139[3] + z_p_prime_base_139[4] = 1
 infeasibility_p_base_140_at_least_one_violated: z_p_base_140[0]
   + z_p_base_140[1] + z_p_base_140[2] + z_p_base_140[3] + z_p_base_140[4]
   = 1
 infeasibility_p_prime_base_140_at_least_one_violated:
   z_p_prime_base_140[0] + z_p_prime_base_140[1] + z_p_prime_base_140[2]
   + z_p_prime_base_140[3] + z_p_prime_base_140[4] = 1
 infeasibility_p_base_141_at_least_one_violated: z_p_base_141[0]
   + z_p_base_141[1] + z_p_base_141[2] + z_p_base_141[3] + z_p_base_141[4]
   = 1
 infeasibility_p_prime_base_141_at_least_one_violated:
   z_p_prime_base_141[0] + z_p_prime_base_141[1] + z_p_prime_base_141[2]
   + z_p_prime_base_141[3] + z_p_prime_base_141[4] = 1
 infeasibility_p_base_142_at_least_one_violated: z_p_base_142[0]
   + z_p_base_142[1] + z_p_base_142[2] + z_p_base_142[3] + z_p_base_142[4]
   = 1
 infeasibility_p_prime_base_142_at_least_one_violated:
   z_p_prime_base_142[0] + z_p_prime_base_142[1] + z_p_prime_base_142[2]
   + z_p_prime_base_142[3] + z_p_prime_base_142[4] = 1
 infeasibility_p_base_143_at_least_one_violated: z_p_base_143[0]
   + z_p_base_143[1] + z_p_base_143[2] + z_p_base_143[3] + z_p_base_143[4]
   = 1
 infeasibility_p_prime_base_143_at_least_one_violated:
   z_p_prime_base_143[0] + z_p_prime_base_143[1] + z_p_prime_base_143[2]
   + z_p_prime_base_143[3] + z_p_prime_base_143[4] = 1
 infeasibility_p_base_144_at_least_one_violated: z_p_base_144[0]
   + z_p_base_144[1] + z_p_base_144[2] + z_p_base_144[3] + z_p_base_144[4]
   = 1
 infeasibility_p_prime_base_144_at_least_one_violated:
   z_p_prime_base_144[0] + z_p_prime_base_144[1] + z_p_prime_base_144[2]
   + z_p_prime_base_144[3] + z_p_prime_base_144[4] = 1
 infeasibility_p_base_145_at_least_one_violated: z_p_base_145[0]
   + z_p_base_145[1] + z_p_base_145[2] + z_p_base_145[3] + z_p_base_145[4]
   = 1
 infeasibility_p_prime_base_145_at_least_one_violated:
   z_p_prime_base_145[0] + z_p_prime_base_145[1] + z_p_prime_base_145[2]
   + z_p_prime_base_145[3] + z_p_prime_base_145[4] = 1
 infeasibility_p_base_146_at_least_one_violated: z_p_base_146[0]
   + z_p_base_146[1] + z_p_base_146[2] + z_p_base_146[3] + z_p_base_146[4]
   = 1
 infeasibility_p_prime_base_146_at_least_one_violated:
   z_p_prime_base_146[0] + z_p_prime_base_146[1] + z_p_prime_base_146[2]
   + z_p_prime_base_146[3] + z_p_prime_base_146[4] = 1
 infeasibility_p_base_147_at_least_one_violated: z_p_base_147[0]
   + z_p_base_147[1] + z_p_base_147[2] + z_p_base_147[3] + z_p_base_147[4]
   = 1
 infeasibility_p_prime_base_147_at_least_one_violated:
   z_p_prime_base_147[0] + z_p_prime_base_147[1] + z_p_prime_base_147[2]
   + z_p_prime_base_147[3] + z_p_prime_base_147[4] = 1
 infeasibility_p_base_148_at_least_one_violated: z_p_base_148[0]
   + z_p_base_148[1] + z_p_base_148[2] + z_p_base_148[3] + z_p_base_148[4]
   = 1
 infeasibility_p_prime_base_148_at_least_one_violated:
   z_p_prime_base_148[0] + z_p_prime_base_148[1] + z_p_prime_base_148[2]
   + z_p_prime_base_148[3] + z_p_prime_base_148[4] = 1
 infeasibility_p_base_149_at_least_one_violated: z_p_base_149[0]
   + z_p_base_149[1] + z_p_base_149[2] + z_p_base_149[3] + z_p_base_149[4]
   = 1
 infeasibility_p_prime_base_149_at_least_one_violated:
   z_p_prime_base_149[0] + z_p_prime_base_149[1] + z_p_prime_base_149[2]
   + z_p_prime_base_149[3] + z_p_prime_base_149[4] = 1
 infeasibility_p_base_150_at_least_one_violated: z_p_base_150[0]
   + z_p_base_150[1] + z_p_base_150[2] + z_p_base_150[3] + z_p_base_150[4]
   = 1
 infeasibility_p_prime_base_150_at_least_one_violated:
   z_p_prime_base_150[0] + z_p_prime_base_150[1] + z_p_prime_base_150[2]
   + z_p_prime_base_150[3] + z_p_prime_base_150[4] = 1
 infeasibility_p_base_151_at_least_one_violated: z_p_base_151[0]
   + z_p_base_151[1] + z_p_base_151[2] + z_p_base_151[3] + z_p_base_151[4]
   = 1
 infeasibility_p_prime_base_151_at_least_one_violated:
   z_p_prime_base_151[0] + z_p_prime_base_151[1] + z_p_prime_base_151[2]
   + z_p_prime_base_151[3] + z_p_prime_base_151[4] = 1
 infeasibility_p_base_152_at_least_one_violated: z_p_base_152[0]
   + z_p_base_152[1] + z_p_base_152[2] + z_p_base_152[3] + z_p_base_152[4]
   = 1
 infeasibility_p_prime_base_152_at_least_one_violated:
   z_p_prime_base_152[0] + z_p_prime_base_152[1] + z_p_prime_base_152[2]
   + z_p_prime_base_152[3] + z_p_prime_base_152[4] = 1
 infeasibility_p_base_153_at_least_one_violated: z_p_base_153[0]
   + z_p_base_153[1] + z_p_base_153[2] + z_p_base_153[3] + z_p_base_153[4]
   = 1
 infeasibility_p_prime_base_153_at_least_one_violated:
   z_p_prime_base_153[0] + z_p_prime_base_153[1] + z_p_prime_base_153[2]
   + z_p_prime_base_153[3] + z_p_prime_base_153[4] = 1
 infeasibility_p_base_154_at_least_one_violated: z_p_base_154[0]
   + z_p_base_154[1] + z_p_base_154[2] + z_p_base_154[3] + z_p_base_154[4]
   = 1
 infeasibility_p_prime_base_154_at_least_one_violated:
   z_p_prime_base_154[0] + z_p_prime_base_154[1] + z_p_prime_base_154[2]
   + z_p_prime_base_154[3] + z_p_prime_base_154[4] = 1
 infeasibility_p_base_155_at_least_one_violated: z_p_base_155[0]
   + z_p_base_155[1] + z_p_base_155[2] + z_p_base_155[3] + z_p_base_155[4]
   = 1
 infeasibility_p_prime_base_155_at_least_one_violated:
   z_p_prime_base_155[0] + z_p_prime_base_155[1] + z_p_prime_base_155[2]
   + z_p_prime_base_155[3] + z_p_prime_base_155[4] = 1
 infeasibility_p_base_156_at_least_one_violated: z_p_base_156[0]
   + z_p_base_156[1] + z_p_base_156[2] + z_p_base_156[3] + z_p_base_156[4]
   = 1
 infeasibility_p_prime_base_156_at_least_one_violated:
   z_p_prime_base_156[0] + z_p_prime_base_156[1] + z_p_prime_base_156[2]
   + z_p_prime_base_156[3] + z_p_prime_base_156[4] = 1
 infeasibility_p_base_157_at_least_one_violated: z_p_base_157[0]
   + z_p_base_157[1] + z_p_base_157[2] + z_p_base_157[3] + z_p_base_157[4]
   = 1
 infeasibility_p_prime_base_157_at_least_one_violated:
   z_p_prime_base_157[0] + z_p_prime_base_157[1] + z_p_prime_base_157[2]
   + z_p_prime_base_157[3] + z_p_prime_base_157[4] = 1
 infeasibility_p_base_158_at_least_one_violated: z_p_base_158[0]
   + z_p_base_158[1] + z_p_base_158[2] + z_p_base_158[3] + z_p_base_158[4]
   = 1
 infeasibility_p_prime_base_158_at_least_one_violated:
   z_p_prime_base_158[0] + z_p_prime_base_158[1] + z_p_prime_base_158[2]
   + z_p_prime_base_158[3] + z_p_prime_base_158[4] = 1
 infeasibility_p_base_159_at_least_one_violated: z_p_base_159[0]
   + z_p_base_159[1] + z_p_base_159[2] + z_p_base_159[3] + z_p_base_159[4]
   = 1
 infeasibility_p_prime_base_159_at_least_one_violated:
   z_p_prime_base_159[0] + z_p_prime_base_159[1] + z_p_prime_base_159[2]
   + z_p_prime_base_159[3] + z_p_prime_base_159[4] = 1
 infeasibility_p_base_160_at_least_one_violated: z_p_base_160[0]
   + z_p_base_160[1] + z_p_base_160[2] + z_p_base_160[3] + z_p_base_160[4]
   = 1
 infeasibility_p_prime_base_160_at_least_one_violated:
   z_p_prime_base_160[0] + z_p_prime_base_160[1] + z_p_prime_base_160[2]
   + z_p_prime_base_160[3] + z_p_prime_base_160[4] = 1
 infeasibility_p_base_161_at_least_one_violated: z_p_base_161[0]
   + z_p_base_161[1] + z_p_base_161[2] + z_p_base_161[3] + z_p_base_161[4]
   = 1
 infeasibility_p_prime_base_161_at_least_one_violated:
   z_p_prime_base_161[0] + z_p_prime_base_161[1] + z_p_prime_base_161[2]
   + z_p_prime_base_161[3] + z_p_prime_base_161[4] = 1
 infeasibility_p_base_B_prime_at_least_one_violated: z_p_base_B_prime[0]
   + z_p_base_B_prime[1] + z_p_base_B_prime[2] + z_p_base_B_prime[3]
   + z_p_base_B_prime[4] >= 1
 infeasibility_p_prime_base_B_at_least_one_violated: z_p_prime_base_B[0]
   + z_p_prime_base_B[1] + z_p_prime_base_B[2] + z_p_prime_base_B[3]
   + z_p_prime_base_B[4] >= 1
 selection_monotonicity_violation: 0.5 alpha[0] + 0.5 alpha[1]
   + 0.5 alpha[2] - 0.5 alpha[3] + 0.5 alpha[4] + 0.5 p[0] - 0.5 p[1]
   - 0.5 p[2] + 0.5 p[3] - 0.5 p[4] <= -5e-06
 infeasibility_p_base_1_row_0_leq_0: 5e-06 z_p_base_1[0] + [
   0.5 p[0] * z_p_base_1[0] + 0.5 p[1] * z_p_base_1[0]
   - 0.5 p[2] * z_p_base_1[0] - 0.5 p[3] * z_p_base_1[0]
   - 0.5 p[4] * z_p_base_1[0] ] <= 0
 infeasibility_p_base_1_row_1_leq_0: 5e-06 z_p_base_1[1] + [
   p[2] * z_p_base_1[1] ] <= 0
 infeasibility_p_base_1_row_2_leq_0: 5e-06 z_p_base_1[2] + [
   0.5 p[0] * z_p_base_1[2] - 0.5 p[1] * z_p_base_1[2]
   - 0.5 p[2] * z_p_base_1[2] + 0.5 p[3] * z_p_base_1[2]
   - 0.5 p[4] * z_p_base_1[2] ] <= 0
 infeasibility_p_base_1_row_3_leq_0: 5e-06 z_p_base_1[3] + [
   p[4] * z_p_base_1[3] ] <= 0
 infeasibility_p_base_1_row_4_leq_0: 5e-06 z_p_base_1[4] + [
   - 0.5 p[0] * z_p_base_1[4] + 0.5 p[1] * z_p_base_1[4]
   + 0.5 p[2] * z_p_base_1[4] + 0.5 p[3] * z_p_base_1[4]
   + 0.5 p[4] * z_p_base_1[4] ] <= 0
 infeasibility_p_prime_base_1_row_0_leq_0: 5e-06 z_p_prime_base_1[0] + [
   0.5 alpha[0] * z_p_prime_base_1[0] - 0.5 alpha[1] * z_p_prime_base_1[0]
   + 0.5 alpha[2] * z_p_prime_base_1[0]
   + 0.5 alpha[3] * z_p_prime_base_1[0]
   - 0.5 alpha[4] * z_p_prime_base_1[0] + 0.5 p[0] * z_p_prime_base_1[0]
   + 0.5 p[1] * z_p_prime_base_1[0] - 0.5 p[2] * z_p_prime_base_1[0]
   - 0.5 p[3] * z_p_prime_base_1[0] - 0.5 p[4] * z_p_prime_base_1[0] ]
   <= 0
 infeasibility_p_prime_base_1_row_1_leq_0: 5e-06 z_p_prime_base_1[1] + [
   - alpha[2] * z_p_prime_base_1[1] + p[2] * z_p_prime_base_1[1] ] <= 0
 infeasibility_p_prime_base_1_row_2_leq_0: 5e-06 z_p_prime_base_1[2] + [
   0.5 alpha[0] * z_p_prime_base_1[2] + 0.5 alpha[1] * z_p_prime_base_1[2]
   + 0.5 alpha[2] * z_p_prime_base_1[2]
   - 0.5 alpha[3] * z_p_prime_base_1[2]
   - 0.5 alpha[4] * z_p_prime_base_1[2] + 0.5 p[0] * z_p_prime_base_1[2]
   - 0.5 p[1] * z_p_prime_base_1[2] - 0.5 p[2] * z_p_prime_base_1[2]
   + 0.5 p[3] * z_p_prime_base_1[2] - 0.5 p[4] * z_p_prime_base_1[2] ]
   <= 0
 infeasibility_p_prime_base_1_row_3_leq_0: 5e-06 z_p_prime_base_1[3] + [
   alpha[4] * z_p_prime_base_1[3] + p[4] * z_p_prime_base_1[3] ] <= 0
 infeasibility_p_prime_base_1_row_4_leq_0: 5e-06 z_p_prime_base_1[4] + [
   - 0.5 alpha[0] * z_p_prime_base_1[4]
   - 0.5 alpha[1] * z_p_prime_base_1[4]
   - 0.5 alpha[2] * z_p_prime_base_1[4]
   - 0.5 alpha[3] * z_p_prime_base_1[4]
   + 0.5 alpha[4] * z_p_prime_base_1[4] - 0.5 p[0] * z_p_prime_base_1[4]
   + 0.5 p[1] * z_p_prime_base_1[4] + 0.5 p[2] * z_p_prime_base_1[4]
   + 0.5 p[3] * z_p_prime_base_1[4] + 0.5 p[4] * z_p_prime_base_1[4] ]
   <= 0
 infeasibility_p_base_2_row_0_leq_0: 5e-06 z_p_base_2[0] + [
   0.5 p[0] * z_p_base_2[0] + 0.5 p[1] * z_p_base_2[0]
   - 0.5 p[2] * z_p_base_2[0] - 0.5 p[3] * z_p_base_2[0]
   - 0.5 p[4] * z_p_base_2[0] ] <= 0
 infeasibility_p_base_2_row_1_leq_0: 5e-06 z_p_base_2[1] + [
   p[2] * z_p_base_2[1] ] <= 0
 infeasibility_p_base_2_row_2_leq_0: 5e-06 z_p_base_2[2] + [
   p[3] * z_p_base_2[2] ] <= 0
 infeasibility_p_base_2_row_3_leq_0: 5e-06 z_p_base_2[3] + [
   0.5 p[0] * z_p_base_2[3] - 0.5 p[1] * z_p_base_2[3]
   - 0.5 p[2] * z_p_base_2[3] - 0.5 p[3] * z_p_base_2[3]
   + 0.5 p[4] * z_p_base_2[3] ] <= 0
 infeasibility_p_base_2_row_4_leq_0: 5e-06 z_p_base_2[4] + [
   - 0.5 p[0] * z_p_base_2[4] + 0.5 p[1] * z_p_base_2[4]
   + 0.5 p[2] * z_p_base_2[4] + 0.5 p[3] * z_p_base_2[4]
   + 0.5 p[4] * z_p_base_2[4] ] <= 0
 infeasibility_p_prime_base_2_row_0_leq_0: 5e-06 z_p_prime_base_2[0] + [
   0.5 alpha[0] * z_p_prime_base_2[0] - 0.5 alpha[1] * z_p_prime_base_2[0]
   + 0.5 alpha[2] * z_p_prime_base_2[0]
   + 0.5 alpha[3] * z_p_prime_base_2[0]
   - 0.5 alpha[4] * z_p_prime_base_2[0] + 0.5 p[0] * z_p_prime_base_2[0]
   + 0.5 p[1] * z_p_prime_base_2[0] - 0.5 p[2] * z_p_prime_base_2[0]
   - 0.5 p[3] * z_p_prime_base_2[0] - 0.5 p[4] * z_p_prime_base_2[0] ]
   <= 0
 infeasibility_p_prime_base_2_row_1_leq_0: 5e-06 z_p_prime_base_2[1] + [
   - alpha[2] * z_p_prime_base_2[1] + p[2] * z_p_prime_base_2[1] ] <= 0
 infeasibility_p_prime_base_2_row_2_leq_0: 5e-06 z_p_prime_base_2[2] + [
   - alpha[3] * z_p_prime_base_2[2] + p[3] * z_p_prime_base_2[2] ] <= 0
 infeasibility_p_prime_base_2_row_3_leq_0: 5e-06 z_p_prime_base_2[3] + [
   0.5 alpha[0] * z_p_prime_base_2[3] + 0.5 alpha[1] * z_p_prime_base_2[3]
   + 0.5 alpha[2] * z_p_prime_base_2[3]
   + 0.5 alpha[3] * z_p_prime_base_2[3]
   + 0.5 alpha[4] * z_p_prime_base_2[3] + 0.5 p[0] * z_p_prime_base_2[3]
   - 0.5 p[1] * z_p_prime_base_2[3] - 0.5 p[2] * z_p_prime_base_2[3]
   - 0.5 p[3] * z_p_prime_base_2[3] + 0.5 p[4] * z_p_prime_base_2[3] ]
   <= 0
 infeasibility_p_prime_base_2_row_4_leq_0: 5e-06 z_p_prime_base_2[4] + [
   - 0.5 alpha[0] * z_p_prime_base_2[4]
   - 0.5 alpha[1] * z_p_prime_base_2[4]
   - 0.5 alpha[2] * z_p_prime_base_2[4]
   - 0.5 alpha[3] * z_p_prime_base_2[4]
   + 0.5 alpha[4] * z_p_prime_base_2[4] - 0.5 p[0] * z_p_prime_base_2[4]
   + 0.5 p[1] * z_p_prime_base_2[4] + 0.5 p[2] * z_p_prime_base_2[4]
   + 0.5 p[3] * z_p_prime_base_2[4] + 0.5 p[4] * z_p_prime_base_2[4] ]
   <= 0
 infeasibility_p_base_3_row_0_leq_0: 5e-06 z_p_base_3[0] + [
   p[1] * z_p_base_3[0] ] <= 0
 infeasibility_p_base_3_row_1_leq_0: 5e-06 z_p_base_3[1] + [
   0.5 p[0] * z_p_base_3[1] - 0.5 p[1] * z_p_base_3[1]
   + 0.5 p[2] * z_p_base_3[1] - 0.5 p[3] * z_p_base_3[1]
   - 0.5 p[4] * z_p_base_3[1] ] <= 0
 infeasibility_p_base_3_row_2_leq_0: 5e-06 z_p_base_3[2] + [
   0.5 p[0] * z_p_base_3[2] - 0.5 p[1] * z_p_base_3[2]
   - 0.5 p[2] * z_p_base_3[2] + 0.5 p[3] * z_p_base_3[2]
   - 0.5 p[4] * z_p_base_3[2] ] <= 0
 infeasibility_p_base_3_row_3_leq_0: 5e-06 z_p_base_3[3] + [
   p[4] * z_p_base_3[3] ] <= 0
 infeasibility_p_base_3_row_4_leq_0: 5e-06 z_p_base_3[4] + [
   - 0.5 p[0] * z_p_base_3[4] + 0.5 p[1] * z_p_base_3[4]
   + 0.5 p[2] * z_p_base_3[4] + 0.5 p[3] * z_p_base_3[4]
   + 0.5 p[4] * z_p_base_3[4] ] <= 0
 infeasibility_p_prime_base_3_row_0_leq_0: 5e-06 z_p_prime_base_3[0] + [
   - alpha[1] * z_p_prime_base_3[0] + p[1] * z_p_prime_base_3[0] ] <= 0
 infeasibility_p_prime_base_3_row_1_leq_0: 5e-06 z_p_prime_base_3[1] + [
   0.5 alpha[0] * z_p_prime_base_3[1] + 0.5 alpha[1] * z_p_prime_base_3[1]
   - 0.5 alpha[2] * z_p_prime_base_3[1]
   + 0.5 alpha[3] * z_p_prime_base_3[1]
   - 0.5 alpha[4] * z_p_prime_base_3[1] + 0.5 p[0] * z_p_prime_base_3[1]
   - 0.5 p[1] * z_p_prime_base_3[1] + 0.5 p[2] * z_p_prime_base_3[1]
   - 0.5 p[3] * z_p_prime_base_3[1] - 0.5 p[4] * z_p_prime_base_3[1] ]
   <= 0
 infeasibility_p_prime_base_3_row_2_leq_0: 5e-06 z_p_prime_base_3[2] + [
   0.5 alpha[0] * z_p_prime_base_3[2] + 0.5 alpha[1] * z_p_prime_base_3[2]
   + 0.5 alpha[2] * z_p_prime_base_3[2]
   - 0.5 alpha[3] * z_p_prime_base_3[2]
   - 0.5 alpha[4] * z_p_prime_base_3[2] + 0.5 p[0] * z_p_prime_base_3[2]
   - 0.5 p[1] * z_p_prime_base_3[2] - 0.5 p[2] * z_p_prime_base_3[2]
   + 0.5 p[3] * z_p_prime_base_3[2] - 0.5 p[4] * z_p_prime_base_3[2] ]
   <= 0
 infeasibility_p_prime_base_3_row_3_leq_0: 5e-06 z_p_prime_base_3[3] + [
   alpha[4] * z_p_prime_base_3[3] + p[4] * z_p_prime_base_3[3] ] <= 0
 infeasibility_p_prime_base_3_row_4_leq_0: 5e-06 z_p_prime_base_3[4] + [
   - 0.5 alpha[0] * z_p_prime_base_3[4]
   - 0.5 alpha[1] * z_p_prime_base_3[4]
   - 0.5 alpha[2] * z_p_prime_base_3[4]
   - 0.5 alpha[3] * z_p_prime_base_3[4]
   + 0.5 alpha[4] * z_p_prime_base_3[4] - 0.5 p[0] * z_p_prime_base_3[4]
   + 0.5 p[1] * z_p_prime_base_3[4] + 0.5 p[2] * z_p_prime_base_3[4]
   + 0.5 p[3] * z_p_prime_base_3[4] + 0.5 p[4] * z_p_prime_base_3[4] ]
   <= 0
 infeasibility_p_base_4_row_0_leq_0: 5e-06 z_p_base_4[0] + [
   p[1] * z_p_base_4[0] ] <= 0
 infeasibility_p_base_4_row_1_leq_0: 5e-06 z_p_base_4[1] + [
   0.5 p[0] * z_p_base_4[1] - 0.5 p[1] * z_p_base_4[1]
   + 0.5 p[2] * z_p_base_4[1] - 0.5 p[3] * z_p_base_4[1]
   - 0.5 p[4] * z_p_base_4[1] ] <= 0
 infeasibility_p_base_4_row_2_leq_0: 5e-06 z_p_base_4[2] + [
   p[3] * z_p_base_4[2] ] <= 0
 infeasibility_p_base_4_row_3_leq_0: 5e-06 z_p_base_4[3] + [
   0.5 p[0] * z_p_base_4[3] - 0.5 p[1] * z_p_base_4[3]
   - 0.5 p[2] * z_p_base_4[3] - 0.5 p[3] * z_p_base_4[3]
   + 0.5 p[4] * z_p_base_4[3] ] <= 0
 infeasibility_p_base_4_row_4_leq_0: 5e-06 z_p_base_4[4] + [
   - 0.5 p[0] * z_p_base_4[4] + 0.5 p[1] * z_p_base_4[4]
   + 0.5 p[2] * z_p_base_4[4] + 0.5 p[3] * z_p_base_4[4]
   + 0.5 p[4] * z_p_base_4[4] ] <= 0
 infeasibility_p_prime_base_4_row_0_leq_0: 5e-06 z_p_prime_base_4[0] + [
   - alpha[1] * z_p_prime_base_4[0] + p[1] * z_p_prime_base_4[0] ] <= 0
 infeasibility_p_prime_base_4_row_1_leq_0: 5e-06 z_p_prime_base_4[1] + [
   0.5 alpha[0] * z_p_prime_base_4[1] + 0.5 alpha[1] * z_p_prime_base_4[1]
   - 0.5 alpha[2] * z_p_prime_base_4[1]
   + 0.5 alpha[3] * z_p_prime_base_4[1]
   - 0.5 alpha[4] * z_p_prime_base_4[1] + 0.5 p[0] * z_p_prime_base_4[1]
   - 0.5 p[1] * z_p_prime_base_4[1] + 0.5 p[2] * z_p_prime_base_4[1]
   - 0.5 p[3] * z_p_prime_base_4[1] - 0.5 p[4] * z_p_prime_base_4[1] ]
   <= 0
 infeasibility_p_prime_base_4_row_2_leq_0: 5e-06 z_p_prime_base_4[2] + [
   - alpha[3] * z_p_prime_base_4[2] + p[3] * z_p_prime_base_4[2] ] <= 0
 infeasibility_p_prime_base_4_row_3_leq_0: 5e-06 z_p_prime_base_4[3] + [
   0.5 alpha[0] * z_p_prime_base_4[3] + 0.5 alpha[1] * z_p_prime_base_4[3]
   + 0.5 alpha[2] * z_p_prime_base_4[3]
   + 0.5 alpha[3] * z_p_prime_base_4[3]
   + 0.5 alpha[4] * z_p_prime_base_4[3] + 0.5 p[0] * z_p_prime_base_4[3]
   - 0.5 p[1] * z_p_prime_base_4[3] - 0.5 p[2] * z_p_prime_base_4[3]
   - 0.5 p[3] * z_p_prime_base_4[3] + 0.5 p[4] * z_p_prime_base_4[3] ]
   <= 0
 infeasibility_p_prime_base_4_row_4_leq_0: 5e-06 z_p_prime_base_4[4] + [
   - 0.5 alpha[0] * z_p_prime_base_4[4]
   - 0.5 alpha[1] * z_p_prime_base_4[4]
   - 0.5 alpha[2] * z_p_prime_base_4[4]
   - 0.5 alpha[3] * z_p_prime_base_4[4]
   + 0.5 alpha[4] * z_p_prime_base_4[4] - 0.5 p[0] * z_p_prime_base_4[4]
   + 0.5 p[1] * z_p_prime_base_4[4] + 0.5 p[2] * z_p_prime_base_4[4]
   + 0.5 p[3] * z_p_prime_base_4[4] + 0.5 p[4] * z_p_prime_base_4[4] ]
   <= 0
 infeasibility_p_base_5_row_0_leq_0: 5e-06 z_p_base_5[0] + [
   p[1] * z_p_base_5[0] ] <= 0
 infeasibility_p_base_5_row_1_leq_0: 5e-06 z_p_base_5[1] + [
   p[2] * z_p_base_5[1] ] <= 0
 infeasibility_p_base_5_row_2_leq_0: 5e-06 z_p_base_5[2] + [
   0.5 p[0] * z_p_base_5[2] - 0.5 p[1] * z_p_base_5[2]
   - 0.5 p[2] * z_p_base_5[2] + 0.5 p[3] * z_p_base_5[2]
   - 0.5 p[4] * z_p_base_5[2] ] <= 0
 infeasibility_p_base_5_row_3_leq_0: 5e-06 z_p_base_5[3] + [
   0.5 p[0] * z_p_base_5[3] - 0.5 p[1] * z_p_base_5[3]
   - 0.5 p[2] * z_p_base_5[3] - 0.5 p[3] * z_p_base_5[3]
   + 0.5 p[4] * z_p_base_5[3] ] <= 0
 infeasibility_p_base_5_row_4_leq_0: 5e-06 z_p_base_5[4] + [
   - 0.5 p[0] * z_p_base_5[4] + 0.5 p[1] * z_p_base_5[4]
   + 0.5 p[2] * z_p_base_5[4] + 0.5 p[3] * z_p_base_5[4]
   + 0.5 p[4] * z_p_base_5[4] ] <= 0
 infeasibility_p_prime_base_5_row_0_leq_0: 5e-06 z_p_prime_base_5[0] + [
   - alpha[1] * z_p_prime_base_5[0] + p[1] * z_p_prime_base_5[0] ] <= 0
 infeasibility_p_prime_base_5_row_1_leq_0: 5e-06 z_p_prime_base_5[1] + [
   - alpha[2] * z_p_prime_base_5[1] + p[2] * z_p_prime_base_5[1] ] <= 0
 infeasibility_p_prime_base_5_row_2_leq_0: 5e-06 z_p_prime_base_5[2] + [
   0.5 alpha[0] * z_p_prime_base_5[2] + 0.5 alpha[1] * z_p_prime_base_5[2]
   + 0.5 alpha[2] * z_p_prime_base_5[2]
   - 0.5 alpha[3] * z_p_prime_base_5[2]
   - 0.5 alpha[4] * z_p_prime_base_5[2] + 0.5 p[0] * z_p_prime_base_5[2]
   - 0.5 p[1] * z_p_prime_base_5[2] - 0.5 p[2] * z_p_prime_base_5[2]
   + 0.5 p[3] * z_p_prime_base_5[2] - 0.5 p[4] * z_p_prime_base_5[2] ]
   <= 0
 infeasibility_p_prime_base_5_row_3_leq_0: 5e-06 z_p_prime_base_5[3] + [
   0.5 alpha[0] * z_p_prime_base_5[3] + 0.5 alpha[1] * z_p_prime_base_5[3]
   + 0.5 alpha[2] * z_p_prime_base_5[3]
   + 0.5 alpha[3] * z_p_prime_base_5[3]
   + 0.5 alpha[4] * z_p_prime_base_5[3] + 0.5 p[0] * z_p_prime_base_5[3]
   - 0.5 p[1] * z_p_prime_base_5[3] - 0.5 p[2] * z_p_prime_base_5[3]
   - 0.5 p[3] * z_p_prime_base_5[3] + 0.5 p[4] * z_p_prime_base_5[3] ]
   <= 0
 infeasibility_p_prime_base_5_row_4_leq_0: 5e-06 z_p_prime_base_5[4] + [
   - 0.5 alpha[0] * z_p_prime_base_5[4]
   - 0.5 alpha[1] * z_p_prime_base_5[4]
   - 0.5 alpha[2] * z_p_prime_base_5[4]
   - 0.5 alpha[3] * z_p_prime_base_5[4]
   + 0.5 alpha[4] * z_p_prime_base_5[4] - 0.5 p[0] * z_p_prime_base_5[4]
   + 0.5 p[1] * z_p_prime_base_5[4] + 0.5 p[2] * z_p_prime_base_5[4]
   + 0.5 p[3] * z_p_prime_base_5[4] + 0.5 p[4] * z_p_prime_base_5[4] ]
   <= 0
 infeasibility_p_base_6_row_0_leq_0: 5e-06 z_p_base_6[0] + [
   0.5 p[0] * z_p_base_6[0] + 0.5 p[1] * z_p_base_6[0]
   - 0.5 p[2] * z_p_base_6[0] - 0.5 p[3] * z_p_base_6[0]
   - 0.5 p[4] * z_p_base_6[0] ] <= 0
 infeasibility_p_base_6_row_1_leq_0: 5e-06 z_p_base_6[1] + [
   0.5 p[0] * z_p_base_6[1] - 0.5 p[1] * z_p_base_6[1]
   + 0.5 p[2] * z_p_base_6[1] - 0.5 p[3] * z_p_base_6[1]
   + 0.5 p[4] * z_p_base_6[1] ] <= 0
 infeasibility_p_base_6_row_2_leq_0: 5e-06 z_p_base_6[2] + [
   p[3] * z_p_base_6[2] ] <= 0
 infeasibility_p_base_6_row_3_leq_0: 5e-06 z_p_base_6[3] + [
   - 0.5 p[0] * z_p_base_6[3] + 0.5 p[1] * z_p_base_6[3]
   + 0.5 p[2] * z_p_base_6[3] + 0.5 p[3] * z_p_base_6[3]
   - 0.5 p[4] * z_p_base_6[3] ] <= 0
 infeasibility_p_base_6_row_4_leq_0: 5e-06 z_p_base_6[4] + [
   p[4] * z_p_base_6[4] ] <= 0
 infeasibility_p_prime_base_6_row_0_leq_0: 5e-06 z_p_prime_base_6[0] + [
   0.5 alpha[0] * z_p_prime_base_6[0] - 0.5 alpha[1] * z_p_prime_base_6[0]
   + 0.5 alpha[2] * z_p_prime_base_6[0]
   + 0.5 alpha[3] * z_p_prime_base_6[0]
   - 0.5 alpha[4] * z_p_prime_base_6[0] + 0.5 p[0] * z_p_prime_base_6[0]
   + 0.5 p[1] * z_p_prime_base_6[0] - 0.5 p[2] * z_p_prime_base_6[0]
   - 0.5 p[3] * z_p_prime_base_6[0] - 0.5 p[4] * z_p_prime_base_6[0] ]
   <= 0
 infeasibility_p_prime_base_6_row_1_leq_0: 5e-06 z_p_prime_base_6[1] + [
   0.5 alpha[0] * z_p_prime_base_6[1] + 0.5 alpha[1] * z_p_prime_base_6[1]
   - 0.5 alpha[2] * z_p_prime_base_6[1]
   + 0.5 alpha[3] * z_p_prime_base_6[1]
   + 0.5 alpha[4] * z_p_prime_base_6[1] + 0.5 p[0] * z_p_prime_base_6[1]
   - 0.5 p[1] * z_p_prime_base_6[1] + 0.5 p[2] * z_p_prime_base_6[1]
   - 0.5 p[3] * z_p_prime_base_6[1] + 0.5 p[4] * z_p_prime_base_6[1] ]
   <= 0
 infeasibility_p_prime_base_6_row_2_leq_0: 5e-06 z_p_prime_base_6[2] + [
   - alpha[3] * z_p_prime_base_6[2] + p[3] * z_p_prime_base_6[2] ] <= 0
 infeasibility_p_prime_base_6_row_3_leq_0: 5e-06 z_p_prime_base_6[3] + [
   - 0.5 alpha[0] * z_p_prime_base_6[3]
   - 0.5 alpha[1] * z_p_prime_base_6[3]
   - 0.5 alpha[2] * z_p_prime_base_6[3]
   - 0.5 alpha[3] * z_p_prime_base_6[3]
   - 0.5 alpha[4] * z_p_prime_base_6[3] - 0.5 p[0] * z_p_prime_base_6[3]
   + 0.5 p[1] * z_p_prime_base_6[3] + 0.5 p[2] * z_p_prime_base_6[3]
   + 0.5 p[3] * z_p_prime_base_6[3] - 0.5 p[4] * z_p_prime_base_6[3] ]
   <= 0
 infeasibility_p_prime_base_6_row_4_leq_0: 5e-06 z_p_prime_base_6[4] + [
   alpha[4] * z_p_prime_base_6[4] + p[4] * z_p_prime_base_6[4] ] <= 0
 infeasibility_p_base_7_row_0_leq_0: 5e-06 z_p_base_7[0] + [
   0.5 p[0] * z_p_base_7[0] + 0.5 p[1] * z_p_base_7[0]
   - 0.5 p[2] * z_p_base_7[0] - 0.5 p[3] * z_p_base_7[0]
   + 0.5 p[4] * z_p_base_7[0] ] <= 0
 infeasibility_p_base_7_row_1_leq_0: 5e-06 z_p_base_7[1] + [
   0.5 p[0] * z_p_base_7[1] - 0.5 p[1] * z_p_base_7[1]
   + 0.5 p[2] * z_p_base_7[1] - 0.5 p[3] * z_p_base_7[1]
   - 0.5 p[4] * z_p_base_7[1] ] <= 0
 infeasibility_p_base_7_row_2_leq_0: 5e-06 z_p_base_7[2] + [
   p[3] * z_p_base_7[2] ] <= 0
 infeasibility_p_base_7_row_3_leq_0: 5e-06 z_p_base_7[3] + [
   - 0.5 p[0] * z_p_base_7[3] + 0.5 p[1] * z_p_base_7[3]
   + 0.5 p[2] * z_p_base_7[3] + 0.5 p[3] * z_p_base_7[3]
   - 0.5 p[4] * z_p_base_7[3] ] <= 0
 infeasibility_p_base_7_row_4_leq_0: 5e-06 z_p_base_7[4] + [
   p[4] * z_p_base_7[4] ] <= 0
 infeasibility_p_prime_base_7_row_0_leq_0: 5e-06 z_p_prime_base_7[0] + [
   0.5 alpha[0] * z_p_prime_base_7[0] - 0.5 alpha[1] * z_p_prime_base_7[0]
   + 0.5 alpha[2] * z_p_prime_base_7[0]
   + 0.5 alpha[3] * z_p_prime_base_7[0]
   + 0.5 alpha[4] * z_p_prime_base_7[0] + 0.5 p[0] * z_p_prime_base_7[0]
   + 0.5 p[1] * z_p_prime_base_7[0] - 0.5 p[2] * z_p_prime_base_7[0]
   - 0.5 p[3] * z_p_prime_base_7[0] + 0.5 p[4] * z_p_prime_base_7[0] ]
   <= 0
 infeasibility_p_prime_base_7_row_1_leq_0: 5e-06 z_p_prime_base_7[1] + [
   0.5 alpha[0] * z_p_prime_base_7[1] + 0.5 alpha[1] * z_p_prime_base_7[1]
   - 0.5 alpha[2] * z_p_prime_base_7[1]
   + 0.5 alpha[3] * z_p_prime_base_7[1]
   - 0.5 alpha[4] * z_p_prime_base_7[1] + 0.5 p[0] * z_p_prime_base_7[1]
   - 0.5 p[1] * z_p_prime_base_7[1] + 0.5 p[2] * z_p_prime_base_7[1]
   - 0.5 p[3] * z_p_prime_base_7[1] - 0.5 p[4] * z_p_prime_base_7[1] ]
   <= 0
 infeasibility_p_prime_base_7_row_2_leq_0: 5e-06 z_p_prime_base_7[2] + [
   - alpha[3] * z_p_prime_base_7[2] + p[3] * z_p_prime_base_7[2] ] <= 0
 infeasibility_p_prime_base_7_row_3_leq_0: 5e-06 z_p_prime_base_7[3] + [
   - 0.5 alpha[0] * z_p_prime_base_7[3]
   - 0.5 alpha[1] * z_p_prime_base_7[3]
   - 0.5 alpha[2] * z_p_prime_base_7[3]
   - 0.5 alpha[3] * z_p_prime_base_7[3]
   - 0.5 alpha[4] * z_p_prime_base_7[3] - 0.5 p[0] * z_p_prime_base_7[3]
   + 0.5 p[1] * z_p_prime_base_7[3] + 0.5 p[2] * z_p_prime_base_7[3]
   + 0.5 p[3] * z_p_prime_base_7[3] - 0.5 p[4] * z_p_prime_base_7[3] ]
   <= 0
 infeasibility_p_prime_base_7_row_4_leq_0: 5e-06 z_p_prime_base_7[4] + [
   alpha[4] * z_p_prime_base_7[4] + p[4] * z_p_prime_base_7[4] ] <= 0
 infeasibility_p_base_8_row_0_leq_0: 5e-06 z_p_base_8[0] + [
   0.5 p[0] * z_p_base_8[0] + 0.5 p[1] * z_p_base_8[0]
   - 0.5 p[2] * z_p_base_8[0] - 0.5 p[3] * z_p_base_8[0]
   + 0.5 p[4] * z_p_base_8[0] ] <= 0
 infeasibility_p_base_8_row_1_leq_0: 5e-06 z_p_base_8[1] + [
   0.5 p[0] * z_p_base_8[1] - 0.5 p[1] * z_p_base_8[1]
   + 0.5 p[2] * z_p_base_8[1] - 0.5 p[3] * z_p_base_8[1]
   + 0.5 p[4] * z_p_base_8[1] ] <= 0
 infeasibility_p_base_8_row_2_leq_0: 5e-06 z_p_base_8[2] + [
   p[3] * z_p_base_8[2] - p[4] * z_p_base_8[2] ] <= 0
 infeasibility_p_base_8_row_3_leq_0: 5e-06 z_p_base_8[3] + [
   - 0.5 p[0] * z_p_base_8[3] + 0.5 p[1] * z_p_base_8[3]
   + 0.5 p[2] * z_p_base_8[3] + 0.5 p[3] * z_p_base_8[3]
   - 0.5 p[4] * z_p_base_8[3] ] <= 0
 infeasibility_p_base_8_row_4_leq_0: 5e-06 z_p_base_8[4] + [
   p[4] * z_p_base_8[4] ] <= 0
 infeasibility_p_prime_base_8_row_0_leq_0: 5e-06 z_p_prime_base_8[0] + [
   0.5 alpha[0] * z_p_prime_base_8[0] - 0.5 alpha[1] * z_p_prime_base_8[0]
   + 0.5 alpha[2] * z_p_prime_base_8[0]
   + 0.5 alpha[3] * z_p_prime_base_8[0]
   + 0.5 alpha[4] * z_p_prime_base_8[0] + 0.5 p[0] * z_p_prime_base_8[0]
   + 0.5 p[1] * z_p_prime_base_8[0] - 0.5 p[2] * z_p_prime_base_8[0]
   - 0.5 p[3] * z_p_prime_base_8[0] + 0.5 p[4] * z_p_prime_base_8[0] ]
   <= 0
 infeasibility_p_prime_base_8_row_1_leq_0: 5e-06 z_p_prime_base_8[1] + [
   0.5 alpha[0] * z_p_prime_base_8[1] + 0.5 alpha[1] * z_p_prime_base_8[1]
   - 0.5 alpha[2] * z_p_prime_base_8[1]
   + 0.5 alpha[3] * z_p_prime_base_8[1]
   + 0.5 alpha[4] * z_p_prime_base_8[1] + 0.5 p[0] * z_p_prime_base_8[1]
   - 0.5 p[1] * z_p_prime_base_8[1] + 0.5 p[2] * z_p_prime_base_8[1]
   - 0.5 p[3] * z_p_prime_base_8[1] + 0.5 p[4] * z_p_prime_base_8[1] ]
   <= 0
 infeasibility_p_prime_base_8_row_2_leq_0: 5e-06 z_p_prime_base_8[2] + [
   - alpha[3] * z_p_prime_base_8[2] - alpha[4] * z_p_prime_base_8[2]
   + p[3] * z_p_prime_base_8[2] - p[4] * z_p_prime_base_8[2] ] <= 0
 infeasibility_p_prime_base_8_row_3_leq_0: 5e-06 z_p_prime_base_8[3] + [
   - 0.5 alpha[0] * z_p_prime_base_8[3]
   - 0.5 alpha[1] * z_p_prime_base_8[3]
   - 0.5 alpha[2] * z_p_prime_base_8[3]
   - 0.5 alpha[3] * z_p_prime_base_8[3]
   - 0.5 alpha[4] * z_p_prime_base_8[3] - 0.5 p[0] * z_p_prime_base_8[3]
   + 0.5 p[1] * z_p_prime_base_8[3] + 0.5 p[2] * z_p_prime_base_8[3]
   + 0.5 p[3] * z_p_prime_base_8[3] - 0.5 p[4] * z_p_prime_base_8[3] ]
   <= 0
 infeasibility_p_prime_base_8_row_4_leq_0: 5e-06 z_p_prime_base_8[4] + [
   alpha[4] * z_p_prime_base_8[4] + p[4] * z_p_prime_base_8[4] ] <= 0
 infeasibility_p_base_9_row_0_leq_0: 5e-06 z_p_base_9[0] + [
   0.5 p[0] * z_p_base_9[0] + 0.5 p[1] * z_p_base_9[0]
   - 0.5 p[2] * z_p_base_9[0] - 0.5 p[3] * z_p_base_9[0]
   - 0.5 p[4] * z_p_base_9[0] ] <= 0
 infeasibility_p_base_9_row_1_leq_0: 5e-06 z_p_base_9[1] + [
   p[2] * z_p_base_9[1] ] <= 0
 infeasibility_p_base_9_row_2_leq_0: 5e-06 z_p_base_9[2] + [
   0.5 p[0] * z_p_base_9[2] - 0.5 p[1] * z_p_base_9[2]
   - 0.5 p[2] * z_p_base_9[2] + 0.5 p[3] * z_p_base_9[2]
   + 0.5 p[4] * z_p_base_9[2] ] <= 0
 infeasibility_p_base_9_row_3_leq_0: 5e-06 z_p_base_9[3] + [
   - 0.5 p[0] * z_p_base_9[3] + 0.5 p[1] * z_p_base_9[3]
   + 0.5 p[2] * z_p_base_9[3] + 0.5 p[3] * z_p_base_9[3]
   - 0.5 p[4] * z_p_base_9[3] ] <= 0
 infeasibility_p_base_9_row_4_leq_0: 5e-06 z_p_base_9[4] + [
   p[4] * z_p_base_9[4] ] <= 0
 infeasibility_p_prime_base_9_row_0_leq_0: 5e-06 z_p_prime_base_9[0] + [
   0.5 alpha[0] * z_p_prime_base_9[0] - 0.5 alpha[1] * z_p_prime_base_9[0]
   + 0.5 alpha[2] * z_p_prime_base_9[0]
   + 0.5 alpha[3] * z_p_prime_base_9[0]
   - 0.5 alpha[4] * z_p_prime_base_9[0] + 0.5 p[0] * z_p_prime_base_9[0]
   + 0.5 p[1] * z_p_prime_base_9[0] - 0.5 p[2] * z_p_prime_base_9[0]
   - 0.5 p[3] * z_p_prime_base_9[0] - 0.5 p[4] * z_p_prime_base_9[0] ]
   <= 0
 infeasibility_p_prime_base_9_row_1_leq_0: 5e-06 z_p_prime_base_9[1] + [
   - alpha[2] * z_p_prime_base_9[1] + p[2] * z_p_prime_base_9[1] ] <= 0
 infeasibility_p_prime_base_9_row_2_leq_0: 5e-06 z_p_prime_base_9[2] + [
   0.5 alpha[0] * z_p_prime_base_9[2] + 0.5 alpha[1] * z_p_prime_base_9[2]
   + 0.5 alpha[2] * z_p_prime_base_9[2]
   - 0.5 alpha[3] * z_p_prime_base_9[2]
   + 0.5 alpha[4] * z_p_prime_base_9[2] + 0.5 p[0] * z_p_prime_base_9[2]
   - 0.5 p[1] * z_p_prime_base_9[2] - 0.5 p[2] * z_p_prime_base_9[2]
   + 0.5 p[3] * z_p_prime_base_9[2] + 0.5 p[4] * z_p_prime_base_9[2] ]
   <= 0
 infeasibility_p_prime_base_9_row_3_leq_0: 5e-06 z_p_prime_base_9[3] + [
   - 0.5 alpha[0] * z_p_prime_base_9[3]
   - 0.5 alpha[1] * z_p_prime_base_9[3]
   - 0.5 alpha[2] * z_p_prime_base_9[3]
   - 0.5 alpha[3] * z_p_prime_base_9[3]
   - 0.5 alpha[4] * z_p_prime_base_9[3] - 0.5 p[0] * z_p_prime_base_9[3]
   + 0.5 p[1] * z_p_prime_base_9[3] + 0.5 p[2] * z_p_prime_base_9[3]
   + 0.5 p[3] * z_p_prime_base_9[3] - 0.5 p[4] * z_p_prime_base_9[3] ]
   <= 0
 infeasibility_p_prime_base_9_row_4_leq_0: 5e-06 z_p_prime_base_9[4] + [
   alpha[4] * z_p_prime_base_9[4] + p[4] * z_p_prime_base_9[4] ] <= 0
 infeasibility_p_base_10_row_0_leq_0: 5e-06 z_p_base_10[0] + [
   0.5 p[0] * z_p_base_10[0] + 0.5 p[1] * z_p_base_10[0]
   - 0.5 p[2] * z_p_base_10[0] - 0.5 p[3] * z_p_base_10[0]
   + 0.5 p[4] * z_p_base_10[0] ] <= 0
 infeasibility_p_base_10_row_1_leq_0: 5e-06 z_p_base_10[1] + [
   p[2] * z_p_base_10[1] - p[4] * z_p_base_10[1] ] <= 0
 infeasibility_p_base_10_row_2_leq_0: 5e-06 z_p_base_10[2] + [
   0.5 p[0] * z_p_base_10[2] - 0.5 p[1] * z_p_base_10[2]
   - 0.5 p[2] * z_p_base_10[2] + 0.5 p[3] * z_p_base_10[2]
   + 0.5 p[4] * z_p_base_10[2] ] <= 0
 infeasibility_p_base_10_row_3_leq_0: 5e-06 z_p_base_10[3] + [
   - 0.5 p[0] * z_p_base_10[3] + 0.5 p[1] * z_p_base_10[3]
   + 0.5 p[2] * z_p_base_10[3] + 0.5 p[3] * z_p_base_10[3]
   - 0.5 p[4] * z_p_base_10[3] ] <= 0
 infeasibility_p_base_10_row_4_leq_0: 5e-06 z_p_base_10[4] + [
   p[4] * z_p_base_10[4] ] <= 0
 infeasibility_p_prime_base_10_row_0_leq_0: 5e-06 z_p_prime_base_10[0] + [
   0.5 alpha[0] * z_p_prime_base_10[0]
   - 0.5 alpha[1] * z_p_prime_base_10[0]
   + 0.5 alpha[2] * z_p_prime_base_10[0]
   + 0.5 alpha[3] * z_p_prime_base_10[0]
   + 0.5 alpha[4] * z_p_prime_base_10[0] + 0.5 p[0] * z_p_prime_base_10[0]
   + 0.5 p[1] * z_p_prime_base_10[0] - 0.5 p[2] * z_p_prime_base_10[0]
   - 0.5 p[3] * z_p_prime_base_10[0] + 0.5 p[4] * z_p_prime_base_10[0] ]
   <= 0
 infeasibility_p_prime_base_10_row_1_leq_0: 5e-06 z_p_prime_base_10[1] + [
   - alpha[2] * z_p_prime_base_10[1] - alpha[4] * z_p_prime_base_10[1]
   + p[2] * z_p_prime_base_10[1] - p[4] * z_p_prime_base_10[1] ] <= 0
 infeasibility_p_prime_base_10_row_2_leq_0: 5e-06 z_p_prime_base_10[2] + [
   0.5 alpha[0] * z_p_prime_base_10[2]
   + 0.5 alpha[1] * z_p_prime_base_10[2]
   + 0.5 alpha[2] * z_p_prime_base_10[2]
   - 0.5 alpha[3] * z_p_prime_base_10[2]
   + 0.5 alpha[4] * z_p_prime_base_10[2] + 0.5 p[0] * z_p_prime_base_10[2]
   - 0.5 p[1] * z_p_prime_base_10[2] - 0.5 p[2] * z_p_prime_base_10[2]
   + 0.5 p[3] * z_p_prime_base_10[2] + 0.5 p[4] * z_p_prime_base_10[2] ]
   <= 0
 infeasibility_p_prime_base_10_row_3_leq_0: 5e-06 z_p_prime_base_10[3] + [
   - 0.5 alpha[0] * z_p_prime_base_10[3]
   - 0.5 alpha[1] * z_p_prime_base_10[3]
   - 0.5 alpha[2] * z_p_prime_base_10[3]
   - 0.5 alpha[3] * z_p_prime_base_10[3]
   - 0.5 alpha[4] * z_p_prime_base_10[3] - 0.5 p[0] * z_p_prime_base_10[3]
   + 0.5 p[1] * z_p_prime_base_10[3] + 0.5 p[2] * z_p_prime_base_10[3]
   + 0.5 p[3] * z_p_prime_base_10[3] - 0.5 p[4] * z_p_prime_base_10[3] ]
   <= 0
 infeasibility_p_prime_base_10_row_4_leq_0: 5e-06 z_p_prime_base_10[4] + [
   alpha[4] * z_p_prime_base_10[4] + p[4] * z_p_prime_base_10[4] ] <= 0
 infeasibility_p_base_11_row_0_leq_0: 5e-06 z_p_base_11[0] + [
   0.5 p[0] * z_p_base_11[0] + 0.5 p[1] * z_p_base_11[0]
   - 0.5 p[2] * z_p_base_11[0] - 0.5 p[3] * z_p_base_11[0]
   + 0.5 p[4] * z_p_base_11[0] ] <= 0
 infeasibility_p_base_11_row_1_leq_0: 5e-06 z_p_base_11[1] + [
   p[2] * z_p_base_11[1] ] <= 0
 infeasibility_p_base_11_row_2_leq_0: 5e-06 z_p_base_11[2] + [
   0.5 p[0] * z_p_base_11[2] - 0.5 p[1] * z_p_base_11[2]
   - 0.5 p[2] * z_p_base_11[2] + 0.5 p[3] * z_p_base_11[2]
   - 0.5 p[4] * z_p_base_11[2] ] <= 0
 infeasibility_p_base_11_row_3_leq_0: 5e-06 z_p_base_11[3] + [
   - 0.5 p[0] * z_p_base_11[3] + 0.5 p[1] * z_p_base_11[3]
   + 0.5 p[2] * z_p_base_11[3] + 0.5 p[3] * z_p_base_11[3]
   - 0.5 p[4] * z_p_base_11[3] ] <= 0
 infeasibility_p_base_11_row_4_leq_0: 5e-06 z_p_base_11[4] + [
   p[4] * z_p_base_11[4] ] <= 0
 infeasibility_p_prime_base_11_row_0_leq_0: 5e-06 z_p_prime_base_11[0] + [
   0.5 alpha[0] * z_p_prime_base_11[0]
   - 0.5 alpha[1] * z_p_prime_base_11[0]
   + 0.5 alpha[2] * z_p_prime_base_11[0]
   + 0.5 alpha[3] * z_p_prime_base_11[0]
   + 0.5 alpha[4] * z_p_prime_base_11[0] + 0.5 p[0] * z_p_prime_base_11[0]
   + 0.5 p[1] * z_p_prime_base_11[0] - 0.5 p[2] * z_p_prime_base_11[0]
   - 0.5 p[3] * z_p_prime_base_11[0] + 0.5 p[4] * z_p_prime_base_11[0] ]
   <= 0
 infeasibility_p_prime_base_11_row_1_leq_0: 5e-06 z_p_prime_base_11[1] + [
   - alpha[2] * z_p_prime_base_11[1] + p[2] * z_p_prime_base_11[1] ] <= 0
 infeasibility_p_prime_base_11_row_2_leq_0: 5e-06 z_p_prime_base_11[2] + [
   0.5 alpha[0] * z_p_prime_base_11[2]
   + 0.5 alpha[1] * z_p_prime_base_11[2]
   + 0.5 alpha[2] * z_p_prime_base_11[2]
   - 0.5 alpha[3] * z_p_prime_base_11[2]
   - 0.5 alpha[4] * z_p_prime_base_11[2] + 0.5 p[0] * z_p_prime_base_11[2]
   - 0.5 p[1] * z_p_prime_base_11[2] - 0.5 p[2] * z_p_prime_base_11[2]
   + 0.5 p[3] * z_p_prime_base_11[2] - 0.5 p[4] * z_p_prime_base_11[2] ]
   <= 0
 infeasibility_p_prime_base_11_row_3_leq_0: 5e-06 z_p_prime_base_11[3] + [
   - 0.5 alpha[0] * z_p_prime_base_11[3]
   - 0.5 alpha[1] * z_p_prime_base_11[3]
   - 0.5 alpha[2] * z_p_prime_base_11[3]
   - 0.5 alpha[3] * z_p_prime_base_11[3]
   - 0.5 alpha[4] * z_p_prime_base_11[3] - 0.5 p[0] * z_p_prime_base_11[3]
   + 0.5 p[1] * z_p_prime_base_11[3] + 0.5 p[2] * z_p_prime_base_11[3]
   + 0.5 p[3] * z_p_prime_base_11[3] - 0.5 p[4] * z_p_prime_base_11[3] ]
   <= 0
 infeasibility_p_prime_base_11_row_4_leq_0: 5e-06 z_p_prime_base_11[4] + [
   alpha[4] * z_p_prime_base_11[4] + p[4] * z_p_prime_base_11[4] ] <= 0
 infeasibility_p_base_12_row_0_leq_0: 5e-06 z_p_base_12[0] + [
   p[1] * z_p_base_12[0] - p[4] * z_p_base_12[0] ] <= 0
 infeasibility_p_base_12_row_1_leq_0: 5e-06 z_p_base_12[1] + [
   0.5 p[0] * z_p_base_12[1] - 0.5 p[1] * z_p_base_12[1]
   + 0.5 p[2] * z_p_base_12[1] - 0.5 p[3] * z_p_base_12[1]
   + 0.5 p[4] * z_p_base_12[1] ] <= 0
 infeasibility_p_base_12_row_2_leq_0: 5e-06 z_p_base_12[2] + [
   0.5 p[0] * z_p_base_12[2] - 0.5 p[1] * z_p_base_12[2]
   - 0.5 p[2] * z_p_base_12[2] + 0.5 p[3] * z_p_base_12[2]
   + 0.5 p[4] * z_p_base_12[2] ] <= 0
 infeasibility_p_base_12_row_3_leq_0: 5e-06 z_p_base_12[3] + [
   p[4] * z_p_base_12[3] ] <= 0
 infeasibility_p_base_12_row_4_leq_0: 5e-06 z_p_base_12[4] + [
   - 0.5 p[0] * z_p_base_12[4] + 0.5 p[1] * z_p_base_12[4]
   + 0.5 p[2] * z_p_base_12[4] + 0.5 p[3] * z_p_base_12[4]
   - 0.5 p[4] * z_p_base_12[4] ] <= 0
 infeasibility_p_prime_base_12_row_0_leq_0: 5e-06 z_p_prime_base_12[0] + [
   - alpha[1] * z_p_prime_base_12[0] - alpha[4] * z_p_prime_base_12[0]
   + p[1] * z_p_prime_base_12[0] - p[4] * z_p_prime_base_12[0] ] <= 0
 infeasibility_p_prime_base_12_row_1_leq_0: 5e-06 z_p_prime_base_12[1] + [
   0.5 alpha[0] * z_p_prime_base_12[1]
   + 0.5 alpha[1] * z_p_prime_base_12[1]
   - 0.5 alpha[2] * z_p_prime_base_12[1]
   + 0.5 alpha[3] * z_p_prime_base_12[1]
   + 0.5 alpha[4] * z_p_prime_base_12[1] + 0.5 p[0] * z_p_prime_base_12[1]
   - 0.5 p[1] * z_p_prime_base_12[1] + 0.5 p[2] * z_p_prime_base_12[1]
   - 0.5 p[3] * z_p_prime_base_12[1] + 0.5 p[4] * z_p_prime_base_12[1] ]
   <= 0
 infeasibility_p_prime_base_12_row_2_leq_0: 5e-06 z_p_prime_base_12[2] + [
   0.5 alpha[0] * z_p_prime_base_12[2]
   + 0.5 alpha[1] * z_p_prime_base_12[2]
   + 0.5 alpha[2] * z_p_prime_base_12[2]
   - 0.5 alpha[3] * z_p_prime_base_12[2]
   + 0.5 alpha[4] * z_p_prime_base_12[2] + 0.5 p[0] * z_p_prime_base_12[2]
   - 0.5 p[1] * z_p_prime_base_12[2] - 0.5 p[2] * z_p_prime_base_12[2]
   + 0.5 p[3] * z_p_prime_base_12[2] + 0.5 p[4] * z_p_prime_base_12[2] ]
   <= 0
 infeasibility_p_prime_base_12_row_3_leq_0: 5e-06 z_p_prime_base_12[3] + [
   alpha[4] * z_p_prime_base_12[3] + p[4] * z_p_prime_base_12[3] ] <= 0
 infeasibility_p_prime_base_12_row_4_leq_0: 5e-06 z_p_prime_base_12[4] + [
   - 0.5 alpha[0] * z_p_prime_base_12[4]
   - 0.5 alpha[1] * z_p_prime_base_12[4]
   - 0.5 alpha[2] * z_p_prime_base_12[4]
   - 0.5 alpha[3] * z_p_prime_base_12[4]
   - 0.5 alpha[4] * z_p_prime_base_12[4] - 0.5 p[0] * z_p_prime_base_12[4]
   + 0.5 p[1] * z_p_prime_base_12[4] + 0.5 p[2] * z_p_prime_base_12[4]
   + 0.5 p[3] * z_p_prime_base_12[4] - 0.5 p[4] * z_p_prime_base_12[4] ]
   <= 0
 infeasibility_p_base_13_row_0_leq_0: 5e-06 z_p_base_13[0] + [
   p[1] * z_p_base_13[0] ] <= 0
 infeasibility_p_base_13_row_1_leq_0: 5e-06 z_p_base_13[1] + [
   0.5 p[0] * z_p_base_13[1] - 0.5 p[1] * z_p_base_13[1]
   + 0.5 p[2] * z_p_base_13[1] - 0.5 p[3] * z_p_base_13[1]
   - 0.5 p[4] * z_p_base_13[1] ] <= 0
 infeasibility_p_base_13_row_2_leq_0: 5e-06 z_p_base_13[2] + [
   0.5 p[0] * z_p_base_13[2] - 0.5 p[1] * z_p_base_13[2]
   - 0.5 p[2] * z_p_base_13[2] + 0.5 p[3] * z_p_base_13[2]
   + 0.5 p[4] * z_p_base_13[2] ] <= 0
 infeasibility_p_base_13_row_3_leq_0: 5e-06 z_p_base_13[3] + [
   - 0.5 p[0] * z_p_base_13[3] + 0.5 p[1] * z_p_base_13[3]
   + 0.5 p[2] * z_p_base_13[3] + 0.5 p[3] * z_p_base_13[3]
   - 0.5 p[4] * z_p_base_13[3] ] <= 0
 infeasibility_p_base_13_row_4_leq_0: 5e-06 z_p_base_13[4] + [
   p[4] * z_p_base_13[4] ] <= 0
 infeasibility_p_prime_base_13_row_0_leq_0: 5e-06 z_p_prime_base_13[0] + [
   - alpha[1] * z_p_prime_base_13[0] + p[1] * z_p_prime_base_13[0] ] <= 0
 infeasibility_p_prime_base_13_row_1_leq_0: 5e-06 z_p_prime_base_13[1] + [
   0.5 alpha[0] * z_p_prime_base_13[1]
   + 0.5 alpha[1] * z_p_prime_base_13[1]
   - 0.5 alpha[2] * z_p_prime_base_13[1]
   + 0.5 alpha[3] * z_p_prime_base_13[1]
   - 0.5 alpha[4] * z_p_prime_base_13[1] + 0.5 p[0] * z_p_prime_base_13[1]
   - 0.5 p[1] * z_p_prime_base_13[1] + 0.5 p[2] * z_p_prime_base_13[1]
   - 0.5 p[3] * z_p_prime_base_13[1] - 0.5 p[4] * z_p_prime_base_13[1] ]
   <= 0
 infeasibility_p_prime_base_13_row_2_leq_0: 5e-06 z_p_prime_base_13[2] + [
   0.5 alpha[0] * z_p_prime_base_13[2]
   + 0.5 alpha[1] * z_p_prime_base_13[2]
   + 0.5 alpha[2] * z_p_prime_base_13[2]
   - 0.5 alpha[3] * z_p_prime_base_13[2]
   + 0.5 alpha[4] * z_p_prime_base_13[2] + 0.5 p[0] * z_p_prime_base_13[2]
   - 0.5 p[1] * z_p_prime_base_13[2] - 0.5 p[2] * z_p_prime_base_13[2]
   + 0.5 p[3] * z_p_prime_base_13[2] + 0.5 p[4] * z_p_prime_base_13[2] ]
   <= 0
 infeasibility_p_prime_base_13_row_3_leq_0: 5e-06 z_p_prime_base_13[3] + [
   - 0.5 alpha[0] * z_p_prime_base_13[3]
   - 0.5 alpha[1] * z_p_prime_base_13[3]
   - 0.5 alpha[2] * z_p_prime_base_13[3]
   - 0.5 alpha[3] * z_p_prime_base_13[3]
   - 0.5 alpha[4] * z_p_prime_base_13[3] - 0.5 p[0] * z_p_prime_base_13[3]
   + 0.5 p[1] * z_p_prime_base_13[3] + 0.5 p[2] * z_p_prime_base_13[3]
   + 0.5 p[3] * z_p_prime_base_13[3] - 0.5 p[4] * z_p_prime_base_13[3] ]
   <= 0
 infeasibility_p_prime_base_13_row_4_leq_0: 5e-06 z_p_prime_base_13[4] + [
   alpha[4] * z_p_prime_base_13[4] + p[4] * z_p_prime_base_13[4] ] <= 0
 infeasibility_p_base_14_row_0_leq_0: 5e-06 z_p_base_14[0] + [
   p[1] * z_p_base_14[0] ] <= 0
 infeasibility_p_base_14_row_1_leq_0: 5e-06 z_p_base_14[1] + [
   0.5 p[0] * z_p_base_14[1] - 0.5 p[1] * z_p_base_14[1]
   + 0.5 p[2] * z_p_base_14[1] - 0.5 p[3] * z_p_base_14[1]
   + 0.5 p[4] * z_p_base_14[1] ] <= 0
 infeasibility_p_base_14_row_2_leq_0: 5e-06 z_p_base_14[2] + [
   0.5 p[0] * z_p_base_14[2] - 0.5 p[1] * z_p_base_14[2]
   - 0.5 p[2] * z_p_base_14[2] + 0.5 p[3] * z_p_base_14[2]
   - 0.5 p[4] * z_p_base_14[2] ] <= 0
 infeasibility_p_base_14_row_3_leq_0: 5e-06 z_p_base_14[3] + [
   - 0.5 p[0] * z_p_base_14[3] + 0.5 p[1] * z_p_base_14[3]
   + 0.5 p[2] * z_p_base_14[3] + 0.5 p[3] * z_p_base_14[3]
   - 0.5 p[4] * z_p_base_14[3] ] <= 0
 infeasibility_p_base_14_row_4_leq_0: 5e-06 z_p_base_14[4] + [
   p[4] * z_p_base_14[4] ] <= 0
 infeasibility_p_prime_base_14_row_0_leq_0: 5e-06 z_p_prime_base_14[0] + [
   - alpha[1] * z_p_prime_base_14[0] + p[1] * z_p_prime_base_14[0] ] <= 0
 infeasibility_p_prime_base_14_row_1_leq_0: 5e-06 z_p_prime_base_14[1] + [
   0.5 alpha[0] * z_p_prime_base_14[1]
   + 0.5 alpha[1] * z_p_prime_base_14[1]
   - 0.5 alpha[2] * z_p_prime_base_14[1]
   + 0.5 alpha[3] * z_p_prime_base_14[1]
   + 0.5 alpha[4] * z_p_prime_base_14[1] + 0.5 p[0] * z_p_prime_base_14[1]
   - 0.5 p[1] * z_p_prime_base_14[1] + 0.5 p[2] * z_p_prime_base_14[1]
   - 0.5 p[3] * z_p_prime_base_14[1] + 0.5 p[4] * z_p_prime_base_14[1] ]
   <= 0
 infeasibility_p_prime_base_14_row_2_leq_0: 5e-06 z_p_prime_base_14[2] + [
   0.5 alpha[0] * z_p_prime_base_14[2]
   + 0.5 alpha[1] * z_p_prime_base_14[2]
   + 0.5 alpha[2] * z_p_prime_base_14[2]
   - 0.5 alpha[3] * z_p_prime_base_14[2]
   - 0.5 alpha[4] * z_p_prime_base_14[2] + 0.5 p[0] * z_p_prime_base_14[2]
   - 0.5 p[1] * z_p_prime_base_14[2] - 0.5 p[2] * z_p_prime_base_14[2]
   + 0.5 p[3] * z_p_prime_base_14[2] - 0.5 p[4] * z_p_prime_base_14[2] ]
   <= 0
 infeasibility_p_prime_base_14_row_3_leq_0: 5e-06 z_p_prime_base_14[3] + [
   - 0.5 alpha[0] * z_p_prime_base_14[3]
   - 0.5 alpha[1] * z_p_prime_base_14[3]
   - 0.5 alpha[2] * z_p_prime_base_14[3]
   - 0.5 alpha[3] * z_p_prime_base_14[3]
   - 0.5 alpha[4] * z_p_prime_base_14[3] - 0.5 p[0] * z_p_prime_base_14[3]
   + 0.5 p[1] * z_p_prime_base_14[3] + 0.5 p[2] * z_p_prime_base_14[3]
   + 0.5 p[3] * z_p_prime_base_14[3] - 0.5 p[4] * z_p_prime_base_14[3] ]
   <= 0
 infeasibility_p_prime_base_14_row_4_leq_0: 5e-06 z_p_prime_base_14[4] + [
   alpha[4] * z_p_prime_base_14[4] + p[4] * z_p_prime_base_14[4] ] <= 0
 infeasibility_p_base_15_row_0_leq_0: 5e-06 z_p_base_15[0] + [
   0.5 p[0] * z_p_base_15[0] + 0.5 p[1] * z_p_base_15[0]
   - 0.5 p[2] * z_p_base_15[0] - 0.5 p[3] * z_p_base_15[0]
   - 0.5 p[4] * z_p_base_15[0] ] <= 0
 infeasibility_p_base_15_row_1_leq_0: 5e-06 z_p_base_15[1] + [
   0.5 p[0] * z_p_base_15[1] - 0.5 p[1] * z_p_base_15[1]
   + 0.5 p[2] * z_p_base_15[1] + 0.5 p[3] * z_p_base_15[1]
   - 0.5 p[4] * z_p_base_15[1] ] <= 0
 infeasibility_p_base_15_row_2_leq_0: 5e-06 z_p_base_15[2] + [
   p[4] * z_p_base_15[2] ] <= 0
 infeasibility_p_base_15_row_3_leq_0: 5e-06 z_p_base_15[3] + [
   - 0.5 p[0] * z_p_base_15[3] + 0.5 p[1] * z_p_base_15[3]
   + 0.5 p[2] * z_p_base_15[3] - 0.5 p[3] * z_p_base_15[3]
   + 0.5 p[4] * z_p_base_15[3] ] <= 0
 infeasibility_p_base_15_row_4_leq_0: 5e-06 z_p_base_15[4] + [
   p[3] * z_p_base_15[4] ] <= 0
 infeasibility_p_prime_base_15_row_0_leq_0: 5e-06 z_p_prime_base_15[0] + [
   0.5 alpha[0] * z_p_prime_base_15[0]
   - 0.5 alpha[1] * z_p_prime_base_15[0]
   + 0.5 alpha[2] * z_p_prime_base_15[0]
   + 0.5 alpha[3] * z_p_prime_base_15[0]
   - 0.5 alpha[4] * z_p_prime_base_15[0] + 0.5 p[0] * z_p_prime_base_15[0]
   + 0.5 p[1] * z_p_prime_base_15[0] - 0.5 p[2] * z_p_prime_base_15[0]
   - 0.5 p[3] * z_p_prime_base_15[0] - 0.5 p[4] * z_p_prime_base_15[0] ]
   <= 0
 infeasibility_p_prime_base_15_row_1_leq_0: 5e-06 z_p_prime_base_15[1] + [
   0.5 alpha[0] * z_p_prime_base_15[1]
   + 0.5 alpha[1] * z_p_prime_base_15[1]
   - 0.5 alpha[2] * z_p_prime_base_15[1]
   - 0.5 alpha[3] * z_p_prime_base_15[1]
   - 0.5 alpha[4] * z_p_prime_base_15[1] + 0.5 p[0] * z_p_prime_base_15[1]
   - 0.5 p[1] * z_p_prime_base_15[1] + 0.5 p[2] * z_p_prime_base_15[1]
   + 0.5 p[3] * z_p_prime_base_15[1] - 0.5 p[4] * z_p_prime_base_15[1] ]
   <= 0
 infeasibility_p_prime_base_15_row_2_leq_0: 5e-06 z_p_prime_base_15[2] + [
   alpha[4] * z_p_prime_base_15[2] + p[4] * z_p_prime_base_15[2] ] <= 0
 infeasibility_p_prime_base_15_row_3_leq_0: 5e-06 z_p_prime_base_15[3] + [
   - 0.5 alpha[0] * z_p_prime_base_15[3]
   - 0.5 alpha[1] * z_p_prime_base_15[3]
   - 0.5 alpha[2] * z_p_prime_base_15[3]
   + 0.5 alpha[3] * z_p_prime_base_15[3]
   + 0.5 alpha[4] * z_p_prime_base_15[3] - 0.5 p[0] * z_p_prime_base_15[3]
   + 0.5 p[1] * z_p_prime_base_15[3] + 0.5 p[2] * z_p_prime_base_15[3]
   - 0.5 p[3] * z_p_prime_base_15[3] + 0.5 p[4] * z_p_prime_base_15[3] ]
   <= 0
 infeasibility_p_prime_base_15_row_4_leq_0: 5e-06 z_p_prime_base_15[4] + [
   - alpha[3] * z_p_prime_base_15[4] + p[3] * z_p_prime_base_15[4] ] <= 0
 infeasibility_p_base_16_row_0_leq_0: 5e-06 z_p_base_16[0] + [
   0.5 p[0] * z_p_base_16[0] + 0.5 p[1] * z_p_base_16[0]
   - 0.5 p[2] * z_p_base_16[0] + 0.5 p[3] * z_p_base_16[0]
   - 0.5 p[4] * z_p_base_16[0] ] <= 0
 infeasibility_p_base_16_row_1_leq_0: 5e-06 z_p_base_16[1] + [
   0.5 p[0] * z_p_base_16[1] - 0.5 p[1] * z_p_base_16[1]
   + 0.5 p[2] * z_p_base_16[1] - 0.5 p[3] * z_p_base_16[1]
   - 0.5 p[4] * z_p_base_16[1] ] <= 0
 infeasibility_p_base_16_row_2_leq_0: 5e-06 z_p_base_16[2] + [
   p[4] * z_p_base_16[2] ] <= 0
 infeasibility_p_base_16_row_3_leq_0: 5e-06 z_p_base_16[3] + [
   - 0.5 p[0] * z_p_base_16[3] + 0.5 p[1] * z_p_base_16[3]
   + 0.5 p[2] * z_p_base_16[3] - 0.5 p[3] * z_p_base_16[3]
   + 0.5 p[4] * z_p_base_16[3] ] <= 0
 infeasibility_p_base_16_row_4_leq_0: 5e-06 z_p_base_16[4] + [
   p[3] * z_p_base_16[4] ] <= 0
 infeasibility_p_prime_base_16_row_0_leq_0: 5e-06 z_p_prime_base_16[0] + [
   0.5 alpha[0] * z_p_prime_base_16[0]
   - 0.5 alpha[1] * z_p_prime_base_16[0]
   + 0.5 alpha[2] * z_p_prime_base_16[0]
   - 0.5 alpha[3] * z_p_prime_base_16[0]
   - 0.5 alpha[4] * z_p_prime_base_16[0] + 0.5 p[0] * z_p_prime_base_16[0]
   + 0.5 p[1] * z_p_prime_base_16[0] - 0.5 p[2] * z_p_prime_base_16[0]
   + 0.5 p[3] * z_p_prime_base_16[0] - 0.5 p[4] * z_p_prime_base_16[0] ]
   <= 0
 infeasibility_p_prime_base_16_row_1_leq_0: 5e-06 z_p_prime_base_16[1] + [
   0.5 alpha[0] * z_p_prime_base_16[1]
   + 0.5 alpha[1] * z_p_prime_base_16[1]
   - 0.5 alpha[2] * z_p_prime_base_16[1]
   + 0.5 alpha[3] * z_p_prime_base_16[1]
   - 0.5 alpha[4] * z_p_prime_base_16[1] + 0.5 p[0] * z_p_prime_base_16[1]
   - 0.5 p[1] * z_p_prime_base_16[1] + 0.5 p[2] * z_p_prime_base_16[1]
   - 0.5 p[3] * z_p_prime_base_16[1] - 0.5 p[4] * z_p_prime_base_16[1] ]
   <= 0
 infeasibility_p_prime_base_16_row_2_leq_0: 5e-06 z_p_prime_base_16[2] + [
   alpha[4] * z_p_prime_base_16[2] + p[4] * z_p_prime_base_16[2] ] <= 0
 infeasibility_p_prime_base_16_row_3_leq_0: 5e-06 z_p_prime_base_16[3] + [
   - 0.5 alpha[0] * z_p_prime_base_16[3]
   - 0.5 alpha[1] * z_p_prime_base_16[3]
   - 0.5 alpha[2] * z_p_prime_base_16[3]
   + 0.5 alpha[3] * z_p_prime_base_16[3]
   + 0.5 alpha[4] * z_p_prime_base_16[3] - 0.5 p[0] * z_p_prime_base_16[3]
   + 0.5 p[1] * z_p_prime_base_16[3] + 0.5 p[2] * z_p_prime_base_16[3]
   - 0.5 p[3] * z_p_prime_base_16[3] + 0.5 p[4] * z_p_prime_base_16[3] ]
   <= 0
 infeasibility_p_prime_base_16_row_4_leq_0: 5e-06 z_p_prime_base_16[4] + [
   - alpha[3] * z_p_prime_base_16[4] + p[3] * z_p_prime_base_16[4] ] <= 0
 infeasibility_p_base_17_row_0_leq_0: 5e-06 z_p_base_17[0] + [
   0.5 p[0] * z_p_base_17[0] + 0.5 p[1] * z_p_base_17[0]
   - 0.5 p[2] * z_p_base_17[0] + 0.5 p[3] * z_p_base_17[0]
   - 0.5 p[4] * z_p_base_17[0] ] <= 0
 infeasibility_p_base_17_row_1_leq_0: 5e-06 z_p_base_17[1] + [
   0.5 p[0] * z_p_base_17[1] - 0.5 p[1] * z_p_base_17[1]
   + 0.5 p[2] * z_p_base_17[1] + 0.5 p[3] * z_p_base_17[1]
   - 0.5 p[4] * z_p_base_17[1] ] <= 0
 infeasibility_p_base_17_row_2_leq_0: 5e-06 z_p_base_17[2] + [
   - p[3] * z_p_base_17[2] + p[4] * z_p_base_17[2] ] <= 0
 infeasibility_p_base_17_row_3_leq_0: 5e-06 z_p_base_17[3] + [
   - 0.5 p[0] * z_p_base_17[3] + 0.5 p[1] * z_p_base_17[3]
   + 0.5 p[2] * z_p_base_17[3] - 0.5 p[3] * z_p_base_17[3]
   + 0.5 p[4] * z_p_base_17[3] ] <= 0
 infeasibility_p_base_17_row_4_leq_0: 5e-06 z_p_base_17[4] + [
   p[3] * z_p_base_17[4] ] <= 0
 infeasibility_p_prime_base_17_row_0_leq_0: 5e-06 z_p_prime_base_17[0] + [
   0.5 alpha[0] * z_p_prime_base_17[0]
   - 0.5 alpha[1] * z_p_prime_base_17[0]
   + 0.5 alpha[2] * z_p_prime_base_17[0]
   - 0.5 alpha[3] * z_p_prime_base_17[0]
   - 0.5 alpha[4] * z_p_prime_base_17[0] + 0.5 p[0] * z_p_prime_base_17[0]
   + 0.5 p[1] * z_p_prime_base_17[0] - 0.5 p[2] * z_p_prime_base_17[0]
   + 0.5 p[3] * z_p_prime_base_17[0] - 0.5 p[4] * z_p_prime_base_17[0] ]
   <= 0
 infeasibility_p_prime_base_17_row_1_leq_0: 5e-06 z_p_prime_base_17[1] + [
   0.5 alpha[0] * z_p_prime_base_17[1]
   + 0.5 alpha[1] * z_p_prime_base_17[1]
   - 0.5 alpha[2] * z_p_prime_base_17[1]
   - 0.5 alpha[3] * z_p_prime_base_17[1]
   - 0.5 alpha[4] * z_p_prime_base_17[1] + 0.5 p[0] * z_p_prime_base_17[1]
   - 0.5 p[1] * z_p_prime_base_17[1] + 0.5 p[2] * z_p_prime_base_17[1]
   + 0.5 p[3] * z_p_prime_base_17[1] - 0.5 p[4] * z_p_prime_base_17[1] ]
   <= 0
 infeasibility_p_prime_base_17_row_2_leq_0: 5e-06 z_p_prime_base_17[2] + [
   alpha[3] * z_p_prime_base_17[2] + alpha[4] * z_p_prime_base_17[2]
   - p[3] * z_p_prime_base_17[2] + p[4] * z_p_prime_base_17[2] ] <= 0
 infeasibility_p_prime_base_17_row_3_leq_0: 5e-06 z_p_prime_base_17[3] + [
   - 0.5 alpha[0] * z_p_prime_base_17[3]
   - 0.5 alpha[1] * z_p_prime_base_17[3]
   - 0.5 alpha[2] * z_p_prime_base_17[3]
   + 0.5 alpha[3] * z_p_prime_base_17[3]
   + 0.5 alpha[4] * z_p_prime_base_17[3] - 0.5 p[0] * z_p_prime_base_17[3]
   + 0.5 p[1] * z_p_prime_base_17[3] + 0.5 p[2] * z_p_prime_base_17[3]
   - 0.5 p[3] * z_p_prime_base_17[3] + 0.5 p[4] * z_p_prime_base_17[3] ]
   <= 0
 infeasibility_p_prime_base_17_row_4_leq_0: 5e-06 z_p_prime_base_17[4] + [
   - alpha[3] * z_p_prime_base_17[4] + p[3] * z_p_prime_base_17[4] ] <= 0
 infeasibility_p_base_18_row_0_leq_0: 5e-06 z_p_base_18[0] + [
   0.5 p[0] * z_p_base_18[0] + 0.5 p[1] * z_p_base_18[0]
   - 0.5 p[2] * z_p_base_18[0] - 0.5 p[3] * z_p_base_18[0]
   - 0.5 p[4] * z_p_base_18[0] ] <= 0
 infeasibility_p_base_18_row_1_leq_0: 5e-06 z_p_base_18[1] + [
   p[2] * z_p_base_18[1] ] <= 0
 infeasibility_p_base_18_row_2_leq_0: 5e-06 z_p_base_18[2] + [
   0.5 p[0] * z_p_base_18[2] - 0.5 p[1] * z_p_base_18[2]
   - 0.5 p[2] * z_p_base_18[2] + 0.5 p[3] * z_p_base_18[2]
   + 0.5 p[4] * z_p_base_18[2] ] <= 0
 infeasibility_p_base_18_row_3_leq_0: 5e-06 z_p_base_18[3] + [
   p[3] * z_p_base_18[3] ] <= 0
 infeasibility_p_base_18_row_4_leq_0: 5e-06 z_p_base_18[4] + [
   - 0.5 p[0] * z_p_base_18[4] + 0.5 p[1] * z_p_base_18[4]
   + 0.5 p[2] * z_p_base_18[4] - 0.5 p[3] * z_p_base_18[4]
   + 0.5 p[4] * z_p_base_18[4] ] <= 0
 infeasibility_p_prime_base_18_row_0_leq_0: 5e-06 z_p_prime_base_18[0] + [
   0.5 alpha[0] * z_p_prime_base_18[0]
   - 0.5 alpha[1] * z_p_prime_base_18[0]
   + 0.5 alpha[2] * z_p_prime_base_18[0]
   + 0.5 alpha[3] * z_p_prime_base_18[0]
   - 0.5 alpha[4] * z_p_prime_base_18[0] + 0.5 p[0] * z_p_prime_base_18[0]
   + 0.5 p[1] * z_p_prime_base_18[0] - 0.5 p[2] * z_p_prime_base_18[0]
   - 0.5 p[3] * z_p_prime_base_18[0] - 0.5 p[4] * z_p_prime_base_18[0] ]
   <= 0
 infeasibility_p_prime_base_18_row_1_leq_0: 5e-06 z_p_prime_base_18[1] + [
   - alpha[2] * z_p_prime_base_18[1] + p[2] * z_p_prime_base_18[1] ] <= 0
 infeasibility_p_prime_base_18_row_2_leq_0: 5e-06 z_p_prime_base_18[2] + [
   0.5 alpha[0] * z_p_prime_base_18[2]
   + 0.5 alpha[1] * z_p_prime_base_18[2]
   + 0.5 alpha[2] * z_p_prime_base_18[2]
   - 0.5 alpha[3] * z_p_prime_base_18[2]
   + 0.5 alpha[4] * z_p_prime_base_18[2] + 0.5 p[0] * z_p_prime_base_18[2]
   - 0.5 p[1] * z_p_prime_base_18[2] - 0.5 p[2] * z_p_prime_base_18[2]
   + 0.5 p[3] * z_p_prime_base_18[2] + 0.5 p[4] * z_p_prime_base_18[2] ]
   <= 0
 infeasibility_p_prime_base_18_row_3_leq_0: 5e-06 z_p_prime_base_18[3] + [
   - alpha[3] * z_p_prime_base_18[3] + p[3] * z_p_prime_base_18[3] ] <= 0
 infeasibility_p_prime_base_18_row_4_leq_0: 5e-06 z_p_prime_base_18[4] + [
   - 0.5 alpha[0] * z_p_prime_base_18[4]
   - 0.5 alpha[1] * z_p_prime_base_18[4]
   - 0.5 alpha[2] * z_p_prime_base_18[4]
   + 0.5 alpha[3] * z_p_prime_base_18[4]
   + 0.5 alpha[4] * z_p_prime_base_18[4] - 0.5 p[0] * z_p_prime_base_18[4]
   + 0.5 p[1] * z_p_prime_base_18[4] + 0.5 p[2] * z_p_prime_base_18[4]
   - 0.5 p[3] * z_p_prime_base_18[4] + 0.5 p[4] * z_p_prime_base_18[4] ]
   <= 0
 infeasibility_p_base_19_row_0_leq_0: 5e-06 z_p_base_19[0] + [
   p[1] * z_p_base_19[0] - p[3] * z_p_base_19[0] ] <= 0
 infeasibility_p_base_19_row_1_leq_0: 5e-06 z_p_base_19[1] + [
   0.5 p[0] * z_p_base_19[1] - 0.5 p[1] * z_p_base_19[1]
   + 0.5 p[2] * z_p_base_19[1] + 0.5 p[3] * z_p_base_19[1]
   - 0.5 p[4] * z_p_base_19[1] ] <= 0
 infeasibility_p_base_19_row_2_leq_0: 5e-06 z_p_base_19[2] + [
   0.5 p[0] * z_p_base_19[2] - 0.5 p[1] * z_p_base_19[2]
   - 0.5 p[2] * z_p_base_19[2] + 0.5 p[3] * z_p_base_19[2]
   + 0.5 p[4] * z_p_base_19[2] ] <= 0
 infeasibility_p_base_19_row_3_leq_0: 5e-06 z_p_base_19[3] + [
   p[3] * z_p_base_19[3] ] <= 0
 infeasibility_p_base_19_row_4_leq_0: 5e-06 z_p_base_19[4] + [
   - 0.5 p[0] * z_p_base_19[4] + 0.5 p[1] * z_p_base_19[4]
   + 0.5 p[2] * z_p_base_19[4] - 0.5 p[3] * z_p_base_19[4]
   + 0.5 p[4] * z_p_base_19[4] ] <= 0
 infeasibility_p_prime_base_19_row_0_leq_0: 5e-06 z_p_prime_base_19[0] + [
   - alpha[1] * z_p_prime_base_19[0] + alpha[3] * z_p_prime_base_19[0]
   + p[1] * z_p_prime_base_19[0] - p[3] * z_p_prime_base_19[0] ] <= 0
 infeasibility_p_prime_base_19_row_1_leq_0: 5e-06 z_p_prime_base_19[1] + [
   0.5 alpha[0] * z_p_prime_base_19[1]
   + 0.5 alpha[1] * z_p_prime_base_19[1]
   - 0.5 alpha[2] * z_p_prime_base_19[1]
   - 0.5 alpha[3] * z_p_prime_base_19[1]
   - 0.5 alpha[4] * z_p_prime_base_19[1] + 0.5 p[0] * z_p_prime_base_19[1]
   - 0.5 p[1] * z_p_prime_base_19[1] + 0.5 p[2] * z_p_prime_base_19[1]
   + 0.5 p[3] * z_p_prime_base_19[1] - 0.5 p[4] * z_p_prime_base_19[1] ]
   <= 0
 infeasibility_p_prime_base_19_row_2_leq_0: 5e-06 z_p_prime_base_19[2] + [
   0.5 alpha[0] * z_p_prime_base_19[2]
   + 0.5 alpha[1] * z_p_prime_base_19[2]
   + 0.5 alpha[2] * z_p_prime_base_19[2]
   - 0.5 alpha[3] * z_p_prime_base_19[2]
   + 0.5 alpha[4] * z_p_prime_base_19[2] + 0.5 p[0] * z_p_prime_base_19[2]
   - 0.5 p[1] * z_p_prime_base_19[2] - 0.5 p[2] * z_p_prime_base_19[2]
   + 0.5 p[3] * z_p_prime_base_19[2] + 0.5 p[4] * z_p_prime_base_19[2] ]
   <= 0
 infeasibility_p_prime_base_19_row_3_leq_0: 5e-06 z_p_prime_base_19[3] + [
   - alpha[3] * z_p_prime_base_19[3] + p[3] * z_p_prime_base_19[3] ] <= 0
 infeasibility_p_prime_base_19_row_4_leq_0: 5e-06 z_p_prime_base_19[4] + [
   - 0.5 alpha[0] * z_p_prime_base_19[4]
   - 0.5 alpha[1] * z_p_prime_base_19[4]
   - 0.5 alpha[2] * z_p_prime_base_19[4]
   + 0.5 alpha[3] * z_p_prime_base_19[4]
   + 0.5 alpha[4] * z_p_prime_base_19[4] - 0.5 p[0] * z_p_prime_base_19[4]
   + 0.5 p[1] * z_p_prime_base_19[4] + 0.5 p[2] * z_p_prime_base_19[4]
   - 0.5 p[3] * z_p_prime_base_19[4] + 0.5 p[4] * z_p_prime_base_19[4] ]
   <= 0
 infeasibility_p_base_20_row_0_leq_0: 5e-06 z_p_base_20[0] + [
   0.5 p[0] * z_p_base_20[0] + 0.5 p[1] * z_p_base_20[0]
   - 0.5 p[2] * z_p_base_20[0] + 0.5 p[3] * z_p_base_20[0]
   - 0.5 p[4] * z_p_base_20[0] ] <= 0
 infeasibility_p_base_20_row_1_leq_0: 5e-06 z_p_base_20[1] + [
   p[2] * z_p_base_20[1] - p[3] * z_p_base_20[1] ] <= 0
 infeasibility_p_base_20_row_2_leq_0: 5e-06 z_p_base_20[2] + [
   0.5 p[0] * z_p_base_20[2] - 0.5 p[1] * z_p_base_20[2]
   - 0.5 p[2] * z_p_base_20[2] + 0.5 p[3] * z_p_base_20[2]
   + 0.5 p[4] * z_p_base_20[2] ] <= 0
 infeasibility_p_base_20_row_3_leq_0: 5e-06 z_p_base_20[3] + [
   - 0.5 p[0] * z_p_base_20[3] + 0.5 p[1] * z_p_base_20[3]
   + 0.5 p[2] * z_p_base_20[3] - 0.5 p[3] * z_p_base_20[3]
   + 0.5 p[4] * z_p_base_20[3] ] <= 0
 infeasibility_p_base_20_row_4_leq_0: 5e-06 z_p_base_20[4] + [
   p[3] * z_p_base_20[4] ] <= 0
 infeasibility_p_prime_base_20_row_0_leq_0: 5e-06 z_p_prime_base_20[0] + [
   0.5 alpha[0] * z_p_prime_base_20[0]
   - 0.5 alpha[1] * z_p_prime_base_20[0]
   + 0.5 alpha[2] * z_p_prime_base_20[0]
   - 0.5 alpha[3] * z_p_prime_base_20[0]
   - 0.5 alpha[4] * z_p_prime_base_20[0] + 0.5 p[0] * z_p_prime_base_20[0]
   + 0.5 p[1] * z_p_prime_base_20[0] - 0.5 p[2] * z_p_prime_base_20[0]
   + 0.5 p[3] * z_p_prime_base_20[0] - 0.5 p[4] * z_p_prime_base_20[0] ]
   <= 0
 infeasibility_p_prime_base_20_row_1_leq_0: 5e-06 z_p_prime_base_20[1] + [
   - alpha[2] * z_p_prime_base_20[1] + alpha[3] * z_p_prime_base_20[1]
   + p[2] * z_p_prime_base_20[1] - p[3] * z_p_prime_base_20[1] ] <= 0
 infeasibility_p_prime_base_20_row_2_leq_0: 5e-06 z_p_prime_base_20[2] + [
   0.5 alpha[0] * z_p_prime_base_20[2]
   + 0.5 alpha[1] * z_p_prime_base_20[2]
   + 0.5 alpha[2] * z_p_prime_base_20[2]
   - 0.5 alpha[3] * z_p_prime_base_20[2]
   + 0.5 alpha[4] * z_p_prime_base_20[2] + 0.5 p[0] * z_p_prime_base_20[2]
   - 0.5 p[1] * z_p_prime_base_20[2] - 0.5 p[2] * z_p_prime_base_20[2]
   + 0.5 p[3] * z_p_prime_base_20[2] + 0.5 p[4] * z_p_prime_base_20[2] ]
   <= 0
 infeasibility_p_prime_base_20_row_3_leq_0: 5e-06 z_p_prime_base_20[3] + [
   - 0.5 alpha[0] * z_p_prime_base_20[3]
   - 0.5 alpha[1] * z_p_prime_base_20[3]
   - 0.5 alpha[2] * z_p_prime_base_20[3]
   + 0.5 alpha[3] * z_p_prime_base_20[3]
   + 0.5 alpha[4] * z_p_prime_base_20[3] - 0.5 p[0] * z_p_prime_base_20[3]
   + 0.5 p[1] * z_p_prime_base_20[3] + 0.5 p[2] * z_p_prime_base_20[3]
   - 0.5 p[3] * z_p_prime_base_20[3] + 0.5 p[4] * z_p_prime_base_20[3] ]
   <= 0
 infeasibility_p_prime_base_20_row_4_leq_0: 5e-06 z_p_prime_base_20[4] + [
   - alpha[3] * z_p_prime_base_20[4] + p[3] * z_p_prime_base_20[4] ] <= 0
 infeasibility_p_base_21_row_0_leq_0: 5e-06 z_p_base_21[0] + [
   0.5 p[0] * z_p_base_21[0] + 0.5 p[1] * z_p_base_21[0]
   - 0.5 p[2] * z_p_base_21[0] + 0.5 p[3] * z_p_base_21[0]
   - 0.5 p[4] * z_p_base_21[0] ] <= 0
 infeasibility_p_base_21_row_1_leq_0: 5e-06 z_p_base_21[1] + [
   p[2] * z_p_base_21[1] ] <= 0
 infeasibility_p_base_21_row_2_leq_0: 5e-06 z_p_base_21[2] + [
   0.5 p[0] * z_p_base_21[2] - 0.5 p[1] * z_p_base_21[2]
   - 0.5 p[2] * z_p_base_21[2] - 0.5 p[3] * z_p_base_21[2]
   + 0.5 p[4] * z_p_base_21[2] ] <= 0
 infeasibility_p_base_21_row_3_leq_0: 5e-06 z_p_base_21[3] + [
   - 0.5 p[0] * z_p_base_21[3] + 0.5 p[1] * z_p_base_21[3]
   + 0.5 p[2] * z_p_base_21[3] - 0.5 p[3] * z_p_base_21[3]
   + 0.5 p[4] * z_p_base_21[3] ] <= 0
 infeasibility_p_base_21_row_4_leq_0: 5e-06 z_p_base_21[4] + [
   p[3] * z_p_base_21[4] ] <= 0
 infeasibility_p_prime_base_21_row_0_leq_0: 5e-06 z_p_prime_base_21[0] + [
   0.5 alpha[0] * z_p_prime_base_21[0]
   - 0.5 alpha[1] * z_p_prime_base_21[0]
   + 0.5 alpha[2] * z_p_prime_base_21[0]
   - 0.5 alpha[3] * z_p_prime_base_21[0]
   - 0.5 alpha[4] * z_p_prime_base_21[0] + 0.5 p[0] * z_p_prime_base_21[0]
   + 0.5 p[1] * z_p_prime_base_21[0] - 0.5 p[2] * z_p_prime_base_21[0]
   + 0.5 p[3] * z_p_prime_base_21[0] - 0.5 p[4] * z_p_prime_base_21[0] ]
   <= 0
 infeasibility_p_prime_base_21_row_1_leq_0: 5e-06 z_p_prime_base_21[1] + [
   - alpha[2] * z_p_prime_base_21[1] + p[2] * z_p_prime_base_21[1] ] <= 0
 infeasibility_p_prime_base_21_row_2_leq_0: 5e-06 z_p_prime_base_21[2] + [
   0.5 alpha[0] * z_p_prime_base_21[2]
   + 0.5 alpha[1] * z_p_prime_base_21[2]
   + 0.5 alpha[2] * z_p_prime_base_21[2]
   + 0.5 alpha[3] * z_p_prime_base_21[2]
   + 0.5 alpha[4] * z_p_prime_base_21[2] + 0.5 p[0] * z_p_prime_base_21[2]
   - 0.5 p[1] * z_p_prime_base_21[2] - 0.5 p[2] * z_p_prime_base_21[2]
   - 0.5 p[3] * z_p_prime_base_21[2] + 0.5 p[4] * z_p_prime_base_21[2] ]
   <= 0
 infeasibility_p_prime_base_21_row_3_leq_0: 5e-06 z_p_prime_base_21[3] + [
   - 0.5 alpha[0] * z_p_prime_base_21[3]
   - 0.5 alpha[1] * z_p_prime_base_21[3]
   - 0.5 alpha[2] * z_p_prime_base_21[3]
   + 0.5 alpha[3] * z_p_prime_base_21[3]
   + 0.5 alpha[4] * z_p_prime_base_21[3] - 0.5 p[0] * z_p_prime_base_21[3]
   + 0.5 p[1] * z_p_prime_base_21[3] + 0.5 p[2] * z_p_prime_base_21[3]
   - 0.5 p[3] * z_p_prime_base_21[3] + 0.5 p[4] * z_p_prime_base_21[3] ]
   <= 0
 infeasibility_p_prime_base_21_row_4_leq_0: 5e-06 z_p_prime_base_21[4] + [
   - alpha[3] * z_p_prime_base_21[4] + p[3] * z_p_prime_base_21[4] ] <= 0
 infeasibility_p_base_22_row_0_leq_0: 5e-06 z_p_base_22[0] + [
   p[1] * z_p_base_22[0] ] <= 0
 infeasibility_p_base_22_row_1_leq_0: 5e-06 z_p_base_22[1] + [
   0.5 p[0] * z_p_base_22[1] - 0.5 p[1] * z_p_base_22[1]
   + 0.5 p[2] * z_p_base_22[1] - 0.5 p[3] * z_p_base_22[1]
   - 0.5 p[4] * z_p_base_22[1] ] <= 0
 infeasibility_p_base_22_row_2_leq_0: 5e-06 z_p_base_22[2] + [
   0.5 p[0] * z_p_base_22[2] - 0.5 p[1] * z_p_base_22[2]
   - 0.5 p[2] * z_p_base_22[2] + 0.5 p[3] * z_p_base_22[2]
   + 0.5 p[4] * z_p_base_22[2] ] <= 0
 infeasibility_p_base_22_row_3_leq_0: 5e-06 z_p_base_22[3] + [
   p[3] * z_p_base_22[3] ] <= 0
 infeasibility_p_base_22_row_4_leq_0: 5e-06 z_p_base_22[4] + [
   - 0.5 p[0] * z_p_base_22[4] + 0.5 p[1] * z_p_base_22[4]
   + 0.5 p[2] * z_p_base_22[4] - 0.5 p[3] * z_p_base_22[4]
   + 0.5 p[4] * z_p_base_22[4] ] <= 0
 infeasibility_p_prime_base_22_row_0_leq_0: 5e-06 z_p_prime_base_22[0] + [
   - alpha[1] * z_p_prime_base_22[0] + p[1] * z_p_prime_base_22[0] ] <= 0
 infeasibility_p_prime_base_22_row_1_leq_0: 5e-06 z_p_prime_base_22[1] + [
   0.5 alpha[0] * z_p_prime_base_22[1]
   + 0.5 alpha[1] * z_p_prime_base_22[1]
   - 0.5 alpha[2] * z_p_prime_base_22[1]
   + 0.5 alpha[3] * z_p_prime_base_22[1]
   - 0.5 alpha[4] * z_p_prime_base_22[1] + 0.5 p[0] * z_p_prime_base_22[1]
   - 0.5 p[1] * z_p_prime_base_22[1] + 0.5 p[2] * z_p_prime_base_22[1]
   - 0.5 p[3] * z_p_prime_base_22[1] - 0.5 p[4] * z_p_prime_base_22[1] ]
   <= 0
 infeasibility_p_prime_base_22_row_2_leq_0: 5e-06 z_p_prime_base_22[2] + [
   0.5 alpha[0] * z_p_prime_base_22[2]
   + 0.5 alpha[1] * z_p_prime_base_22[2]
   + 0.5 alpha[2] * z_p_prime_base_22[2]
   - 0.5 alpha[3] * z_p_prime_base_22[2]
   + 0.5 alpha[4] * z_p_prime_base_22[2] + 0.5 p[0] * z_p_prime_base_22[2]
   - 0.5 p[1] * z_p_prime_base_22[2] - 0.5 p[2] * z_p_prime_base_22[2]
   + 0.5 p[3] * z_p_prime_base_22[2] + 0.5 p[4] * z_p_prime_base_22[2] ]
   <= 0
 infeasibility_p_prime_base_22_row_3_leq_0: 5e-06 z_p_prime_base_22[3] + [
   - alpha[3] * z_p_prime_base_22[3] + p[3] * z_p_prime_base_22[3] ] <= 0
 infeasibility_p_prime_base_22_row_4_leq_0: 5e-06 z_p_prime_base_22[4] + [
   - 0.5 alpha[0] * z_p_prime_base_22[4]
   - 0.5 alpha[1] * z_p_prime_base_22[4]
   - 0.5 alpha[2] * z_p_prime_base_22[4]
   + 0.5 alpha[3] * z_p_prime_base_22[4]
   + 0.5 alpha[4] * z_p_prime_base_22[4] - 0.5 p[0] * z_p_prime_base_22[4]
   + 0.5 p[1] * z_p_prime_base_22[4] + 0.5 p[2] * z_p_prime_base_22[4]
   - 0.5 p[3] * z_p_prime_base_22[4] + 0.5 p[4] * z_p_prime_base_22[4] ]
   <= 0
 infeasibility_p_base_23_row_0_leq_0: 5e-06 z_p_base_23[0] + [
   p[1] * z_p_base_23[0] ] <= 0
 infeasibility_p_base_23_row_1_leq_0: 5e-06 z_p_base_23[1] + [
   0.5 p[0] * z_p_base_23[1] - 0.5 p[1] * z_p_base_23[1]
   + 0.5 p[2] * z_p_base_23[1] + 0.5 p[3] * z_p_base_23[1]
   - 0.5 p[4] * z_p_base_23[1] ] <= 0
 infeasibility_p_base_23_row_2_leq_0: 5e-06 z_p_base_23[2] + [
   0.5 p[0] * z_p_base_23[2] - 0.5 p[1] * z_p_base_23[2]
   - 0.5 p[2] * z_p_base_23[2] - 0.5 p[3] * z_p_base_23[2]
   + 0.5 p[4] * z_p_base_23[2] ] <= 0
 infeasibility_p_base_23_row_3_leq_0: 5e-06 z_p_base_23[3] + [
   - 0.5 p[0] * z_p_base_23[3] + 0.5 p[1] * z_p_base_23[3]
   + 0.5 p[2] * z_p_base_23[3] - 0.5 p[3] * z_p_base_23[3]
   + 0.5 p[4] * z_p_base_23[3] ] <= 0
 infeasibility_p_base_23_row_4_leq_0: 5e-06 z_p_base_23[4] + [
   p[3] * z_p_base_23[4] ] <= 0
 infeasibility_p_prime_base_23_row_0_leq_0: 5e-06 z_p_prime_base_23[0] + [
   - alpha[1] * z_p_prime_base_23[0] + p[1] * z_p_prime_base_23[0] ] <= 0
 infeasibility_p_prime_base_23_row_1_leq_0: 5e-06 z_p_prime_base_23[1] + [
   0.5 alpha[0] * z_p_prime_base_23[1]
   + 0.5 alpha[1] * z_p_prime_base_23[1]
   - 0.5 alpha[2] * z_p_prime_base_23[1]
   - 0.5 alpha[3] * z_p_prime_base_23[1]
   - 0.5 alpha[4] * z_p_prime_base_23[1] + 0.5 p[0] * z_p_prime_base_23[1]
   - 0.5 p[1] * z_p_prime_base_23[1] + 0.5 p[2] * z_p_prime_base_23[1]
   + 0.5 p[3] * z_p_prime_base_23[1] - 0.5 p[4] * z_p_prime_base_23[1] ]
   <= 0
 infeasibility_p_prime_base_23_row_2_leq_0: 5e-06 z_p_prime_base_23[2] + [
   0.5 alpha[0] * z_p_prime_base_23[2]
   + 0.5 alpha[1] * z_p_prime_base_23[2]
   + 0.5 alpha[2] * z_p_prime_base_23[2]
   + 0.5 alpha[3] * z_p_prime_base_23[2]
   + 0.5 alpha[4] * z_p_prime_base_23[2] + 0.5 p[0] * z_p_prime_base_23[2]
   - 0.5 p[1] * z_p_prime_base_23[2] - 0.5 p[2] * z_p_prime_base_23[2]
   - 0.5 p[3] * z_p_prime_base_23[2] + 0.5 p[4] * z_p_prime_base_23[2] ]
   <= 0
 infeasibility_p_prime_base_23_row_3_leq_0: 5e-06 z_p_prime_base_23[3] + [
   - 0.5 alpha[0] * z_p_prime_base_23[3]
   - 0.5 alpha[1] * z_p_prime_base_23[3]
   - 0.5 alpha[2] * z_p_prime_base_23[3]
   + 0.5 alpha[3] * z_p_prime_base_23[3]
   + 0.5 alpha[4] * z_p_prime_base_23[3] - 0.5 p[0] * z_p_prime_base_23[3]
   + 0.5 p[1] * z_p_prime_base_23[3] + 0.5 p[2] * z_p_prime_base_23[3]
   - 0.5 p[3] * z_p_prime_base_23[3] + 0.5 p[4] * z_p_prime_base_23[3] ]
   <= 0
 infeasibility_p_prime_base_23_row_4_leq_0: 5e-06 z_p_prime_base_23[4] + [
   - alpha[3] * z_p_prime_base_23[4] + p[3] * z_p_prime_base_23[4] ] <= 0
 infeasibility_p_base_24_row_0_leq_0: 5e-06 z_p_base_24[0] + [
   0.5 p[0] * z_p_base_24[0] + 0.5 p[1] * z_p_base_24[0]
   - 0.5 p[2] * z_p_base_24[0] - 0.5 p[3] * z_p_base_24[0]
   - 0.5 p[4] * z_p_base_24[0] ] <= 0
 infeasibility_p_base_24_row_1_leq_0: 5e-06 z_p_base_24[1] + [
   0.5 p[0] * z_p_base_24[1] - 0.5 p[1] * z_p_base_24[1]
   + 0.5 p[2] * z_p_base_24[1] + 0.5 p[3] * z_p_base_24[1]
   + 0.5 p[4] * z_p_base_24[1] ] <= 0
 infeasibility_p_base_24_row_2_leq_0: 5e-06 z_p_base_24[2] + [
   - 0.5 p[0] * z_p_base_24[2] + 0.5 p[1] * z_p_base_24[2]
   + 0.5 p[2] * z_p_base_24[2] - 0.5 p[3] * z_p_base_24[2]
   - 0.5 p[4] * z_p_base_24[2] ] <= 0
 infeasibility_p_base_24_row_3_leq_0: 5e-06 z_p_base_24[3] + [
   p[3] * z_p_base_24[3] ] <= 0
 infeasibility_p_base_24_row_4_leq_0: 5e-06 z_p_base_24[4] + [
   p[4] * z_p_base_24[4] ] <= 0
 infeasibility_p_prime_base_24_row_0_leq_0: 5e-06 z_p_prime_base_24[0] + [
   0.5 alpha[0] * z_p_prime_base_24[0]
   - 0.5 alpha[1] * z_p_prime_base_24[0]
   + 0.5 alpha[2] * z_p_prime_base_24[0]
   + 0.5 alpha[3] * z_p_prime_base_24[0]
   - 0.5 alpha[4] * z_p_prime_base_24[0] + 0.5 p[0] * z_p_prime_base_24[0]
   + 0.5 p[1] * z_p_prime_base_24[0] - 0.5 p[2] * z_p_prime_base_24[0]
   - 0.5 p[3] * z_p_prime_base_24[0] - 0.5 p[4] * z_p_prime_base_24[0] ]
   <= 0
 infeasibility_p_prime_base_24_row_1_leq_0: 5e-06 z_p_prime_base_24[1] + [
   0.5 alpha[0] * z_p_prime_base_24[1]
   + 0.5 alpha[1] * z_p_prime_base_24[1]
   - 0.5 alpha[2] * z_p_prime_base_24[1]
   - 0.5 alpha[3] * z_p_prime_base_24[1]
   + 0.5 alpha[4] * z_p_prime_base_24[1] + 0.5 p[0] * z_p_prime_base_24[1]
   - 0.5 p[1] * z_p_prime_base_24[1] + 0.5 p[2] * z_p_prime_base_24[1]
   + 0.5 p[3] * z_p_prime_base_24[1] + 0.5 p[4] * z_p_prime_base_24[1] ]
   <= 0
 infeasibility_p_prime_base_24_row_2_leq_0: 5e-06 z_p_prime_base_24[2] + [
   - 0.5 alpha[0] * z_p_prime_base_24[2]
   - 0.5 alpha[1] * z_p_prime_base_24[2]
   - 0.5 alpha[2] * z_p_prime_base_24[2]
   + 0.5 alpha[3] * z_p_prime_base_24[2]
   - 0.5 alpha[4] * z_p_prime_base_24[2] - 0.5 p[0] * z_p_prime_base_24[2]
   + 0.5 p[1] * z_p_prime_base_24[2] + 0.5 p[2] * z_p_prime_base_24[2]
   - 0.5 p[3] * z_p_prime_base_24[2] - 0.5 p[4] * z_p_prime_base_24[2] ]
   <= 0
 infeasibility_p_prime_base_24_row_3_leq_0: 5e-06 z_p_prime_base_24[3] + [
   - alpha[3] * z_p_prime_base_24[3] + p[3] * z_p_prime_base_24[3] ] <= 0
 infeasibility_p_prime_base_24_row_4_leq_0: 5e-06 z_p_prime_base_24[4] + [
   alpha[4] * z_p_prime_base_24[4] + p[4] * z_p_prime_base_24[4] ] <= 0
 infeasibility_p_base_25_row_0_leq_0: 5e-06 z_p_base_25[0] + [
   0.5 p[0] * z_p_base_25[0] + 0.5 p[1] * z_p_base_25[0]
   - 0.5 p[2] * z_p_base_25[0] - 0.5 p[3] * z_p_base_25[0]
   + 0.5 p[4] * z_p_base_25[0] ] <= 0
 infeasibility_p_base_25_row_1_leq_0: 5e-06 z_p_base_25[1] + [
   0.5 p[0] * z_p_base_25[1] - 0.5 p[1] * z_p_base_25[1]
   + 0.5 p[2] * z_p_base_25[1] + 0.5 p[3] * z_p_base_25[1]
   - 0.5 p[4] * z_p_base_25[1] ] <= 0
 infeasibility_p_base_25_row_2_leq_0: 5e-06 z_p_base_25[2] + [
   - 0.5 p[0] * z_p_base_25[2] + 0.5 p[1] * z_p_base_25[2]
   + 0.5 p[2] * z_p_base_25[2] - 0.5 p[3] * z_p_base_25[2]
   - 0.5 p[4] * z_p_base_25[2] ] <= 0
 infeasibility_p_base_25_row_3_leq_0: 5e-06 z_p_base_25[3] + [
   p[3] * z_p_base_25[3] ] <= 0
 infeasibility_p_base_25_row_4_leq_0: 5e-06 z_p_base_25[4] + [
   p[4] * z_p_base_25[4] ] <= 0
 infeasibility_p_prime_base_25_row_0_leq_0: 5e-06 z_p_prime_base_25[0] + [
   0.5 alpha[0] * z_p_prime_base_25[0]
   - 0.5 alpha[1] * z_p_prime_base_25[0]
   + 0.5 alpha[2] * z_p_prime_base_25[0]
   + 0.5 alpha[3] * z_p_prime_base_25[0]
   + 0.5 alpha[4] * z_p_prime_base_25[0] + 0.5 p[0] * z_p_prime_base_25[0]
   + 0.5 p[1] * z_p_prime_base_25[0] - 0.5 p[2] * z_p_prime_base_25[0]
   - 0.5 p[3] * z_p_prime_base_25[0] + 0.5 p[4] * z_p_prime_base_25[0] ]
   <= 0
 infeasibility_p_prime_base_25_row_1_leq_0: 5e-06 z_p_prime_base_25[1] + [
   0.5 alpha[0] * z_p_prime_base_25[1]
   + 0.5 alpha[1] * z_p_prime_base_25[1]
   - 0.5 alpha[2] * z_p_prime_base_25[1]
   - 0.5 alpha[3] * z_p_prime_base_25[1]
   - 0.5 alpha[4] * z_p_prime_base_25[1] + 0.5 p[0] * z_p_prime_base_25[1]
   - 0.5 p[1] * z_p_prime_base_25[1] + 0.5 p[2] * z_p_prime_base_25[1]
   + 0.5 p[3] * z_p_prime_base_25[1] - 0.5 p[4] * z_p_prime_base_25[1] ]
   <= 0
 infeasibility_p_prime_base_25_row_2_leq_0: 5e-06 z_p_prime_base_25[2] + [
   - 0.5 alpha[0] * z_p_prime_base_25[2]
   - 0.5 alpha[1] * z_p_prime_base_25[2]
   - 0.5 alpha[2] * z_p_prime_base_25[2]
   + 0.5 alpha[3] * z_p_prime_base_25[2]
   - 0.5 alpha[4] * z_p_prime_base_25[2] - 0.5 p[0] * z_p_prime_base_25[2]
   + 0.5 p[1] * z_p_prime_base_25[2] + 0.5 p[2] * z_p_prime_base_25[2]
   - 0.5 p[3] * z_p_prime_base_25[2] - 0.5 p[4] * z_p_prime_base_25[2] ]
   <= 0
 infeasibility_p_prime_base_25_row_3_leq_0: 5e-06 z_p_prime_base_25[3] + [
   - alpha[3] * z_p_prime_base_25[3] + p[3] * z_p_prime_base_25[3] ] <= 0
 infeasibility_p_prime_base_25_row_4_leq_0: 5e-06 z_p_prime_base_25[4] + [
   alpha[4] * z_p_prime_base_25[4] + p[4] * z_p_prime_base_25[4] ] <= 0
 infeasibility_p_base_26_row_0_leq_0: 5e-06 z_p_base_26[0] + [
   0.5 p[0] * z_p_base_26[0] + 0.5 p[1] * z_p_base_26[0]
   - 0.5 p[2] * z_p_base_26[0] - 0.5 p[3] * z_p_base_26[0]
   + 0.5 p[4] * z_p_base_26[0] ] <= 0
 infeasibility_p_base_26_row_1_leq_0: 5e-06 z_p_base_26[1] + [
   0.5 p[0] * z_p_base_26[1] - 0.5 p[1] * z_p_base_26[1]
   + 0.5 p[2] * z_p_base_26[1] + 0.5 p[3] * z_p_base_26[1]
   - 0.5 p[4] * z_p_base_26[1] ] <= 0
 infeasibility_p_base_26_row_2_leq_0: 5e-06 z_p_base_26[2] + [
   - 0.5 p[0] * z_p_base_26[2] + 0.5 p[1] * z_p_base_26[2]
   + 0.5 p[2] * z_p_base_26[2] - 0.5 p[3] * z_p_base_26[2]
   + 0.5 p[4] * z_p_base_26[2] ] <= 0
 infeasibility_p_base_26_row_3_leq_0: 5e-06 z_p_base_26[3] + [
   p[3] * z_p_base_26[3] - p[4] * z_p_base_26[3] ] <= 0
 infeasibility_p_base_26_row_4_leq_0: 5e-06 z_p_base_26[4] + [
   p[4] * z_p_base_26[4] ] <= 0
 infeasibility_p_prime_base_26_row_0_leq_0: 5e-06 z_p_prime_base_26[0] + [
   0.5 alpha[0] * z_p_prime_base_26[0]
   - 0.5 alpha[1] * z_p_prime_base_26[0]
   + 0.5 alpha[2] * z_p_prime_base_26[0]
   + 0.5 alpha[3] * z_p_prime_base_26[0]
   + 0.5 alpha[4] * z_p_prime_base_26[0] + 0.5 p[0] * z_p_prime_base_26[0]
   + 0.5 p[1] * z_p_prime_base_26[0] - 0.5 p[2] * z_p_prime_base_26[0]
   - 0.5 p[3] * z_p_prime_base_26[0] + 0.5 p[4] * z_p_prime_base_26[0] ]
   <= 0
 infeasibility_p_prime_base_26_row_1_leq_0: 5e-06 z_p_prime_base_26[1] + [
   0.5 alpha[0] * z_p_prime_base_26[1]
   + 0.5 alpha[1] * z_p_prime_base_26[1]
   - 0.5 alpha[2] * z_p_prime_base_26[1]
   - 0.5 alpha[3] * z_p_prime_base_26[1]
   - 0.5 alpha[4] * z_p_prime_base_26[1] + 0.5 p[0] * z_p_prime_base_26[1]
   - 0.5 p[1] * z_p_prime_base_26[1] + 0.5 p[2] * z_p_prime_base_26[1]
   + 0.5 p[3] * z_p_prime_base_26[1] - 0.5 p[4] * z_p_prime_base_26[1] ]
   <= 0
 infeasibility_p_prime_base_26_row_2_leq_0: 5e-06 z_p_prime_base_26[2] + [
   - 0.5 alpha[0] * z_p_prime_base_26[2]
   - 0.5 alpha[1] * z_p_prime_base_26[2]
   - 0.5 alpha[2] * z_p_prime_base_26[2]
   + 0.5 alpha[3] * z_p_prime_base_26[2]
   + 0.5 alpha[4] * z_p_prime_base_26[2] - 0.5 p[0] * z_p_prime_base_26[2]
   + 0.5 p[1] * z_p_prime_base_26[2] + 0.5 p[2] * z_p_prime_base_26[2]
   - 0.5 p[3] * z_p_prime_base_26[2] + 0.5 p[4] * z_p_prime_base_26[2] ]
   <= 0
 infeasibility_p_prime_base_26_row_3_leq_0: 5e-06 z_p_prime_base_26[3] + [
   - alpha[3] * z_p_prime_base_26[3] - alpha[4] * z_p_prime_base_26[3]
   + p[3] * z_p_prime_base_26[3] - p[4] * z_p_prime_base_26[3] ] <= 0
 infeasibility_p_prime_base_26_row_4_leq_0: 5e-06 z_p_prime_base_26[4] + [
   alpha[4] * z_p_prime_base_26[4] + p[4] * z_p_prime_base_26[4] ] <= 0
 infeasibility_p_base_27_row_0_leq_0: 5e-06 z_p_base_27[0] + [
   0.5 p[0] * z_p_base_27[0] + 0.5 p[1] * z_p_base_27[0]
   - 0.5 p[2] * z_p_base_27[0] + 0.5 p[3] * z_p_base_27[0]
   - 0.5 p[4] * z_p_base_27[0] ] <= 0
 infeasibility_p_base_27_row_1_leq_0: 5e-06 z_p_base_27[1] + [
   0.5 p[0] * z_p_base_27[1] - 0.5 p[1] * z_p_base_27[1]
   + 0.5 p[2] * z_p_base_27[1] - 0.5 p[3] * z_p_base_27[1]
   + 0.5 p[4] * z_p_base_27[1] ] <= 0
 infeasibility_p_base_27_row_2_leq_0: 5e-06 z_p_base_27[2] + [
   - 0.5 p[0] * z_p_base_27[2] + 0.5 p[1] * z_p_base_27[2]
   + 0.5 p[2] * z_p_base_27[2] - 0.5 p[3] * z_p_base_27[2]
   - 0.5 p[4] * z_p_base_27[2] ] <= 0
 infeasibility_p_base_27_row_3_leq_0: 5e-06 z_p_base_27[3] + [
   p[4] * z_p_base_27[3] ] <= 0
 infeasibility_p_base_27_row_4_leq_0: 5e-06 z_p_base_27[4] + [
   p[3] * z_p_base_27[4] ] <= 0
 infeasibility_p_prime_base_27_row_0_leq_0: 5e-06 z_p_prime_base_27[0] + [
   0.5 alpha[0] * z_p_prime_base_27[0]
   - 0.5 alpha[1] * z_p_prime_base_27[0]
   + 0.5 alpha[2] * z_p_prime_base_27[0]
   - 0.5 alpha[3] * z_p_prime_base_27[0]
   - 0.5 alpha[4] * z_p_prime_base_27[0] + 0.5 p[0] * z_p_prime_base_27[0]
   + 0.5 p[1] * z_p_prime_base_27[0] - 0.5 p[2] * z_p_prime_base_27[0]
   + 0.5 p[3] * z_p_prime_base_27[0] - 0.5 p[4] * z_p_prime_base_27[0] ]
   <= 0
 infeasibility_p_prime_base_27_row_1_leq_0: 5e-06 z_p_prime_base_27[1] + [
   0.5 alpha[0] * z_p_prime_base_27[1]
   + 0.5 alpha[1] * z_p_prime_base_27[1]
   - 0.5 alpha[2] * z_p_prime_base_27[1]
   + 0.5 alpha[3] * z_p_prime_base_27[1]
   + 0.5 alpha[4] * z_p_prime_base_27[1] + 0.5 p[0] * z_p_prime_base_27[1]
   - 0.5 p[1] * z_p_prime_base_27[1] + 0.5 p[2] * z_p_prime_base_27[1]
   - 0.5 p[3] * z_p_prime_base_27[1] + 0.5 p[4] * z_p_prime_base_27[1] ]
   <= 0
 infeasibility_p_prime_base_27_row_2_leq_0: 5e-06 z_p_prime_base_27[2] + [
   - 0.5 alpha[0] * z_p_prime_base_27[2]
   - 0.5 alpha[1] * z_p_prime_base_27[2]
   - 0.5 alpha[2] * z_p_prime_base_27[2]
   + 0.5 alpha[3] * z_p_prime_base_27[2]
   - 0.5 alpha[4] * z_p_prime_base_27[2] - 0.5 p[0] * z_p_prime_base_27[2]
   + 0.5 p[1] * z_p_prime_base_27[2] + 0.5 p[2] * z_p_prime_base_27[2]
   - 0.5 p[3] * z_p_prime_base_27[2] - 0.5 p[4] * z_p_prime_base_27[2] ]
   <= 0
 infeasibility_p_prime_base_27_row_3_leq_0: 5e-06 z_p_prime_base_27[3] + [
   alpha[4] * z_p_prime_base_27[3] + p[4] * z_p_prime_base_27[3] ] <= 0
 infeasibility_p_prime_base_27_row_4_leq_0: 5e-06 z_p_prime_base_27[4] + [
   - alpha[3] * z_p_prime_base_27[4] + p[3] * z_p_prime_base_27[4] ] <= 0
 infeasibility_p_base_28_row_0_leq_0: 5e-06 z_p_base_28[0] + [
   0.5 p[0] * z_p_base_28[0] + 0.5 p[1] * z_p_base_28[0]
   - 0.5 p[2] * z_p_base_28[0] + 0.5 p[3] * z_p_base_28[0]
   - 0.5 p[4] * z_p_base_28[0] ] <= 0
 infeasibility_p_base_28_row_1_leq_0: 5e-06 z_p_base_28[1] + [
   0.5 p[0] * z_p_base_28[1] - 0.5 p[1] * z_p_base_28[1]
   + 0.5 p[2] * z_p_base_28[1] - 0.5 p[3] * z_p_base_28[1]
   + 0.5 p[4] * z_p_base_28[1] ] <= 0
 infeasibility_p_base_28_row_2_leq_0: 5e-06 z_p_base_28[2] + [
   - 0.5 p[0] * z_p_base_28[2] + 0.5 p[1] * z_p_base_28[2]
   + 0.5 p[2] * z_p_base_28[2] + 0.5 p[3] * z_p_base_28[2]
   - 0.5 p[4] * z_p_base_28[2] ] <= 0
 infeasibility_p_base_28_row_3_leq_0: 5e-06 z_p_base_28[3] + [
   - p[3] * z_p_base_28[3] + p[4] * z_p_base_28[3] ] <= 0
 infeasibility_p_base_28_row_4_leq_0: 5e-06 z_p_base_28[4] + [
   p[3] * z_p_base_28[4] ] <= 0
 infeasibility_p_prime_base_28_row_0_leq_0: 5e-06 z_p_prime_base_28[0] + [
   0.5 alpha[0] * z_p_prime_base_28[0]
   - 0.5 alpha[1] * z_p_prime_base_28[0]
   + 0.5 alpha[2] * z_p_prime_base_28[0]
   - 0.5 alpha[3] * z_p_prime_base_28[0]
   - 0.5 alpha[4] * z_p_prime_base_28[0] + 0.5 p[0] * z_p_prime_base_28[0]
   + 0.5 p[1] * z_p_prime_base_28[0] - 0.5 p[2] * z_p_prime_base_28[0]
   + 0.5 p[3] * z_p_prime_base_28[0] - 0.5 p[4] * z_p_prime_base_28[0] ]
   <= 0
 infeasibility_p_prime_base_28_row_1_leq_0: 5e-06 z_p_prime_base_28[1] + [
   0.5 alpha[0] * z_p_prime_base_28[1]
   + 0.5 alpha[1] * z_p_prime_base_28[1]
   - 0.5 alpha[2] * z_p_prime_base_28[1]
   + 0.5 alpha[3] * z_p_prime_base_28[1]
   + 0.5 alpha[4] * z_p_prime_base_28[1] + 0.5 p[0] * z_p_prime_base_28[1]
   - 0.5 p[1] * z_p_prime_base_28[1] + 0.5 p[2] * z_p_prime_base_28[1]
   - 0.5 p[3] * z_p_prime_base_28[1] + 0.5 p[4] * z_p_prime_base_28[1] ]
   <= 0
 infeasibility_p_prime_base_28_row_2_leq_0: 5e-06 z_p_prime_base_28[2] + [
   - 0.5 alpha[0] * z_p_prime_base_28[2]
   - 0.5 alpha[1] * z_p_prime_base_28[2]
   - 0.5 alpha[2] * z_p_prime_base_28[2]
   - 0.5 alpha[3] * z_p_prime_base_28[2]
   - 0.5 alpha[4] * z_p_prime_base_28[2] - 0.5 p[0] * z_p_prime_base_28[2]
   + 0.5 p[1] * z_p_prime_base_28[2] + 0.5 p[2] * z_p_prime_base_28[2]
   + 0.5 p[3] * z_p_prime_base_28[2] - 0.5 p[4] * z_p_prime_base_28[2] ]
   <= 0
 infeasibility_p_prime_base_28_row_3_leq_0: 5e-06 z_p_prime_base_28[3] + [
   alpha[3] * z_p_prime_base_28[3] + alpha[4] * z_p_prime_base_28[3]
   - p[3] * z_p_prime_base_28[3] + p[4] * z_p_prime_base_28[3] ] <= 0
 infeasibility_p_prime_base_28_row_4_leq_0: 5e-06 z_p_prime_base_28[4] + [
   - alpha[3] * z_p_prime_base_28[4] + p[3] * z_p_prime_base_28[4] ] <= 0
 infeasibility_p_base_29_row_0_leq_0: 5e-06 z_p_base_29[0] + [
   0.5 p[0] * z_p_base_29[0] + 0.5 p[1] * z_p_base_29[0]
   - 0.5 p[2] * z_p_base_29[0] + 0.5 p[3] * z_p_base_29[0]
   + 0.5 p[4] * z_p_base_29[0] ] <= 0
 infeasibility_p_base_29_row_1_leq_0: 5e-06 z_p_base_29[1] + [
   0.5 p[0] * z_p_base_29[1] - 0.5 p[1] * z_p_base_29[1]
   + 0.5 p[2] * z_p_base_29[1] - 0.5 p[3] * z_p_base_29[1]
   - 0.5 p[4] * z_p_base_29[1] ] <= 0
 infeasibility_p_base_29_row_2_leq_0: 5e-06 z_p_base_29[2] + [
   - 0.5 p[0] * z_p_base_29[2] + 0.5 p[1] * z_p_base_29[2]
   + 0.5 p[2] * z_p_base_29[2] - 0.5 p[3] * z_p_base_29[2]
   - 0.5 p[4] * z_p_base_29[2] ] <= 0
 infeasibility_p_base_29_row_3_leq_0: 5e-06 z_p_base_29[3] + [
   p[3] * z_p_base_29[3] ] <= 0
 infeasibility_p_base_29_row_4_leq_0: 5e-06 z_p_base_29[4] + [
   p[4] * z_p_base_29[4] ] <= 0
 infeasibility_p_prime_base_29_row_0_leq_0: 5e-06 z_p_prime_base_29[0] + [
   0.5 alpha[0] * z_p_prime_base_29[0]
   - 0.5 alpha[1] * z_p_prime_base_29[0]
   + 0.5 alpha[2] * z_p_prime_base_29[0]
   - 0.5 alpha[3] * z_p_prime_base_29[0]
   + 0.5 alpha[4] * z_p_prime_base_29[0] + 0.5 p[0] * z_p_prime_base_29[0]
   + 0.5 p[1] * z_p_prime_base_29[0] - 0.5 p[2] * z_p_prime_base_29[0]
   + 0.5 p[3] * z_p_prime_base_29[0] + 0.5 p[4] * z_p_prime_base_29[0] ]
   <= 0
 infeasibility_p_prime_base_29_row_1_leq_0: 5e-06 z_p_prime_base_29[1] + [
   0.5 alpha[0] * z_p_prime_base_29[1]
   + 0.5 alpha[1] * z_p_prime_base_29[1]
   - 0.5 alpha[2] * z_p_prime_base_29[1]
   + 0.5 alpha[3] * z_p_prime_base_29[1]
   - 0.5 alpha[4] * z_p_prime_base_29[1] + 0.5 p[0] * z_p_prime_base_29[1]
   - 0.5 p[1] * z_p_prime_base_29[1] + 0.5 p[2] * z_p_prime_base_29[1]
   - 0.5 p[3] * z_p_prime_base_29[1] - 0.5 p[4] * z_p_prime_base_29[1] ]
   <= 0
 infeasibility_p_prime_base_29_row_2_leq_0: 5e-06 z_p_prime_base_29[2] + [
   - 0.5 alpha[0] * z_p_prime_base_29[2]
   - 0.5 alpha[1] * z_p_prime_base_29[2]
   - 0.5 alpha[2] * z_p_prime_base_29[2]
   + 0.5 alpha[3] * z_p_prime_base_29[2]
   - 0.5 alpha[4] * z_p_prime_base_29[2] - 0.5 p[0] * z_p_prime_base_29[2]
   + 0.5 p[1] * z_p_prime_base_29[2] + 0.5 p[2] * z_p_prime_base_29[2]
   - 0.5 p[3] * z_p_prime_base_29[2] - 0.5 p[4] * z_p_prime_base_29[2] ]
   <= 0
 infeasibility_p_prime_base_29_row_3_leq_0: 5e-06 z_p_prime_base_29[3] + [
   - alpha[3] * z_p_prime_base_29[3] + p[3] * z_p_prime_base_29[3] ] <= 0
 infeasibility_p_prime_base_29_row_4_leq_0: 5e-06 z_p_prime_base_29[4] + [
   alpha[4] * z_p_prime_base_29[4] + p[4] * z_p_prime_base_29[4] ] <= 0
 infeasibility_p_base_30_row_0_leq_0: 5e-06 z_p_base_30[0] + [
   0.5 p[0] * z_p_base_30[0] + 0.5 p[1] * z_p_base_30[0]
   - 0.5 p[2] * z_p_base_30[0] + 0.5 p[3] * z_p_base_30[0]
   - 0.5 p[4] * z_p_base_30[0] ] <= 0
 infeasibility_p_base_30_row_1_leq_0: 5e-06 z_p_base_30[1] + [
   0.5 p[0] * z_p_base_30[1] - 0.5 p[1] * z_p_base_30[1]
   + 0.5 p[2] * z_p_base_30[1] - 0.5 p[3] * z_p_base_30[1]
   + 0.5 p[4] * z_p_base_30[1] ] <= 0
 infeasibility_p_base_30_row_2_leq_0: 5e-06 z_p_base_30[2] + [
   - 0.5 p[0] * z_p_base_30[2] + 0.5 p[1] * z_p_base_30[2]
   + 0.5 p[2] * z_p_base_30[2] - 0.5 p[3] * z_p_base_30[2]
   + 0.5 p[4] * z_p_base_30[2] ] <= 0
 infeasibility_p_base_30_row_3_leq_0: 5e-06 z_p_base_30[3] + [
   p[3] * z_p_base_30[3] - p[4] * z_p_base_30[3] ] <= 0
 infeasibility_p_base_30_row_4_leq_0: 5e-06 z_p_base_30[4] + [
   p[4] * z_p_base_30[4] ] <= 0
 infeasibility_p_prime_base_30_row_0_leq_0: 5e-06 z_p_prime_base_30[0] + [
   0.5 alpha[0] * z_p_prime_base_30[0]
   - 0.5 alpha[1] * z_p_prime_base_30[0]
   + 0.5 alpha[2] * z_p_prime_base_30[0]
   - 0.5 alpha[3] * z_p_prime_base_30[0]
   - 0.5 alpha[4] * z_p_prime_base_30[0] + 0.5 p[0] * z_p_prime_base_30[0]
   + 0.5 p[1] * z_p_prime_base_30[0] - 0.5 p[2] * z_p_prime_base_30[0]
   + 0.5 p[3] * z_p_prime_base_30[0] - 0.5 p[4] * z_p_prime_base_30[0] ]
   <= 0
 infeasibility_p_prime_base_30_row_1_leq_0: 5e-06 z_p_prime_base_30[1] + [
   0.5 alpha[0] * z_p_prime_base_30[1]
   + 0.5 alpha[1] * z_p_prime_base_30[1]
   - 0.5 alpha[2] * z_p_prime_base_30[1]
   + 0.5 alpha[3] * z_p_prime_base_30[1]
   + 0.5 alpha[4] * z_p_prime_base_30[1] + 0.5 p[0] * z_p_prime_base_30[1]
   - 0.5 p[1] * z_p_prime_base_30[1] + 0.5 p[2] * z_p_prime_base_30[1]
   - 0.5 p[3] * z_p_prime_base_30[1] + 0.5 p[4] * z_p_prime_base_30[1] ]
   <= 0
 infeasibility_p_prime_base_30_row_2_leq_0: 5e-06 z_p_prime_base_30[2] + [
   - 0.5 alpha[0] * z_p_prime_base_30[2]
   - 0.5 alpha[1] * z_p_prime_base_30[2]
   - 0.5 alpha[2] * z_p_prime_base_30[2]
   + 0.5 alpha[3] * z_p_prime_base_30[2]
   + 0.5 alpha[4] * z_p_prime_base_30[2] - 0.5 p[0] * z_p_prime_base_30[2]
   + 0.5 p[1] * z_p_prime_base_30[2] + 0.5 p[2] * z_p_prime_base_30[2]
   - 0.5 p[3] * z_p_prime_base_30[2] + 0.5 p[4] * z_p_prime_base_30[2] ]
   <= 0
 infeasibility_p_prime_base_30_row_3_leq_0: 5e-06 z_p_prime_base_30[3] + [
   - alpha[3] * z_p_prime_base_30[3] - alpha[4] * z_p_prime_base_30[3]
   + p[3] * z_p_prime_base_30[3] - p[4] * z_p_prime_base_30[3] ] <= 0
 infeasibility_p_prime_base_30_row_4_leq_0: 5e-06 z_p_prime_base_30[4] + [
   alpha[4] * z_p_prime_base_30[4] + p[4] * z_p_prime_base_30[4] ] <= 0
 infeasibility_p_base_31_row_0_leq_0: 5e-06 z_p_base_31[0] + [
   0.5 p[0] * z_p_base_31[0] + 0.5 p[1] * z_p_base_31[0]
   - 0.5 p[2] * z_p_base_31[0] - 0.5 p[3] * z_p_base_31[0]
   + 0.5 p[4] * z_p_base_31[0] ] <= 0
 infeasibility_p_base_31_row_1_leq_0: 5e-06 z_p_base_31[1] + [
   0.5 p[0] * z_p_base_31[1] - 0.5 p[1] * z_p_base_31[1]
   + 0.5 p[2] * z_p_base_31[1] + 0.5 p[3] * z_p_base_31[1]
   - 0.5 p[4] * z_p_base_31[1] ] <= 0
 infeasibility_p_base_31_row_2_leq_0: 5e-06 z_p_base_31[2] + [
   - 0.5 p[0] * z_p_base_31[2] + 0.5 p[1] * z_p_base_31[2]
   + 0.5 p[2] * z_p_base_31[2] + 0.5 p[3] * z_p_base_31[2]
   - 0.5 p[4] * z_p_base_31[2] ] <= 0
 infeasibility_p_base_31_row_3_leq_0: 5e-06 z_p_base_31[3] + [
   - p[3] * z_p_base_31[3] + p[4] * z_p_base_31[3] ] <= 0
 infeasibility_p_base_31_row_4_leq_0: 5e-06 z_p_base_31[4] + [
   p[3] * z_p_base_31[4] ] <= 0
 infeasibility_p_prime_base_31_row_0_leq_0: 5e-06 z_p_prime_base_31[0] + [
   0.5 alpha[0] * z_p_prime_base_31[0]
   - 0.5 alpha[1] * z_p_prime_base_31[0]
   + 0.5 alpha[2] * z_p_prime_base_31[0]
   + 0.5 alpha[3] * z_p_prime_base_31[0]
   + 0.5 alpha[4] * z_p_prime_base_31[0] + 0.5 p[0] * z_p_prime_base_31[0]
   + 0.5 p[1] * z_p_prime_base_31[0] - 0.5 p[2] * z_p_prime_base_31[0]
   - 0.5 p[3] * z_p_prime_base_31[0] + 0.5 p[4] * z_p_prime_base_31[0] ]
   <= 0
 infeasibility_p_prime_base_31_row_1_leq_0: 5e-06 z_p_prime_base_31[1] + [
   0.5 alpha[0] * z_p_prime_base_31[1]
   + 0.5 alpha[1] * z_p_prime_base_31[1]
   - 0.5 alpha[2] * z_p_prime_base_31[1]
   - 0.5 alpha[3] * z_p_prime_base_31[1]
   - 0.5 alpha[4] * z_p_prime_base_31[1] + 0.5 p[0] * z_p_prime_base_31[1]
   - 0.5 p[1] * z_p_prime_base_31[1] + 0.5 p[2] * z_p_prime_base_31[1]
   + 0.5 p[3] * z_p_prime_base_31[1] - 0.5 p[4] * z_p_prime_base_31[1] ]
   <= 0
 infeasibility_p_prime_base_31_row_2_leq_0: 5e-06 z_p_prime_base_31[2] + [
   - 0.5 alpha[0] * z_p_prime_base_31[2]
   - 0.5 alpha[1] * z_p_prime_base_31[2]
   - 0.5 alpha[2] * z_p_prime_base_31[2]
   - 0.5 alpha[3] * z_p_prime_base_31[2]
   - 0.5 alpha[4] * z_p_prime_base_31[2] - 0.5 p[0] * z_p_prime_base_31[2]
   + 0.5 p[1] * z_p_prime_base_31[2] + 0.5 p[2] * z_p_prime_base_31[2]
   + 0.5 p[3] * z_p_prime_base_31[2] - 0.5 p[4] * z_p_prime_base_31[2] ]
   <= 0
 infeasibility_p_prime_base_31_row_3_leq_0: 5e-06 z_p_prime_base_31[3] + [
   alpha[3] * z_p_prime_base_31[3] + alpha[4] * z_p_prime_base_31[3]
   - p[3] * z_p_prime_base_31[3] + p[4] * z_p_prime_base_31[3] ] <= 0
 infeasibility_p_prime_base_31_row_4_leq_0: 5e-06 z_p_prime_base_31[4] + [
   - alpha[3] * z_p_prime_base_31[4] + p[3] * z_p_prime_base_31[4] ] <= 0
 infeasibility_p_base_32_row_0_leq_0: 5e-06 z_p_base_32[0] + [
   p[0] * z_p_base_32[0] - p[2] * z_p_base_32[0] ] <= 0
 infeasibility_p_base_32_row_1_leq_0: 5e-06 z_p_base_32[1] + [
   p[2] * z_p_base_32[1] ] <= 0
 infeasibility_p_base_32_row_2_leq_0: 5e-06 z_p_base_32[2] + [
   - 0.5 p[0] * z_p_base_32[2] + 0.5 p[1] * z_p_base_32[2]
   + 0.5 p[2] * z_p_base_32[2] + 0.5 p[3] * z_p_base_32[2]
   - 0.5 p[4] * z_p_base_32[2] ] <= 0
 infeasibility_p_base_32_row_3_leq_0: 5e-06 z_p_base_32[3] + [
   - 0.5 p[0] * z_p_base_32[3] + 0.5 p[1] * z_p_base_32[3]
   + 0.5 p[2] * z_p_base_32[3] - 0.5 p[3] * z_p_base_32[3]
   + 0.5 p[4] * z_p_base_32[3] ] <= 0
 infeasibility_p_base_32_row_4_leq_0: 5e-06 z_p_base_32[4] + [
   0.5 p[0] * z_p_base_32[4] - 0.5 p[1] * z_p_base_32[4]
   - 0.5 p[2] * z_p_base_32[4] + 0.5 p[3] * z_p_base_32[4]
   + 0.5 p[4] * z_p_base_32[4] ] <= 0
 infeasibility_p_prime_base_32_row_0_leq_0: 5e-06 z_p_prime_base_32[0] + [
   alpha[0] * z_p_prime_base_32[0] + alpha[2] * z_p_prime_base_32[0]
   + p[0] * z_p_prime_base_32[0] - p[2] * z_p_prime_base_32[0] ] <= 0
 infeasibility_p_prime_base_32_row_1_leq_0: 5e-06 z_p_prime_base_32[1] + [
   - alpha[2] * z_p_prime_base_32[1] + p[2] * z_p_prime_base_32[1] ] <= 0
 infeasibility_p_prime_base_32_row_2_leq_0: 5e-06 z_p_prime_base_32[2] + [
   - 0.5 alpha[0] * z_p_prime_base_32[2]
   - 0.5 alpha[1] * z_p_prime_base_32[2]
   - 0.5 alpha[2] * z_p_prime_base_32[2]
   - 0.5 alpha[3] * z_p_prime_base_32[2]
   - 0.5 alpha[4] * z_p_prime_base_32[2] - 0.5 p[0] * z_p_prime_base_32[2]
   + 0.5 p[1] * z_p_prime_base_32[2] + 0.5 p[2] * z_p_prime_base_32[2]
   + 0.5 p[3] * z_p_prime_base_32[2] - 0.5 p[4] * z_p_prime_base_32[2] ]
   <= 0
 infeasibility_p_prime_base_32_row_3_leq_0: 5e-06 z_p_prime_base_32[3] + [
   - 0.5 alpha[0] * z_p_prime_base_32[3]
   - 0.5 alpha[1] * z_p_prime_base_32[3]
   - 0.5 alpha[2] * z_p_prime_base_32[3]
   + 0.5 alpha[3] * z_p_prime_base_32[3]
   + 0.5 alpha[4] * z_p_prime_base_32[3] - 0.5 p[0] * z_p_prime_base_32[3]
   + 0.5 p[1] * z_p_prime_base_32[3] + 0.5 p[2] * z_p_prime_base_32[3]
   - 0.5 p[3] * z_p_prime_base_32[3] + 0.5 p[4] * z_p_prime_base_32[3] ]
   <= 0
 infeasibility_p_prime_base_32_row_4_leq_0: 5e-06 z_p_prime_base_32[4] + [
   0.5 alpha[0] * z_p_prime_base_32[4]
   + 0.5 alpha[1] * z_p_prime_base_32[4]
   + 0.5 alpha[2] * z_p_prime_base_32[4]
   - 0.5 alpha[3] * z_p_prime_base_32[4]
   + 0.5 alpha[4] * z_p_prime_base_32[4] + 0.5 p[0] * z_p_prime_base_32[4]
   - 0.5 p[1] * z_p_prime_base_32[4] - 0.5 p[2] * z_p_prime_base_32[4]
   + 0.5 p[3] * z_p_prime_base_32[4] + 0.5 p[4] * z_p_prime_base_32[4] ]
   <= 0
 infeasibility_p_base_33_row_0_leq_0: 5e-06 z_p_base_33[0] + [
   0.5 p[0] * z_p_base_33[0] + 0.5 p[1] * z_p_base_33[0]
   - 0.5 p[2] * z_p_base_33[0] - 0.5 p[3] * z_p_base_33[0]
   + 0.5 p[4] * z_p_base_33[0] ] <= 0
 infeasibility_p_base_33_row_1_leq_0: 5e-06 z_p_base_33[1] + [
   0.5 p[0] * z_p_base_33[1] - 0.5 p[1] * z_p_base_33[1]
   + 0.5 p[2] * z_p_base_33[1] + 0.5 p[3] * z_p_base_33[1]
   - 0.5 p[4] * z_p_base_33[1] ] <= 0
 infeasibility_p_base_33_row_2_leq_0: 5e-06 z_p_base_33[2] + [
   - 0.5 p[0] * z_p_base_33[2] + 0.5 p[1] * z_p_base_33[2]
   + 0.5 p[2] * z_p_base_33[2] + 0.5 p[3] * z_p_base_33[2]
   - 0.5 p[4] * z_p_base_33[2] ] <= 0
 infeasibility_p_base_33_row_3_leq_0: 5e-06 z_p_base_33[3] + [
   - 0.5 p[0] * z_p_base_33[3] + 0.5 p[1] * z_p_base_33[3]
   + 0.5 p[2] * z_p_base_33[3] - 0.5 p[3] * z_p_base_33[3]
   + 0.5 p[4] * z_p_base_33[3] ] <= 0
 infeasibility_p_base_33_row_4_leq_0: 5e-06 z_p_base_33[4] + [
   0.5 p[0] * z_p_base_33[4] - 0.5 p[1] * z_p_base_33[4]
   - 0.5 p[2] * z_p_base_33[4] + 0.5 p[3] * z_p_base_33[4]
   + 0.5 p[4] * z_p_base_33[4] ] <= 0
 infeasibility_p_prime_base_33_row_0_leq_0: 5e-06 z_p_prime_base_33[0] + [
   0.5 alpha[0] * z_p_prime_base_33[0]
   - 0.5 alpha[1] * z_p_prime_base_33[0]
   + 0.5 alpha[2] * z_p_prime_base_33[0]
   + 0.5 alpha[3] * z_p_prime_base_33[0]
   + 0.5 alpha[4] * z_p_prime_base_33[0] + 0.5 p[0] * z_p_prime_base_33[0]
   + 0.5 p[1] * z_p_prime_base_33[0] - 0.5 p[2] * z_p_prime_base_33[0]
   - 0.5 p[3] * z_p_prime_base_33[0] + 0.5 p[4] * z_p_prime_base_33[0] ]
   <= 0
 infeasibility_p_prime_base_33_row_1_leq_0: 5e-06 z_p_prime_base_33[1] + [
   0.5 alpha[0] * z_p_prime_base_33[1]
   + 0.5 alpha[1] * z_p_prime_base_33[1]
   - 0.5 alpha[2] * z_p_prime_base_33[1]
   - 0.5 alpha[3] * z_p_prime_base_33[1]
   - 0.5 alpha[4] * z_p_prime_base_33[1] + 0.5 p[0] * z_p_prime_base_33[1]
   - 0.5 p[1] * z_p_prime_base_33[1] + 0.5 p[2] * z_p_prime_base_33[1]
   + 0.5 p[3] * z_p_prime_base_33[1] - 0.5 p[4] * z_p_prime_base_33[1] ]
   <= 0
 infeasibility_p_prime_base_33_row_2_leq_0: 5e-06 z_p_prime_base_33[2] + [
   - 0.5 alpha[0] * z_p_prime_base_33[2]
   - 0.5 alpha[1] * z_p_prime_base_33[2]
   - 0.5 alpha[2] * z_p_prime_base_33[2]
   - 0.5 alpha[3] * z_p_prime_base_33[2]
   - 0.5 alpha[4] * z_p_prime_base_33[2] - 0.5 p[0] * z_p_prime_base_33[2]
   + 0.5 p[1] * z_p_prime_base_33[2] + 0.5 p[2] * z_p_prime_base_33[2]
   + 0.5 p[3] * z_p_prime_base_33[2] - 0.5 p[4] * z_p_prime_base_33[2] ]
   <= 0
 infeasibility_p_prime_base_33_row_3_leq_0: 5e-06 z_p_prime_base_33[3] + [
   - 0.5 alpha[0] * z_p_prime_base_33[3]
   - 0.5 alpha[1] * z_p_prime_base_33[3]
   - 0.5 alpha[2] * z_p_prime_base_33[3]
   + 0.5 alpha[3] * z_p_prime_base_33[3]
   + 0.5 alpha[4] * z_p_prime_base_33[3] - 0.5 p[0] * z_p_prime_base_33[3]
   + 0.5 p[1] * z_p_prime_base_33[3] + 0.5 p[2] * z_p_prime_base_33[3]
   - 0.5 p[3] * z_p_prime_base_33[3] + 0.5 p[4] * z_p_prime_base_33[3] ]
   <= 0
 infeasibility_p_prime_base_33_row_4_leq_0: 5e-06 z_p_prime_base_33[4] + [
   0.5 alpha[0] * z_p_prime_base_33[4]
   + 0.5 alpha[1] * z_p_prime_base_33[4]
   + 0.5 alpha[2] * z_p_prime_base_33[4]
   - 0.5 alpha[3] * z_p_prime_base_33[4]
   + 0.5 alpha[4] * z_p_prime_base_33[4] + 0.5 p[0] * z_p_prime_base_33[4]
   - 0.5 p[1] * z_p_prime_base_33[4] - 0.5 p[2] * z_p_prime_base_33[4]
   + 0.5 p[3] * z_p_prime_base_33[4] + 0.5 p[4] * z_p_prime_base_33[4] ]
   <= 0
 infeasibility_p_base_34_row_0_leq_0: 5e-06 z_p_base_34[0] + [
   0.5 p[0] * z_p_base_34[0] + 0.5 p[1] * z_p_base_34[0]
   - 0.5 p[2] * z_p_base_34[0] + 0.5 p[3] * z_p_base_34[0]
   - 0.5 p[4] * z_p_base_34[0] ] <= 0
 infeasibility_p_base_34_row_1_leq_0: 5e-06 z_p_base_34[1] + [
   0.5 p[0] * z_p_base_34[1] - 0.5 p[1] * z_p_base_34[1]
   + 0.5 p[2] * z_p_base_34[1] - 0.5 p[3] * z_p_base_34[1]
   + 0.5 p[4] * z_p_base_34[1] ] <= 0
 infeasibility_p_base_34_row_2_leq_0: 5e-06 z_p_base_34[2] + [
   - 0.5 p[0] * z_p_base_34[2] + 0.5 p[1] * z_p_base_34[2]
   + 0.5 p[2] * z_p_base_34[2] - 0.5 p[3] * z_p_base_34[2]
   + 0.5 p[4] * z_p_base_34[2] ] <= 0
 infeasibility_p_base_34_row_3_leq_0: 5e-06 z_p_base_34[3] + [
   - 0.5 p[0] * z_p_base_34[3] + 0.5 p[1] * z_p_base_34[3]
   + 0.5 p[2] * z_p_base_34[3] + 0.5 p[3] * z_p_base_34[3]
   - 0.5 p[4] * z_p_base_34[3] ] <= 0
 infeasibility_p_base_34_row_4_leq_0: 5e-06 z_p_base_34[4] + [
   0.5 p[0] * z_p_base_34[4] - 0.5 p[1] * z_p_base_34[4]
   - 0.5 p[2] * z_p_base_34[4] + 0.5 p[3] * z_p_base_34[4]
   + 0.5 p[4] * z_p_base_34[4] ] <= 0
 infeasibility_p_prime_base_34_row_0_leq_0: 5e-06 z_p_prime_base_34[0] + [
   0.5 alpha[0] * z_p_prime_base_34[0]
   - 0.5 alpha[1] * z_p_prime_base_34[0]
   + 0.5 alpha[2] * z_p_prime_base_34[0]
   - 0.5 alpha[3] * z_p_prime_base_34[0]
   - 0.5 alpha[4] * z_p_prime_base_34[0] + 0.5 p[0] * z_p_prime_base_34[0]
   + 0.5 p[1] * z_p_prime_base_34[0] - 0.5 p[2] * z_p_prime_base_34[0]
   + 0.5 p[3] * z_p_prime_base_34[0] - 0.5 p[4] * z_p_prime_base_34[0] ]
   <= 0
 infeasibility_p_prime_base_34_row_1_leq_0: 5e-06 z_p_prime_base_34[1] + [
   0.5 alpha[0] * z_p_prime_base_34[1]
   + 0.5 alpha[1] * z_p_prime_base_34[1]
   - 0.5 alpha[2] * z_p_prime_base_34[1]
   + 0.5 alpha[3] * z_p_prime_base_34[1]
   + 0.5 alpha[4] * z_p_prime_base_34[1] + 0.5 p[0] * z_p_prime_base_34[1]
   - 0.5 p[1] * z_p_prime_base_34[1] + 0.5 p[2] * z_p_prime_base_34[1]
   - 0.5 p[3] * z_p_prime_base_34[1] + 0.5 p[4] * z_p_prime_base_34[1] ]
   <= 0
 infeasibility_p_prime_base_34_row_2_leq_0: 5e-06 z_p_prime_base_34[2] + [
   - 0.5 alpha[0] * z_p_prime_base_34[2]
   - 0.5 alpha[1] * z_p_prime_base_34[2]
   - 0.5 alpha[2] * z_p_prime_base_34[2]
   + 0.5 alpha[3] * z_p_prime_base_34[2]
   + 0.5 alpha[4] * z_p_prime_base_34[2] - 0.5 p[0] * z_p_prime_base_34[2]
   + 0.5 p[1] * z_p_prime_base_34[2] + 0.5 p[2] * z_p_prime_base_34[2]
   - 0.5 p[3] * z_p_prime_base_34[2] + 0.5 p[4] * z_p_prime_base_34[2] ]
   <= 0
 infeasibility_p_prime_base_34_row_3_leq_0: 5e-06 z_p_prime_base_34[3] + [
   - 0.5 alpha[0] * z_p_prime_base_34[3]
   - 0.5 alpha[1] * z_p_prime_base_34[3]
   - 0.5 alpha[2] * z_p_prime_base_34[3]
   - 0.5 alpha[3] * z_p_prime_base_34[3]
   - 0.5 alpha[4] * z_p_prime_base_34[3] - 0.5 p[0] * z_p_prime_base_34[3]
   + 0.5 p[1] * z_p_prime_base_34[3] + 0.5 p[2] * z_p_prime_base_34[3]
   + 0.5 p[3] * z_p_prime_base_34[3] - 0.5 p[4] * z_p_prime_base_34[3] ]
   <= 0
 infeasibility_p_prime_base_34_row_4_leq_0: 5e-06 z_p_prime_base_34[4] + [
   0.5 alpha[0] * z_p_prime_base_34[4]
   + 0.5 alpha[1] * z_p_prime_base_34[4]
   + 0.5 alpha[2] * z_p_prime_base_34[4]
   - 0.5 alpha[3] * z_p_prime_base_34[4]
   + 0.5 alpha[4] * z_p_prime_base_34[4] + 0.5 p[0] * z_p_prime_base_34[4]
   - 0.5 p[1] * z_p_prime_base_34[4] - 0.5 p[2] * z_p_prime_base_34[4]
   + 0.5 p[3] * z_p_prime_base_34[4] + 0.5 p[4] * z_p_prime_base_34[4] ]
   <= 0
 infeasibility_p_base_35_row_0_leq_0: 5e-06 z_p_base_35[0] + [
   p[1] * z_p_base_35[0] ] <= 0
 infeasibility_p_base_35_row_1_leq_0: 5e-06 z_p_base_35[1] + [
   p[0] * z_p_base_35[1] - p[1] * z_p_base_35[1] ] <= 0
 infeasibility_p_base_35_row_2_leq_0: 5e-06 z_p_base_35[2] + [
   - 0.5 p[0] * z_p_base_35[2] + 0.5 p[1] * z_p_base_35[2]
   + 0.5 p[2] * z_p_base_35[2] + 0.5 p[3] * z_p_base_35[2]
   - 0.5 p[4] * z_p_base_35[2] ] <= 0
 infeasibility_p_base_35_row_3_leq_0: 5e-06 z_p_base_35[3] + [
   - 0.5 p[0] * z_p_base_35[3] + 0.5 p[1] * z_p_base_35[3]
   + 0.5 p[2] * z_p_base_35[3] - 0.5 p[3] * z_p_base_35[3]
   + 0.5 p[4] * z_p_base_35[3] ] <= 0
 infeasibility_p_base_35_row_4_leq_0: 5e-06 z_p_base_35[4] + [
   0.5 p[0] * z_p_base_35[4] - 0.5 p[1] * z_p_base_35[4]
   - 0.5 p[2] * z_p_base_35[4] + 0.5 p[3] * z_p_base_35[4]
   + 0.5 p[4] * z_p_base_35[4] ] <= 0
 infeasibility_p_prime_base_35_row_0_leq_0: 5e-06 z_p_prime_base_35[0] + [
   - alpha[1] * z_p_prime_base_35[0] + p[1] * z_p_prime_base_35[0] ] <= 0
 infeasibility_p_prime_base_35_row_1_leq_0: 5e-06 z_p_prime_base_35[1] + [
   alpha[0] * z_p_prime_base_35[1] + alpha[1] * z_p_prime_base_35[1]
   + p[0] * z_p_prime_base_35[1] - p[1] * z_p_prime_base_35[1] ] <= 0
 infeasibility_p_prime_base_35_row_2_leq_0: 5e-06 z_p_prime_base_35[2] + [
   - 0.5 alpha[0] * z_p_prime_base_35[2]
   - 0.5 alpha[1] * z_p_prime_base_35[2]
   - 0.5 alpha[2] * z_p_prime_base_35[2]
   - 0.5 alpha[3] * z_p_prime_base_35[2]
   - 0.5 alpha[4] * z_p_prime_base_35[2] - 0.5 p[0] * z_p_prime_base_35[2]
   + 0.5 p[1] * z_p_prime_base_35[2] + 0.5 p[2] * z_p_prime_base_35[2]
   + 0.5 p[3] * z_p_prime_base_35[2] - 0.5 p[4] * z_p_prime_base_35[2] ]
   <= 0
 infeasibility_p_prime_base_35_row_3_leq_0: 5e-06 z_p_prime_base_35[3] + [
   - 0.5 alpha[0] * z_p_prime_base_35[3]
   - 0.5 alpha[1] * z_p_prime_base_35[3]
   - 0.5 alpha[2] * z_p_prime_base_35[3]
   + 0.5 alpha[3] * z_p_prime_base_35[3]
   + 0.5 alpha[4] * z_p_prime_base_35[3] - 0.5 p[0] * z_p_prime_base_35[3]
   + 0.5 p[1] * z_p_prime_base_35[3] + 0.5 p[2] * z_p_prime_base_35[3]
   - 0.5 p[3] * z_p_prime_base_35[3] + 0.5 p[4] * z_p_prime_base_35[3] ]
   <= 0
 infeasibility_p_prime_base_35_row_4_leq_0: 5e-06 z_p_prime_base_35[4] + [
   0.5 alpha[0] * z_p_prime_base_35[4]
   + 0.5 alpha[1] * z_p_prime_base_35[4]
   + 0.5 alpha[2] * z_p_prime_base_35[4]
   - 0.5 alpha[3] * z_p_prime_base_35[4]
   + 0.5 alpha[4] * z_p_prime_base_35[4] + 0.5 p[0] * z_p_prime_base_35[4]
   - 0.5 p[1] * z_p_prime_base_35[4] - 0.5 p[2] * z_p_prime_base_35[4]
   + 0.5 p[3] * z_p_prime_base_35[4] + 0.5 p[4] * z_p_prime_base_35[4] ]
   <= 0
 infeasibility_p_base_36_row_0_leq_0: 5e-06 z_p_base_36[0] + [
   0.5 p[0] * z_p_base_36[0] + 0.5 p[1] * z_p_base_36[0]
   - 0.5 p[2] * z_p_base_36[0] - 0.5 p[3] * z_p_base_36[0]
   - 0.5 p[4] * z_p_base_36[0] ] <= 0
 infeasibility_p_base_36_row_1_leq_0: 5e-06 z_p_base_36[1] + [
   0.5 p[0] * z_p_base_36[1] - 0.5 p[1] * z_p_base_36[1]
   + 0.5 p[2] * z_p_base_36[1] + 0.5 p[3] * z_p_base_36[1]
   - 0.5 p[4] * z_p_base_36[1] ] <= 0
 infeasibility_p_base_36_row_2_leq_0: 5e-06 z_p_base_36[2] + [
   p[4] * z_p_base_36[2] ] <= 0
 infeasibility_p_base_36_row_3_leq_0: 5e-06 z_p_base_36[3] + [
   p[2] * z_p_base_36[3] ] <= 0
 infeasibility_p_base_36_row_4_leq_0: 5e-06 z_p_base_36[4] + [
   - 0.5 p[0] * z_p_base_36[4] + 0.5 p[1] * z_p_base_36[4]
   - 0.5 p[2] * z_p_base_36[4] + 0.5 p[3] * z_p_base_36[4]
   + 0.5 p[4] * z_p_base_36[4] ] <= 0
 infeasibility_p_prime_base_36_row_0_leq_0: 5e-06 z_p_prime_base_36[0] + [
   0.5 alpha[0] * z_p_prime_base_36[0]
   - 0.5 alpha[1] * z_p_prime_base_36[0]
   + 0.5 alpha[2] * z_p_prime_base_36[0]
   + 0.5 alpha[3] * z_p_prime_base_36[0]
   - 0.5 alpha[4] * z_p_prime_base_36[0] + 0.5 p[0] * z_p_prime_base_36[0]
   + 0.5 p[1] * z_p_prime_base_36[0] - 0.5 p[2] * z_p_prime_base_36[0]
   - 0.5 p[3] * z_p_prime_base_36[0] - 0.5 p[4] * z_p_prime_base_36[0] ]
   <= 0
 infeasibility_p_prime_base_36_row_1_leq_0: 5e-06 z_p_prime_base_36[1] + [
   0.5 alpha[0] * z_p_prime_base_36[1]
   + 0.5 alpha[1] * z_p_prime_base_36[1]
   - 0.5 alpha[2] * z_p_prime_base_36[1]
   - 0.5 alpha[3] * z_p_prime_base_36[1]
   - 0.5 alpha[4] * z_p_prime_base_36[1] + 0.5 p[0] * z_p_prime_base_36[1]
   - 0.5 p[1] * z_p_prime_base_36[1] + 0.5 p[2] * z_p_prime_base_36[1]
   + 0.5 p[3] * z_p_prime_base_36[1] - 0.5 p[4] * z_p_prime_base_36[1] ]
   <= 0
 infeasibility_p_prime_base_36_row_2_leq_0: 5e-06 z_p_prime_base_36[2] + [
   alpha[4] * z_p_prime_base_36[2] + p[4] * z_p_prime_base_36[2] ] <= 0
 infeasibility_p_prime_base_36_row_3_leq_0: 5e-06 z_p_prime_base_36[3] + [
   - alpha[2] * z_p_prime_base_36[3] + p[2] * z_p_prime_base_36[3] ] <= 0
 infeasibility_p_prime_base_36_row_4_leq_0: 5e-06 z_p_prime_base_36[4] + [
   - 0.5 alpha[0] * z_p_prime_base_36[4]
   - 0.5 alpha[1] * z_p_prime_base_36[4]
   + 0.5 alpha[2] * z_p_prime_base_36[4]
   - 0.5 alpha[3] * z_p_prime_base_36[4]
   + 0.5 alpha[4] * z_p_prime_base_36[4] - 0.5 p[0] * z_p_prime_base_36[4]
   + 0.5 p[1] * z_p_prime_base_36[4] - 0.5 p[2] * z_p_prime_base_36[4]
   + 0.5 p[3] * z_p_prime_base_36[4] + 0.5 p[4] * z_p_prime_base_36[4] ]
   <= 0
 infeasibility_p_base_37_row_0_leq_0: 5e-06 z_p_base_37[0] + [
   0.5 p[0] * z_p_base_37[0] + 0.5 p[1] * z_p_base_37[0]
   - 0.5 p[2] * z_p_base_37[0] - 0.5 p[3] * z_p_base_37[0]
   - 0.5 p[4] * z_p_base_37[0] ] <= 0
 infeasibility_p_base_37_row_1_leq_0: 5e-06 z_p_base_37[1] + [
   p[3] * z_p_base_37[1] ] <= 0
 infeasibility_p_base_37_row_2_leq_0: 5e-06 z_p_base_37[2] + [
   0.5 p[0] * z_p_base_37[2] - 0.5 p[1] * z_p_base_37[2]
   + 0.5 p[2] * z_p_base_37[2] - 0.5 p[3] * z_p_base_37[2]
   + 0.5 p[4] * z_p_base_37[2] ] <= 0
 infeasibility_p_base_37_row_3_leq_0: 5e-06 z_p_base_37[3] + [
   p[2] * z_p_base_37[3] ] <= 0
 infeasibility_p_base_37_row_4_leq_0: 5e-06 z_p_base_37[4] + [
   - 0.5 p[0] * z_p_base_37[4] + 0.5 p[1] * z_p_base_37[4]
   - 0.5 p[2] * z_p_base_37[4] + 0.5 p[3] * z_p_base_37[4]
   + 0.5 p[4] * z_p_base_37[4] ] <= 0
 infeasibility_p_prime_base_37_row_0_leq_0: 5e-06 z_p_prime_base_37[0] + [
   0.5 alpha[0] * z_p_prime_base_37[0]
   - 0.5 alpha[1] * z_p_prime_base_37[0]
   + 0.5 alpha[2] * z_p_prime_base_37[0]
   + 0.5 alpha[3] * z_p_prime_base_37[0]
   - 0.5 alpha[4] * z_p_prime_base_37[0] + 0.5 p[0] * z_p_prime_base_37[0]
   + 0.5 p[1] * z_p_prime_base_37[0] - 0.5 p[2] * z_p_prime_base_37[0]
   - 0.5 p[3] * z_p_prime_base_37[0] - 0.5 p[4] * z_p_prime_base_37[0] ]
   <= 0
 infeasibility_p_prime_base_37_row_1_leq_0: 5e-06 z_p_prime_base_37[1] + [
   - alpha[3] * z_p_prime_base_37[1] + p[3] * z_p_prime_base_37[1] ] <= 0
 infeasibility_p_prime_base_37_row_2_leq_0: 5e-06 z_p_prime_base_37[2] + [
   0.5 alpha[0] * z_p_prime_base_37[2]
   + 0.5 alpha[1] * z_p_prime_base_37[2]
   - 0.5 alpha[2] * z_p_prime_base_37[2]
   + 0.5 alpha[3] * z_p_prime_base_37[2]
   + 0.5 alpha[4] * z_p_prime_base_37[2] + 0.5 p[0] * z_p_prime_base_37[2]
   - 0.5 p[1] * z_p_prime_base_37[2] + 0.5 p[2] * z_p_prime_base_37[2]
   - 0.5 p[3] * z_p_prime_base_37[2] + 0.5 p[4] * z_p_prime_base_37[2] ]
   <= 0
 infeasibility_p_prime_base_37_row_3_leq_0: 5e-06 z_p_prime_base_37[3] + [
   - alpha[2] * z_p_prime_base_37[3] + p[2] * z_p_prime_base_37[3] ] <= 0
 infeasibility_p_prime_base_37_row_4_leq_0: 5e-06 z_p_prime_base_37[4] + [
   - 0.5 alpha[0] * z_p_prime_base_37[4]
   - 0.5 alpha[1] * z_p_prime_base_37[4]
   + 0.5 alpha[2] * z_p_prime_base_37[4]
   - 0.5 alpha[3] * z_p_prime_base_37[4]
   + 0.5 alpha[4] * z_p_prime_base_37[4] - 0.5 p[0] * z_p_prime_base_37[4]
   + 0.5 p[1] * z_p_prime_base_37[4] - 0.5 p[2] * z_p_prime_base_37[4]
   + 0.5 p[3] * z_p_prime_base_37[4] + 0.5 p[4] * z_p_prime_base_37[4] ]
   <= 0
 infeasibility_p_base_38_row_0_leq_0: 5e-06 z_p_base_38[0] + [
   p[1] * z_p_base_38[0] - p[2] * z_p_base_38[0] ] <= 0
 infeasibility_p_base_38_row_1_leq_0: 5e-06 z_p_base_38[1] + [
   0.5 p[0] * z_p_base_38[1] - 0.5 p[1] * z_p_base_38[1]
   + 0.5 p[2] * z_p_base_38[1] + 0.5 p[3] * z_p_base_38[1]
   - 0.5 p[4] * z_p_base_38[1] ] <= 0
 infeasibility_p_base_38_row_2_leq_0: 5e-06 z_p_base_38[2] + [
   0.5 p[0] * z_p_base_38[2] - 0.5 p[1] * z_p_base_38[2]
   + 0.5 p[2] * z_p_base_38[2] - 0.5 p[3] * z_p_base_38[2]
   + 0.5 p[4] * z_p_base_38[2] ] <= 0
 infeasibility_p_base_38_row_3_leq_0: 5e-06 z_p_base_38[3] + [
   p[2] * z_p_base_38[3] ] <= 0
 infeasibility_p_base_38_row_4_leq_0: 5e-06 z_p_base_38[4] + [
   - 0.5 p[0] * z_p_base_38[4] + 0.5 p[1] * z_p_base_38[4]
   - 0.5 p[2] * z_p_base_38[4] + 0.5 p[3] * z_p_base_38[4]
   + 0.5 p[4] * z_p_base_38[4] ] <= 0
 infeasibility_p_prime_base_38_row_0_leq_0: 5e-06 z_p_prime_base_38[0] + [
   - alpha[1] * z_p_prime_base_38[0] + alpha[2] * z_p_prime_base_38[0]
   + p[1] * z_p_prime_base_38[0] - p[2] * z_p_prime_base_38[0] ] <= 0
 infeasibility_p_prime_base_38_row_1_leq_0: 5e-06 z_p_prime_base_38[1] + [
   0.5 alpha[0] * z_p_prime_base_38[1]
   + 0.5 alpha[1] * z_p_prime_base_38[1]
   - 0.5 alpha[2] * z_p_prime_base_38[1]
   - 0.5 alpha[3] * z_p_prime_base_38[1]
   - 0.5 alpha[4] * z_p_prime_base_38[1] + 0.5 p[0] * z_p_prime_base_38[1]
   - 0.5 p[1] * z_p_prime_base_38[1] + 0.5 p[2] * z_p_prime_base_38[1]
   + 0.5 p[3] * z_p_prime_base_38[1] - 0.5 p[4] * z_p_prime_base_38[1] ]
   <= 0
 infeasibility_p_prime_base_38_row_2_leq_0: 5e-06 z_p_prime_base_38[2] + [
   0.5 alpha[0] * z_p_prime_base_38[2]
   + 0.5 alpha[1] * z_p_prime_base_38[2]
   - 0.5 alpha[2] * z_p_prime_base_38[2]
   + 0.5 alpha[3] * z_p_prime_base_38[2]
   + 0.5 alpha[4] * z_p_prime_base_38[2] + 0.5 p[0] * z_p_prime_base_38[2]
   - 0.5 p[1] * z_p_prime_base_38[2] + 0.5 p[2] * z_p_prime_base_38[2]
   - 0.5 p[3] * z_p_prime_base_38[2] + 0.5 p[4] * z_p_prime_base_38[2] ]
   <= 0
 infeasibility_p_prime_base_38_row_3_leq_0: 5e-06 z_p_prime_base_38[3] + [
   - alpha[2] * z_p_prime_base_38[3] + p[2] * z_p_prime_base_38[3] ] <= 0
 infeasibility_p_prime_base_38_row_4_leq_0: 5e-06 z_p_prime_base_38[4] + [
   - 0.5 alpha[0] * z_p_prime_base_38[4]
   - 0.5 alpha[1] * z_p_prime_base_38[4]
   + 0.5 alpha[2] * z_p_prime_base_38[4]
   - 0.5 alpha[3] * z_p_prime_base_38[4]
   + 0.5 alpha[4] * z_p_prime_base_38[4] - 0.5 p[0] * z_p_prime_base_38[4]
   + 0.5 p[1] * z_p_prime_base_38[4] - 0.5 p[2] * z_p_prime_base_38[4]
   + 0.5 p[3] * z_p_prime_base_38[4] + 0.5 p[4] * z_p_prime_base_38[4] ]
   <= 0
 infeasibility_p_base_39_row_0_leq_0: 5e-06 z_p_base_39[0] + [
   0.5 p[0] * z_p_base_39[0] + 0.5 p[1] * z_p_base_39[0]
   + 0.5 p[2] * z_p_base_39[0] - 0.5 p[3] * z_p_base_39[0]
   - 0.5 p[4] * z_p_base_39[0] ] <= 0
 infeasibility_p_base_39_row_1_leq_0: 5e-06 z_p_base_39[1] + [
   0.5 p[0] * z_p_base_39[1] - 0.5 p[1] * z_p_base_39[1]
   - 0.5 p[2] * z_p_base_39[1] + 0.5 p[3] * z_p_base_39[1]
   - 0.5 p[4] * z_p_base_39[1] ] <= 0
 infeasibility_p_base_39_row_2_leq_0: 5e-06 z_p_base_39[2] + [
   p[4] * z_p_base_39[2] ] <= 0
 infeasibility_p_base_39_row_3_leq_0: 5e-06 z_p_base_39[3] + [
   - 0.5 p[0] * z_p_base_39[3] + 0.5 p[1] * z_p_base_39[3]
   - 0.5 p[2] * z_p_base_39[3] + 0.5 p[3] * z_p_base_39[3]
   + 0.5 p[4] * z_p_base_39[3] ] <= 0
 infeasibility_p_base_39_row_4_leq_0: 5e-06 z_p_base_39[4] + [
   p[2] * z_p_base_39[4] ] <= 0
 infeasibility_p_prime_base_39_row_0_leq_0: 5e-06 z_p_prime_base_39[0] + [
   0.5 alpha[0] * z_p_prime_base_39[0]
   - 0.5 alpha[1] * z_p_prime_base_39[0]
   - 0.5 alpha[2] * z_p_prime_base_39[0]
   + 0.5 alpha[3] * z_p_prime_base_39[0]
   - 0.5 alpha[4] * z_p_prime_base_39[0] + 0.5 p[0] * z_p_prime_base_39[0]
   + 0.5 p[1] * z_p_prime_base_39[0] + 0.5 p[2] * z_p_prime_base_39[0]
   - 0.5 p[3] * z_p_prime_base_39[0] - 0.5 p[4] * z_p_prime_base_39[0] ]
   <= 0
 infeasibility_p_prime_base_39_row_1_leq_0: 5e-06 z_p_prime_base_39[1] + [
   0.5 alpha[0] * z_p_prime_base_39[1]
   + 0.5 alpha[1] * z_p_prime_base_39[1]
   + 0.5 alpha[2] * z_p_prime_base_39[1]
   - 0.5 alpha[3] * z_p_prime_base_39[1]
   - 0.5 alpha[4] * z_p_prime_base_39[1] + 0.5 p[0] * z_p_prime_base_39[1]
   - 0.5 p[1] * z_p_prime_base_39[1] - 0.5 p[2] * z_p_prime_base_39[1]
   + 0.5 p[3] * z_p_prime_base_39[1] - 0.5 p[4] * z_p_prime_base_39[1] ]
   <= 0
 infeasibility_p_prime_base_39_row_2_leq_0: 5e-06 z_p_prime_base_39[2] + [
   alpha[4] * z_p_prime_base_39[2] + p[4] * z_p_prime_base_39[2] ] <= 0
 infeasibility_p_prime_base_39_row_3_leq_0: 5e-06 z_p_prime_base_39[3] + [
   - 0.5 alpha[0] * z_p_prime_base_39[3]
   - 0.5 alpha[1] * z_p_prime_base_39[3]
   + 0.5 alpha[2] * z_p_prime_base_39[3]
   - 0.5 alpha[3] * z_p_prime_base_39[3]
   + 0.5 alpha[4] * z_p_prime_base_39[3] - 0.5 p[0] * z_p_prime_base_39[3]
   + 0.5 p[1] * z_p_prime_base_39[3] - 0.5 p[2] * z_p_prime_base_39[3]
   + 0.5 p[3] * z_p_prime_base_39[3] + 0.5 p[4] * z_p_prime_base_39[3] ]
   <= 0
 infeasibility_p_prime_base_39_row_4_leq_0: 5e-06 z_p_prime_base_39[4] + [
   - alpha[2] * z_p_prime_base_39[4] + p[2] * z_p_prime_base_39[4] ] <= 0
 infeasibility_p_base_40_row_0_leq_0: 5e-06 z_p_base_40[0] + [
   0.5 p[0] * z_p_base_40[0] + 0.5 p[1] * z_p_base_40[0]
   + 0.5 p[2] * z_p_base_40[0] - 0.5 p[3] * z_p_base_40[0]
   - 0.5 p[4] * z_p_base_40[0] ] <= 0
 infeasibility_p_base_40_row_1_leq_0: 5e-06 z_p_base_40[1] + [
   0.5 p[0] * z_p_base_40[1] - 0.5 p[1] * z_p_base_40[1]
   + 0.5 p[2] * z_p_base_40[1] + 0.5 p[3] * z_p_base_40[1]
   - 0.5 p[4] * z_p_base_40[1] ] <= 0
 infeasibility_p_base_40_row_2_leq_0: 5e-06 z_p_base_40[2] + [
   - p[2] * z_p_base_40[2] + p[4] * z_p_base_40[2] ] <= 0
 infeasibility_p_base_40_row_3_leq_0: 5e-06 z_p_base_40[3] + [
   - 0.5 p[0] * z_p_base_40[3] + 0.5 p[1] * z_p_base_40[3]
   - 0.5 p[2] * z_p_base_40[3] + 0.5 p[3] * z_p_base_40[3]
   + 0.5 p[4] * z_p_base_40[3] ] <= 0
 infeasibility_p_base_40_row_4_leq_0: 5e-06 z_p_base_40[4] + [
   p[2] * z_p_base_40[4] ] <= 0
 infeasibility_p_prime_base_40_row_0_leq_0: 5e-06 z_p_prime_base_40[0] + [
   0.5 alpha[0] * z_p_prime_base_40[0]
   - 0.5 alpha[1] * z_p_prime_base_40[0]
   - 0.5 alpha[2] * z_p_prime_base_40[0]
   + 0.5 alpha[3] * z_p_prime_base_40[0]
   - 0.5 alpha[4] * z_p_prime_base_40[0] + 0.5 p[0] * z_p_prime_base_40[0]
   + 0.5 p[1] * z_p_prime_base_40[0] + 0.5 p[2] * z_p_prime_base_40[0]
   - 0.5 p[3] * z_p_prime_base_40[0] - 0.5 p[4] * z_p_prime_base_40[0] ]
   <= 0
 infeasibility_p_prime_base_40_row_1_leq_0: 5e-06 z_p_prime_base_40[1] + [
   0.5 alpha[0] * z_p_prime_base_40[1]
   + 0.5 alpha[1] * z_p_prime_base_40[1]
   - 0.5 alpha[2] * z_p_prime_base_40[1]
   - 0.5 alpha[3] * z_p_prime_base_40[1]
   - 0.5 alpha[4] * z_p_prime_base_40[1] + 0.5 p[0] * z_p_prime_base_40[1]
   - 0.5 p[1] * z_p_prime_base_40[1] + 0.5 p[2] * z_p_prime_base_40[1]
   + 0.5 p[3] * z_p_prime_base_40[1] - 0.5 p[4] * z_p_prime_base_40[1] ]
   <= 0
 infeasibility_p_prime_base_40_row_2_leq_0: 5e-06 z_p_prime_base_40[2] + [
   alpha[2] * z_p_prime_base_40[2] + alpha[4] * z_p_prime_base_40[2]
   - p[2] * z_p_prime_base_40[2] + p[4] * z_p_prime_base_40[2] ] <= 0
 infeasibility_p_prime_base_40_row_3_leq_0: 5e-06 z_p_prime_base_40[3] + [
   - 0.5 alpha[0] * z_p_prime_base_40[3]
   - 0.5 alpha[1] * z_p_prime_base_40[3]
   + 0.5 alpha[2] * z_p_prime_base_40[3]
   - 0.5 alpha[3] * z_p_prime_base_40[3]
   + 0.5 alpha[4] * z_p_prime_base_40[3] - 0.5 p[0] * z_p_prime_base_40[3]
   + 0.5 p[1] * z_p_prime_base_40[3] - 0.5 p[2] * z_p_prime_base_40[3]
   + 0.5 p[3] * z_p_prime_base_40[3] + 0.5 p[4] * z_p_prime_base_40[3] ]
   <= 0
 infeasibility_p_prime_base_40_row_4_leq_0: 5e-06 z_p_prime_base_40[4] + [
   - alpha[2] * z_p_prime_base_40[4] + p[2] * z_p_prime_base_40[4] ] <= 0
 infeasibility_p_base_41_row_0_leq_0: 5e-06 z_p_base_41[0] + [
   0.5 p[0] * z_p_base_41[0] + 0.5 p[1] * z_p_base_41[0]
   + 0.5 p[2] * z_p_base_41[0] - 0.5 p[3] * z_p_base_41[0]
   - 0.5 p[4] * z_p_base_41[0] ] <= 0
 infeasibility_p_base_41_row_1_leq_0: 5e-06 z_p_base_41[1] + [
   - p[2] * z_p_base_41[1] + p[3] * z_p_base_41[1] ] <= 0
 infeasibility_p_base_41_row_2_leq_0: 5e-06 z_p_base_41[2] + [
   0.5 p[0] * z_p_base_41[2] - 0.5 p[1] * z_p_base_41[2]
   + 0.5 p[2] * z_p_base_41[2] - 0.5 p[3] * z_p_base_41[2]
   + 0.5 p[4] * z_p_base_41[2] ] <= 0
 infeasibility_p_base_41_row_3_leq_0: 5e-06 z_p_base_41[3] + [
   - 0.5 p[0] * z_p_base_41[3] + 0.5 p[1] * z_p_base_41[3]
   - 0.5 p[2] * z_p_base_41[3] + 0.5 p[3] * z_p_base_41[3]
   + 0.5 p[4] * z_p_base_41[3] ] <= 0
 infeasibility_p_base_41_row_4_leq_0: 5e-06 z_p_base_41[4] + [
   p[2] * z_p_base_41[4] ] <= 0
 infeasibility_p_prime_base_41_row_0_leq_0: 5e-06 z_p_prime_base_41[0] + [
   0.5 alpha[0] * z_p_prime_base_41[0]
   - 0.5 alpha[1] * z_p_prime_base_41[0]
   - 0.5 alpha[2] * z_p_prime_base_41[0]
   + 0.5 alpha[3] * z_p_prime_base_41[0]
   - 0.5 alpha[4] * z_p_prime_base_41[0] + 0.5 p[0] * z_p_prime_base_41[0]
   + 0.5 p[1] * z_p_prime_base_41[0] + 0.5 p[2] * z_p_prime_base_41[0]
   - 0.5 p[3] * z_p_prime_base_41[0] - 0.5 p[4] * z_p_prime_base_41[0] ]
   <= 0
 infeasibility_p_prime_base_41_row_1_leq_0: 5e-06 z_p_prime_base_41[1] + [
   alpha[2] * z_p_prime_base_41[1] - alpha[3] * z_p_prime_base_41[1]
   - p[2] * z_p_prime_base_41[1] + p[3] * z_p_prime_base_41[1] ] <= 0
 infeasibility_p_prime_base_41_row_2_leq_0: 5e-06 z_p_prime_base_41[2] + [
   0.5 alpha[0] * z_p_prime_base_41[2]
   + 0.5 alpha[1] * z_p_prime_base_41[2]
   - 0.5 alpha[2] * z_p_prime_base_41[2]
   + 0.5 alpha[3] * z_p_prime_base_41[2]
   + 0.5 alpha[4] * z_p_prime_base_41[2] + 0.5 p[0] * z_p_prime_base_41[2]
   - 0.5 p[1] * z_p_prime_base_41[2] + 0.5 p[2] * z_p_prime_base_41[2]
   - 0.5 p[3] * z_p_prime_base_41[2] + 0.5 p[4] * z_p_prime_base_41[2] ]
   <= 0
 infeasibility_p_prime_base_41_row_3_leq_0: 5e-06 z_p_prime_base_41[3] + [
   - 0.5 alpha[0] * z_p_prime_base_41[3]
   - 0.5 alpha[1] * z_p_prime_base_41[3]
   + 0.5 alpha[2] * z_p_prime_base_41[3]
   - 0.5 alpha[3] * z_p_prime_base_41[3]
   + 0.5 alpha[4] * z_p_prime_base_41[3] - 0.5 p[0] * z_p_prime_base_41[3]
   + 0.5 p[1] * z_p_prime_base_41[3] - 0.5 p[2] * z_p_prime_base_41[3]
   + 0.5 p[3] * z_p_prime_base_41[3] + 0.5 p[4] * z_p_prime_base_41[3] ]
   <= 0
 infeasibility_p_prime_base_41_row_4_leq_0: 5e-06 z_p_prime_base_41[4] + [
   - alpha[2] * z_p_prime_base_41[4] + p[2] * z_p_prime_base_41[4] ] <= 0
 infeasibility_p_base_42_row_0_leq_0: 5e-06 z_p_base_42[0] + [
   0.5 p[0] * z_p_base_42[0] + 0.5 p[1] * z_p_base_42[0]
   + 0.5 p[2] * z_p_base_42[0] - 0.5 p[3] * z_p_base_42[0]
   - 0.5 p[4] * z_p_base_42[0] ] <= 0
 infeasibility_p_base_42_row_1_leq_0: 5e-06 z_p_base_42[1] + [
   p[3] * z_p_base_42[1] ] <= 0
 infeasibility_p_base_42_row_2_leq_0: 5e-06 z_p_base_42[2] + [
   0.5 p[0] * z_p_base_42[2] - 0.5 p[1] * z_p_base_42[2]
   - 0.5 p[2] * z_p_base_42[2] - 0.5 p[3] * z_p_base_42[2]
   + 0.5 p[4] * z_p_base_42[2] ] <= 0
 infeasibility_p_base_42_row_3_leq_0: 5e-06 z_p_base_42[3] + [
   - 0.5 p[0] * z_p_base_42[3] + 0.5 p[1] * z_p_base_42[3]
   - 0.5 p[2] * z_p_base_42[3] + 0.5 p[3] * z_p_base_42[3]
   + 0.5 p[4] * z_p_base_42[3] ] <= 0
 infeasibility_p_base_42_row_4_leq_0: 5e-06 z_p_base_42[4] + [
   p[2] * z_p_base_42[4] ] <= 0
 infeasibility_p_prime_base_42_row_0_leq_0: 5e-06 z_p_prime_base_42[0] + [
   0.5 alpha[0] * z_p_prime_base_42[0]
   - 0.5 alpha[1] * z_p_prime_base_42[0]
   - 0.5 alpha[2] * z_p_prime_base_42[0]
   + 0.5 alpha[3] * z_p_prime_base_42[0]
   - 0.5 alpha[4] * z_p_prime_base_42[0] + 0.5 p[0] * z_p_prime_base_42[0]
   + 0.5 p[1] * z_p_prime_base_42[0] + 0.5 p[2] * z_p_prime_base_42[0]
   - 0.5 p[3] * z_p_prime_base_42[0] - 0.5 p[4] * z_p_prime_base_42[0] ]
   <= 0
 infeasibility_p_prime_base_42_row_1_leq_0: 5e-06 z_p_prime_base_42[1] + [
   - alpha[3] * z_p_prime_base_42[1] + p[3] * z_p_prime_base_42[1] ] <= 0
 infeasibility_p_prime_base_42_row_2_leq_0: 5e-06 z_p_prime_base_42[2] + [
   0.5 alpha[0] * z_p_prime_base_42[2]
   + 0.5 alpha[1] * z_p_prime_base_42[2]
   + 0.5 alpha[2] * z_p_prime_base_42[2]
   + 0.5 alpha[3] * z_p_prime_base_42[2]
   + 0.5 alpha[4] * z_p_prime_base_42[2] + 0.5 p[0] * z_p_prime_base_42[2]
   - 0.5 p[1] * z_p_prime_base_42[2] - 0.5 p[2] * z_p_prime_base_42[2]
   - 0.5 p[3] * z_p_prime_base_42[2] + 0.5 p[4] * z_p_prime_base_42[2] ]
   <= 0
 infeasibility_p_prime_base_42_row_3_leq_0: 5e-06 z_p_prime_base_42[3] + [
   - 0.5 alpha[0] * z_p_prime_base_42[3]
   - 0.5 alpha[1] * z_p_prime_base_42[3]
   + 0.5 alpha[2] * z_p_prime_base_42[3]
   - 0.5 alpha[3] * z_p_prime_base_42[3]
   + 0.5 alpha[4] * z_p_prime_base_42[3] - 0.5 p[0] * z_p_prime_base_42[3]
   + 0.5 p[1] * z_p_prime_base_42[3] - 0.5 p[2] * z_p_prime_base_42[3]
   + 0.5 p[3] * z_p_prime_base_42[3] + 0.5 p[4] * z_p_prime_base_42[3] ]
   <= 0
 infeasibility_p_prime_base_42_row_4_leq_0: 5e-06 z_p_prime_base_42[4] + [
   - alpha[2] * z_p_prime_base_42[4] + p[2] * z_p_prime_base_42[4] ] <= 0
 infeasibility_p_base_43_row_0_leq_0: 5e-06 z_p_base_43[0] + [
   p[1] * z_p_base_43[0] ] <= 0
 infeasibility_p_base_43_row_1_leq_0: 5e-06 z_p_base_43[1] + [
   0.5 p[0] * z_p_base_43[1] - 0.5 p[1] * z_p_base_43[1]
   - 0.5 p[2] * z_p_base_43[1] + 0.5 p[3] * z_p_base_43[1]
   - 0.5 p[4] * z_p_base_43[1] ] <= 0
 infeasibility_p_base_43_row_2_leq_0: 5e-06 z_p_base_43[2] + [
   0.5 p[0] * z_p_base_43[2] - 0.5 p[1] * z_p_base_43[2]
   + 0.5 p[2] * z_p_base_43[2] - 0.5 p[3] * z_p_base_43[2]
   + 0.5 p[4] * z_p_base_43[2] ] <= 0
 infeasibility_p_base_43_row_3_leq_0: 5e-06 z_p_base_43[3] + [
   p[2] * z_p_base_43[3] ] <= 0
 infeasibility_p_base_43_row_4_leq_0: 5e-06 z_p_base_43[4] + [
   - 0.5 p[0] * z_p_base_43[4] + 0.5 p[1] * z_p_base_43[4]
   - 0.5 p[2] * z_p_base_43[4] + 0.5 p[3] * z_p_base_43[4]
   + 0.5 p[4] * z_p_base_43[4] ] <= 0
 infeasibility_p_prime_base_43_row_0_leq_0: 5e-06 z_p_prime_base_43[0] + [
   - alpha[1] * z_p_prime_base_43[0] + p[1] * z_p_prime_base_43[0] ] <= 0
 infeasibility_p_prime_base_43_row_1_leq_0: 5e-06 z_p_prime_base_43[1] + [
   0.5 alpha[0] * z_p_prime_base_43[1]
   + 0.5 alpha[1] * z_p_prime_base_43[1]
   + 0.5 alpha[2] * z_p_prime_base_43[1]
   - 0.5 alpha[3] * z_p_prime_base_43[1]
   - 0.5 alpha[4] * z_p_prime_base_43[1] + 0.5 p[0] * z_p_prime_base_43[1]
   - 0.5 p[1] * z_p_prime_base_43[1] - 0.5 p[2] * z_p_prime_base_43[1]
   + 0.5 p[3] * z_p_prime_base_43[1] - 0.5 p[4] * z_p_prime_base_43[1] ]
   <= 0
 infeasibility_p_prime_base_43_row_2_leq_0: 5e-06 z_p_prime_base_43[2] + [
   0.5 alpha[0] * z_p_prime_base_43[2]
   + 0.5 alpha[1] * z_p_prime_base_43[2]
   - 0.5 alpha[2] * z_p_prime_base_43[2]
   + 0.5 alpha[3] * z_p_prime_base_43[2]
   + 0.5 alpha[4] * z_p_prime_base_43[2] + 0.5 p[0] * z_p_prime_base_43[2]
   - 0.5 p[1] * z_p_prime_base_43[2] + 0.5 p[2] * z_p_prime_base_43[2]
   - 0.5 p[3] * z_p_prime_base_43[2] + 0.5 p[4] * z_p_prime_base_43[2] ]
   <= 0
 infeasibility_p_prime_base_43_row_3_leq_0: 5e-06 z_p_prime_base_43[3] + [
   - alpha[2] * z_p_prime_base_43[3] + p[2] * z_p_prime_base_43[3] ] <= 0
 infeasibility_p_prime_base_43_row_4_leq_0: 5e-06 z_p_prime_base_43[4] + [
   - 0.5 alpha[0] * z_p_prime_base_43[4]
   - 0.5 alpha[1] * z_p_prime_base_43[4]
   + 0.5 alpha[2] * z_p_prime_base_43[4]
   - 0.5 alpha[3] * z_p_prime_base_43[4]
   + 0.5 alpha[4] * z_p_prime_base_43[4] - 0.5 p[0] * z_p_prime_base_43[4]
   + 0.5 p[1] * z_p_prime_base_43[4] - 0.5 p[2] * z_p_prime_base_43[4]
   + 0.5 p[3] * z_p_prime_base_43[4] + 0.5 p[4] * z_p_prime_base_43[4] ]
   <= 0
 infeasibility_p_base_44_row_0_leq_0: 5e-06 z_p_base_44[0] + [
   p[1] * z_p_base_44[0] ] <= 0
 infeasibility_p_base_44_row_1_leq_0: 5e-06 z_p_base_44[1] + [
   0.5 p[0] * z_p_base_44[1] - 0.5 p[1] * z_p_base_44[1]
   + 0.5 p[2] * z_p_base_44[1] + 0.5 p[3] * z_p_base_44[1]
   - 0.5 p[4] * z_p_base_44[1] ] <= 0
 infeasibility_p_base_44_row_2_leq_0: 5e-06 z_p_base_44[2] + [
   0.5 p[0] * z_p_base_44[2] - 0.5 p[1] * z_p_base_44[2]
   - 0.5 p[2] * z_p_base_44[2] - 0.5 p[3] * z_p_base_44[2]
   + 0.5 p[4] * z_p_base_44[2] ] <= 0
 infeasibility_p_base_44_row_3_leq_0: 5e-06 z_p_base_44[3] + [
   p[2] * z_p_base_44[3] ] <= 0
 infeasibility_p_base_44_row_4_leq_0: 5e-06 z_p_base_44[4] + [
   - 0.5 p[0] * z_p_base_44[4] + 0.5 p[1] * z_p_base_44[4]
   - 0.5 p[2] * z_p_base_44[4] + 0.5 p[3] * z_p_base_44[4]
   + 0.5 p[4] * z_p_base_44[4] ] <= 0
 infeasibility_p_prime_base_44_row_0_leq_0: 5e-06 z_p_prime_base_44[0] + [
   - alpha[1] * z_p_prime_base_44[0] + p[1] * z_p_prime_base_44[0] ] <= 0
 infeasibility_p_prime_base_44_row_1_leq_0: 5e-06 z_p_prime_base_44[1] + [
   0.5 alpha[0] * z_p_prime_base_44[1]
   + 0.5 alpha[1] * z_p_prime_base_44[1]
   - 0.5 alpha[2] * z_p_prime_base_44[1]
   - 0.5 alpha[3] * z_p_prime_base_44[1]
   - 0.5 alpha[4] * z_p_prime_base_44[1] + 0.5 p[0] * z_p_prime_base_44[1]
   - 0.5 p[1] * z_p_prime_base_44[1] + 0.5 p[2] * z_p_prime_base_44[1]
   + 0.5 p[3] * z_p_prime_base_44[1] - 0.5 p[4] * z_p_prime_base_44[1] ]
   <= 0
 infeasibility_p_prime_base_44_row_2_leq_0: 5e-06 z_p_prime_base_44[2] + [
   0.5 alpha[0] * z_p_prime_base_44[2]
   + 0.5 alpha[1] * z_p_prime_base_44[2]
   + 0.5 alpha[2] * z_p_prime_base_44[2]
   + 0.5 alpha[3] * z_p_prime_base_44[2]
   + 0.5 alpha[4] * z_p_prime_base_44[2] + 0.5 p[0] * z_p_prime_base_44[2]
   - 0.5 p[1] * z_p_prime_base_44[2] - 0.5 p[2] * z_p_prime_base_44[2]
   - 0.5 p[3] * z_p_prime_base_44[2] + 0.5 p[4] * z_p_prime_base_44[2] ]
   <= 0
 infeasibility_p_prime_base_44_row_3_leq_0: 5e-06 z_p_prime_base_44[3] + [
   - alpha[2] * z_p_prime_base_44[3] + p[2] * z_p_prime_base_44[3] ] <= 0
 infeasibility_p_prime_base_44_row_4_leq_0: 5e-06 z_p_prime_base_44[4] + [
   - 0.5 alpha[0] * z_p_prime_base_44[4]
   - 0.5 alpha[1] * z_p_prime_base_44[4]
   + 0.5 alpha[2] * z_p_prime_base_44[4]
   - 0.5 alpha[3] * z_p_prime_base_44[4]
   + 0.5 alpha[4] * z_p_prime_base_44[4] - 0.5 p[0] * z_p_prime_base_44[4]
   + 0.5 p[1] * z_p_prime_base_44[4] - 0.5 p[2] * z_p_prime_base_44[4]
   + 0.5 p[3] * z_p_prime_base_44[4] + 0.5 p[4] * z_p_prime_base_44[4] ]
   <= 0
 infeasibility_p_base_45_row_0_leq_0: 5e-06 z_p_base_45[0] + [
   0.5 p[0] * z_p_base_45[0] + 0.5 p[1] * z_p_base_45[0]
   - 0.5 p[2] * z_p_base_45[0] - 0.5 p[3] * z_p_base_45[0]
   - 0.5 p[4] * z_p_base_45[0] ] <= 0
 infeasibility_p_base_45_row_1_leq_0: 5e-06 z_p_base_45[1] + [
   0.5 p[0] * z_p_base_45[1] - 0.5 p[1] * z_p_base_45[1]
   + 0.5 p[2] * z_p_base_45[1] + 0.5 p[3] * z_p_base_45[1]
   + 0.5 p[4] * z_p_base_45[1] ] <= 0
 infeasibility_p_base_45_row_2_leq_0: 5e-06 z_p_base_45[2] + [
   p[2] * z_p_base_45[2] ] <= 0
 infeasibility_p_base_45_row_3_leq_0: 5e-06 z_p_base_45[3] + [
   - 0.5 p[0] * z_p_base_45[3] + 0.5 p[1] * z_p_base_45[3]
   - 0.5 p[2] * z_p_base_45[3] + 0.5 p[3] * z_p_base_45[3]
   - 0.5 p[4] * z_p_base_45[3] ] <= 0
 infeasibility_p_base_45_row_4_leq_0: 5e-06 z_p_base_45[4] + [
   p[4] * z_p_base_45[4] ] <= 0
 infeasibility_p_prime_base_45_row_0_leq_0: 5e-06 z_p_prime_base_45[0] + [
   0.5 alpha[0] * z_p_prime_base_45[0]
   - 0.5 alpha[1] * z_p_prime_base_45[0]
   + 0.5 alpha[2] * z_p_prime_base_45[0]
   + 0.5 alpha[3] * z_p_prime_base_45[0]
   - 0.5 alpha[4] * z_p_prime_base_45[0] + 0.5 p[0] * z_p_prime_base_45[0]
   + 0.5 p[1] * z_p_prime_base_45[0] - 0.5 p[2] * z_p_prime_base_45[0]
   - 0.5 p[3] * z_p_prime_base_45[0] - 0.5 p[4] * z_p_prime_base_45[0] ]
   <= 0
 infeasibility_p_prime_base_45_row_1_leq_0: 5e-06 z_p_prime_base_45[1] + [
   0.5 alpha[0] * z_p_prime_base_45[1]
   + 0.5 alpha[1] * z_p_prime_base_45[1]
   - 0.5 alpha[2] * z_p_prime_base_45[1]
   - 0.5 alpha[3] * z_p_prime_base_45[1]
   + 0.5 alpha[4] * z_p_prime_base_45[1] + 0.5 p[0] * z_p_prime_base_45[1]
   - 0.5 p[1] * z_p_prime_base_45[1] + 0.5 p[2] * z_p_prime_base_45[1]
   + 0.5 p[3] * z_p_prime_base_45[1] + 0.5 p[4] * z_p_prime_base_45[1] ]
   <= 0
 infeasibility_p_prime_base_45_row_2_leq_0: 5e-06 z_p_prime_base_45[2] + [
   - alpha[2] * z_p_prime_base_45[2] + p[2] * z_p_prime_base_45[2] ] <= 0
 infeasibility_p_prime_base_45_row_3_leq_0: 5e-06 z_p_prime_base_45[3] + [
   - 0.5 alpha[0] * z_p_prime_base_45[3]
   - 0.5 alpha[1] * z_p_prime_base_45[3]
   + 0.5 alpha[2] * z_p_prime_base_45[3]
   - 0.5 alpha[3] * z_p_prime_base_45[3]
   - 0.5 alpha[4] * z_p_prime_base_45[3] - 0.5 p[0] * z_p_prime_base_45[3]
   + 0.5 p[1] * z_p_prime_base_45[3] - 0.5 p[2] * z_p_prime_base_45[3]
   + 0.5 p[3] * z_p_prime_base_45[3] - 0.5 p[4] * z_p_prime_base_45[3] ]
   <= 0
 infeasibility_p_prime_base_45_row_4_leq_0: 5e-06 z_p_prime_base_45[4] + [
   alpha[4] * z_p_prime_base_45[4] + p[4] * z_p_prime_base_45[4] ] <= 0
 infeasibility_p_base_46_row_0_leq_0: 5e-06 z_p_base_46[0] + [
   0.5 p[0] * z_p_base_46[0] + 0.5 p[1] * z_p_base_46[0]
   - 0.5 p[2] * z_p_base_46[0] - 0.5 p[3] * z_p_base_46[0]
   + 0.5 p[4] * z_p_base_46[0] ] <= 0
 infeasibility_p_base_46_row_1_leq_0: 5e-06 z_p_base_46[1] + [
   0.5 p[0] * z_p_base_46[1] - 0.5 p[1] * z_p_base_46[1]
   + 0.5 p[2] * z_p_base_46[1] + 0.5 p[3] * z_p_base_46[1]
   - 0.5 p[4] * z_p_base_46[1] ] <= 0
 infeasibility_p_base_46_row_2_leq_0: 5e-06 z_p_base_46[2] + [
   p[2] * z_p_base_46[2] - p[4] * z_p_base_46[2] ] <= 0
 infeasibility_p_base_46_row_3_leq_0: 5e-06 z_p_base_46[3] + [
   - 0.5 p[0] * z_p_base_46[3] + 0.5 p[1] * z_p_base_46[3]
   - 0.5 p[2] * z_p_base_46[3] + 0.5 p[3] * z_p_base_46[3]
   + 0.5 p[4] * z_p_base_46[3] ] <= 0
 infeasibility_p_base_46_row_4_leq_0: 5e-06 z_p_base_46[4] + [
   p[4] * z_p_base_46[4] ] <= 0
 infeasibility_p_prime_base_46_row_0_leq_0: 5e-06 z_p_prime_base_46[0] + [
   0.5 alpha[0] * z_p_prime_base_46[0]
   - 0.5 alpha[1] * z_p_prime_base_46[0]
   + 0.5 alpha[2] * z_p_prime_base_46[0]
   + 0.5 alpha[3] * z_p_prime_base_46[0]
   + 0.5 alpha[4] * z_p_prime_base_46[0] + 0.5 p[0] * z_p_prime_base_46[0]
   + 0.5 p[1] * z_p_prime_base_46[0] - 0.5 p[2] * z_p_prime_base_46[0]
   - 0.5 p[3] * z_p_prime_base_46[0] + 0.5 p[4] * z_p_prime_base_46[0] ]
   <= 0
 infeasibility_p_prime_base_46_row_1_leq_0: 5e-06 z_p_prime_base_46[1] + [
   0.5 alpha[0] * z_p_prime_base_46[1]
   + 0.5 alpha[1] * z_p_prime_base_46[1]
   - 0.5 alpha[2] * z_p_prime_base_46[1]
   - 0.5 alpha[3] * z_p_prime_base_46[1]
   - 0.5 alpha[4] * z_p_prime_base_46[1] + 0.5 p[0] * z_p_prime_base_46[1]
   - 0.5 p[1] * z_p_prime_base_46[1] + 0.5 p[2] * z_p_prime_base_46[1]
   + 0.5 p[3] * z_p_prime_base_46[1] - 0.5 p[4] * z_p_prime_base_46[1] ]
   <= 0
 infeasibility_p_prime_base_46_row_2_leq_0: 5e-06 z_p_prime_base_46[2] + [
   - alpha[2] * z_p_prime_base_46[2] - alpha[4] * z_p_prime_base_46[2]
   + p[2] * z_p_prime_base_46[2] - p[4] * z_p_prime_base_46[2] ] <= 0
 infeasibility_p_prime_base_46_row_3_leq_0: 5e-06 z_p_prime_base_46[3] + [
   - 0.5 alpha[0] * z_p_prime_base_46[3]
   - 0.5 alpha[1] * z_p_prime_base_46[3]
   + 0.5 alpha[2] * z_p_prime_base_46[3]
   - 0.5 alpha[3] * z_p_prime_base_46[3]
   + 0.5 alpha[4] * z_p_prime_base_46[3] - 0.5 p[0] * z_p_prime_base_46[3]
   + 0.5 p[1] * z_p_prime_base_46[3] - 0.5 p[2] * z_p_prime_base_46[3]
   + 0.5 p[3] * z_p_prime_base_46[3] + 0.5 p[4] * z_p_prime_base_46[3] ]
   <= 0
 infeasibility_p_prime_base_46_row_4_leq_0: 5e-06 z_p_prime_base_46[4] + [
   alpha[4] * z_p_prime_base_46[4] + p[4] * z_p_prime_base_46[4] ] <= 0
 infeasibility_p_base_47_row_0_leq_0: 5e-06 z_p_base_47[0] + [
   0.5 p[0] * z_p_base_47[0] + 0.5 p[1] * z_p_base_47[0]
   - 0.5 p[2] * z_p_base_47[0] - 0.5 p[3] * z_p_base_47[0]
   + 0.5 p[4] * z_p_base_47[0] ] <= 0
 infeasibility_p_base_47_row_1_leq_0: 5e-06 z_p_base_47[1] + [
   0.5 p[0] * z_p_base_47[1] - 0.5 p[1] * z_p_base_47[1]
   + 0.5 p[2] * z_p_base_47[1] + 0.5 p[3] * z_p_base_47[1]
   - 0.5 p[4] * z_p_base_47[1] ] <= 0
 infeasibility_p_base_47_row_2_leq_0: 5e-06 z_p_base_47[2] + [
   p[2] * z_p_base_47[2] ] <= 0
 infeasibility_p_base_47_row_3_leq_0: 5e-06 z_p_base_47[3] + [
   - 0.5 p[0] * z_p_base_47[3] + 0.5 p[1] * z_p_base_47[3]
   - 0.5 p[2] * z_p_base_47[3] + 0.5 p[3] * z_p_base_47[3]
   - 0.5 p[4] * z_p_base_47[3] ] <= 0
 infeasibility_p_base_47_row_4_leq_0: 5e-06 z_p_base_47[4] + [
   p[4] * z_p_base_47[4] ] <= 0
 infeasibility_p_prime_base_47_row_0_leq_0: 5e-06 z_p_prime_base_47[0] + [
   0.5 alpha[0] * z_p_prime_base_47[0]
   - 0.5 alpha[1] * z_p_prime_base_47[0]
   + 0.5 alpha[2] * z_p_prime_base_47[0]
   + 0.5 alpha[3] * z_p_prime_base_47[0]
   + 0.5 alpha[4] * z_p_prime_base_47[0] + 0.5 p[0] * z_p_prime_base_47[0]
   + 0.5 p[1] * z_p_prime_base_47[0] - 0.5 p[2] * z_p_prime_base_47[0]
   - 0.5 p[3] * z_p_prime_base_47[0] + 0.5 p[4] * z_p_prime_base_47[0] ]
   <= 0
 infeasibility_p_prime_base_47_row_1_leq_0: 5e-06 z_p_prime_base_47[1] + [
   0.5 alpha[0] * z_p_prime_base_47[1]
   + 0.5 alpha[1] * z_p_prime_base_47[1]
   - 0.5 alpha[2] * z_p_prime_base_47[1]
   - 0.5 alpha[3] * z_p_prime_base_47[1]
   - 0.5 alpha[4] * z_p_prime_base_47[1] + 0.5 p[0] * z_p_prime_base_47[1]
   - 0.5 p[1] * z_p_prime_base_47[1] + 0.5 p[2] * z_p_prime_base_47[1]
   + 0.5 p[3] * z_p_prime_base_47[1] - 0.5 p[4] * z_p_prime_base_47[1] ]
   <= 0
 infeasibility_p_prime_base_47_row_2_leq_0: 5e-06 z_p_prime_base_47[2] + [
   - alpha[2] * z_p_prime_base_47[2] + p[2] * z_p_prime_base_47[2] ] <= 0
 infeasibility_p_prime_base_47_row_3_leq_0: 5e-06 z_p_prime_base_47[3] + [
   - 0.5 alpha[0] * z_p_prime_base_47[3]
   - 0.5 alpha[1] * z_p_prime_base_47[3]
   + 0.5 alpha[2] * z_p_prime_base_47[3]
   - 0.5 alpha[3] * z_p_prime_base_47[3]
   - 0.5 alpha[4] * z_p_prime_base_47[3] - 0.5 p[0] * z_p_prime_base_47[3]
   + 0.5 p[1] * z_p_prime_base_47[3] - 0.5 p[2] * z_p_prime_base_47[3]
   + 0.5 p[3] * z_p_prime_base_47[3] - 0.5 p[4] * z_p_prime_base_47[3] ]
   <= 0
 infeasibility_p_prime_base_47_row_4_leq_0: 5e-06 z_p_prime_base_47[4] + [
   alpha[4] * z_p_prime_base_47[4] + p[4] * z_p_prime_base_47[4] ] <= 0
 infeasibility_p_base_48_row_0_leq_0: 5e-06 z_p_base_48[0] + [
   p[0] * z_p_base_48[0] - p[3] * z_p_base_48[0] ] <= 0
 infeasibility_p_base_48_row_1_leq_0: 5e-06 z_p_base_48[1] + [
   p[3] * z_p_base_48[1] ] <= 0
 infeasibility_p_base_48_row_2_leq_0: 5e-06 z_p_base_48[2] + [
   - 0.5 p[0] * z_p_base_48[2] + 0.5 p[1] * z_p_base_48[2]
   + 0.5 p[2] * z_p_base_48[2] + 0.5 p[3] * z_p_base_48[2]
   - 0.5 p[4] * z_p_base_48[2] ] <= 0
 infeasibility_p_base_48_row_3_leq_0: 5e-06 z_p_base_48[3] + [
   - 0.5 p[0] * z_p_base_48[3] + 0.5 p[1] * z_p_base_48[3]
   - 0.5 p[2] * z_p_base_48[3] + 0.5 p[3] * z_p_base_48[3]
   + 0.5 p[4] * z_p_base_48[3] ] <= 0
 infeasibility_p_base_48_row_4_leq_0: 5e-06 z_p_base_48[4] + [
   0.5 p[0] * z_p_base_48[4] - 0.5 p[1] * z_p_base_48[4]
   + 0.5 p[2] * z_p_base_48[4] - 0.5 p[3] * z_p_base_48[4]
   + 0.5 p[4] * z_p_base_48[4] ] <= 0
 infeasibility_p_prime_base_48_row_0_leq_0: 5e-06 z_p_prime_base_48[0] + [
   alpha[0] * z_p_prime_base_48[0] + alpha[3] * z_p_prime_base_48[0]
   + p[0] * z_p_prime_base_48[0] - p[3] * z_p_prime_base_48[0] ] <= 0
 infeasibility_p_prime_base_48_row_1_leq_0: 5e-06 z_p_prime_base_48[1] + [
   - alpha[3] * z_p_prime_base_48[1] + p[3] * z_p_prime_base_48[1] ] <= 0
 infeasibility_p_prime_base_48_row_2_leq_0: 5e-06 z_p_prime_base_48[2] + [
   - 0.5 alpha[0] * z_p_prime_base_48[2]
   - 0.5 alpha[1] * z_p_prime_base_48[2]
   - 0.5 alpha[2] * z_p_prime_base_48[2]
   - 0.5 alpha[3] * z_p_prime_base_48[2]
   - 0.5 alpha[4] * z_p_prime_base_48[2] - 0.5 p[0] * z_p_prime_base_48[2]
   + 0.5 p[1] * z_p_prime_base_48[2] + 0.5 p[2] * z_p_prime_base_48[2]
   + 0.5 p[3] * z_p_prime_base_48[2] - 0.5 p[4] * z_p_prime_base_48[2] ]
   <= 0
 infeasibility_p_prime_base_48_row_3_leq_0: 5e-06 z_p_prime_base_48[3] + [
   - 0.5 alpha[0] * z_p_prime_base_48[3]
   - 0.5 alpha[1] * z_p_prime_base_48[3]
   + 0.5 alpha[2] * z_p_prime_base_48[3]
   - 0.5 alpha[3] * z_p_prime_base_48[3]
   + 0.5 alpha[4] * z_p_prime_base_48[3] - 0.5 p[0] * z_p_prime_base_48[3]
   + 0.5 p[1] * z_p_prime_base_48[3] - 0.5 p[2] * z_p_prime_base_48[3]
   + 0.5 p[3] * z_p_prime_base_48[3] + 0.5 p[4] * z_p_prime_base_48[3] ]
   <= 0
 infeasibility_p_prime_base_48_row_4_leq_0: 5e-06 z_p_prime_base_48[4] + [
   0.5 alpha[0] * z_p_prime_base_48[4]
   + 0.5 alpha[1] * z_p_prime_base_48[4]
   - 0.5 alpha[2] * z_p_prime_base_48[4]
   + 0.5 alpha[3] * z_p_prime_base_48[4]
   + 0.5 alpha[4] * z_p_prime_base_48[4] + 0.5 p[0] * z_p_prime_base_48[4]
   - 0.5 p[1] * z_p_prime_base_48[4] + 0.5 p[2] * z_p_prime_base_48[4]
   - 0.5 p[3] * z_p_prime_base_48[4] + 0.5 p[4] * z_p_prime_base_48[4] ]
   <= 0
 infeasibility_p_base_49_row_0_leq_0: 5e-06 z_p_base_49[0] + [
   0.5 p[0] * z_p_base_49[0] + 0.5 p[1] * z_p_base_49[0]
   - 0.5 p[2] * z_p_base_49[0] - 0.5 p[3] * z_p_base_49[0]
   + 0.5 p[4] * z_p_base_49[0] ] <= 0
 infeasibility_p_base_49_row_1_leq_0: 5e-06 z_p_base_49[1] + [
   0.5 p[0] * z_p_base_49[1] - 0.5 p[1] * z_p_base_49[1]
   + 0.5 p[2] * z_p_base_49[1] + 0.5 p[3] * z_p_base_49[1]
   - 0.5 p[4] * z_p_base_49[1] ] <= 0
 infeasibility_p_base_49_row_2_leq_0: 5e-06 z_p_base_49[2] + [
   - 0.5 p[0] * z_p_base_49[2] + 0.5 p[1] * z_p_base_49[2]
   + 0.5 p[2] * z_p_base_49[2] + 0.5 p[3] * z_p_base_49[2]
   - 0.5 p[4] * z_p_base_49[2] ] <= 0
 infeasibility_p_base_49_row_3_leq_0: 5e-06 z_p_base_49[3] + [
   0.5 p[0] * z_p_base_49[3] - 0.5 p[1] * z_p_base_49[3]
   + 0.5 p[2] * z_p_base_49[3] - 0.5 p[3] * z_p_base_49[3]
   + 0.5 p[4] * z_p_base_49[3] ] <= 0
 infeasibility_p_base_49_row_4_leq_0: 5e-06 z_p_base_49[4] + [
   - 0.5 p[0] * z_p_base_49[4] + 0.5 p[1] * z_p_base_49[4]
   - 0.5 p[2] * z_p_base_49[4] + 0.5 p[3] * z_p_base_49[4]
   + 0.5 p[4] * z_p_base_49[4] ] <= 0
 infeasibility_p_prime_base_49_row_0_leq_0: 5e-06 z_p_prime_base_49[0] + [
   0.5 alpha[0] * z_p_prime_base_49[0]
   - 0.5 alpha[1] * z_p_prime_base_49[0]
   + 0.5 alpha[2] * z_p_prime_base_49[0]
   + 0.5 alpha[3] * z_p_prime_base_49[0]
   + 0.5 alpha[4] * z_p_prime_base_49[0] + 0.5 p[0] * z_p_prime_base_49[0]
   + 0.5 p[1] * z_p_prime_base_49[0] - 0.5 p[2] * z_p_prime_base_49[0]
   - 0.5 p[3] * z_p_prime_base_49[0] + 0.5 p[4] * z_p_prime_base_49[0] ]
   <= 0
 infeasibility_p_prime_base_49_row_1_leq_0: 5e-06 z_p_prime_base_49[1] + [
   0.5 alpha[0] * z_p_prime_base_49[1]
   + 0.5 alpha[1] * z_p_prime_base_49[1]
   - 0.5 alpha[2] * z_p_prime_base_49[1]
   - 0.5 alpha[3] * z_p_prime_base_49[1]
   - 0.5 alpha[4] * z_p_prime_base_49[1] + 0.5 p[0] * z_p_prime_base_49[1]
   - 0.5 p[1] * z_p_prime_base_49[1] + 0.5 p[2] * z_p_prime_base_49[1]
   + 0.5 p[3] * z_p_prime_base_49[1] - 0.5 p[4] * z_p_prime_base_49[1] ]
   <= 0
 infeasibility_p_prime_base_49_row_2_leq_0: 5e-06 z_p_prime_base_49[2] + [
   - 0.5 alpha[0] * z_p_prime_base_49[2]
   - 0.5 alpha[1] * z_p_prime_base_49[2]
   - 0.5 alpha[2] * z_p_prime_base_49[2]
   - 0.5 alpha[3] * z_p_prime_base_49[2]
   - 0.5 alpha[4] * z_p_prime_base_49[2] - 0.5 p[0] * z_p_prime_base_49[2]
   + 0.5 p[1] * z_p_prime_base_49[2] + 0.5 p[2] * z_p_prime_base_49[2]
   + 0.5 p[3] * z_p_prime_base_49[2] - 0.5 p[4] * z_p_prime_base_49[2] ]
   <= 0
 infeasibility_p_prime_base_49_row_3_leq_0: 5e-06 z_p_prime_base_49[3] + [
   0.5 alpha[0] * z_p_prime_base_49[3]
   + 0.5 alpha[1] * z_p_prime_base_49[3]
   - 0.5 alpha[2] * z_p_prime_base_49[3]
   + 0.5 alpha[3] * z_p_prime_base_49[3]
   + 0.5 alpha[4] * z_p_prime_base_49[3] + 0.5 p[0] * z_p_prime_base_49[3]
   - 0.5 p[1] * z_p_prime_base_49[3] + 0.5 p[2] * z_p_prime_base_49[3]
   - 0.5 p[3] * z_p_prime_base_49[3] + 0.5 p[4] * z_p_prime_base_49[3] ]
   <= 0
 infeasibility_p_prime_base_49_row_4_leq_0: 5e-06 z_p_prime_base_49[4] + [
   - 0.5 alpha[0] * z_p_prime_base_49[4]
   - 0.5 alpha[1] * z_p_prime_base_49[4]
   + 0.5 alpha[2] * z_p_prime_base_49[4]
   - 0.5 alpha[3] * z_p_prime_base_49[4]
   + 0.5 alpha[4] * z_p_prime_base_49[4] - 0.5 p[0] * z_p_prime_base_49[4]
   + 0.5 p[1] * z_p_prime_base_49[4] - 0.5 p[2] * z_p_prime_base_49[4]
   + 0.5 p[3] * z_p_prime_base_49[4] + 0.5 p[4] * z_p_prime_base_49[4] ]
   <= 0
 infeasibility_p_base_50_row_0_leq_0: 5e-06 z_p_base_50[0] + [
   0.5 p[0] * z_p_base_50[0] + 0.5 p[1] * z_p_base_50[0]
   + 0.5 p[2] * z_p_base_50[0] - 0.5 p[3] * z_p_base_50[0]
   - 0.5 p[4] * z_p_base_50[0] ] <= 0
 infeasibility_p_base_50_row_1_leq_0: 5e-06 z_p_base_50[1] + [
   0.5 p[0] * z_p_base_50[1] - 0.5 p[1] * z_p_base_50[1]
   - 0.5 p[2] * z_p_base_50[1] + 0.5 p[3] * z_p_base_50[1]
   + 0.5 p[4] * z_p_base_50[1] ] <= 0
 infeasibility_p_base_50_row_2_leq_0: 5e-06 z_p_base_50[2] + [
   - 0.5 p[0] * z_p_base_50[2] + 0.5 p[1] * z_p_base_50[2]
   - 0.5 p[2] * z_p_base_50[2] + 0.5 p[3] * z_p_base_50[2]
   - 0.5 p[4] * z_p_base_50[2] ] <= 0
 infeasibility_p_base_50_row_3_leq_0: 5e-06 z_p_base_50[3] + [
   p[4] * z_p_base_50[3] ] <= 0
 infeasibility_p_base_50_row_4_leq_0: 5e-06 z_p_base_50[4] + [
   p[2] * z_p_base_50[4] ] <= 0
 infeasibility_p_prime_base_50_row_0_leq_0: 5e-06 z_p_prime_base_50[0] + [
   0.5 alpha[0] * z_p_prime_base_50[0]
   - 0.5 alpha[1] * z_p_prime_base_50[0]
   - 0.5 alpha[2] * z_p_prime_base_50[0]
   + 0.5 alpha[3] * z_p_prime_base_50[0]
   - 0.5 alpha[4] * z_p_prime_base_50[0] + 0.5 p[0] * z_p_prime_base_50[0]
   + 0.5 p[1] * z_p_prime_base_50[0] + 0.5 p[2] * z_p_prime_base_50[0]
   - 0.5 p[3] * z_p_prime_base_50[0] - 0.5 p[4] * z_p_prime_base_50[0] ]
   <= 0
 infeasibility_p_prime_base_50_row_1_leq_0: 5e-06 z_p_prime_base_50[1] + [
   0.5 alpha[0] * z_p_prime_base_50[1]
   + 0.5 alpha[1] * z_p_prime_base_50[1]
   + 0.5 alpha[2] * z_p_prime_base_50[1]
   - 0.5 alpha[3] * z_p_prime_base_50[1]
   + 0.5 alpha[4] * z_p_prime_base_50[1] + 0.5 p[0] * z_p_prime_base_50[1]
   - 0.5 p[1] * z_p_prime_base_50[1] - 0.5 p[2] * z_p_prime_base_50[1]
   + 0.5 p[3] * z_p_prime_base_50[1] + 0.5 p[4] * z_p_prime_base_50[1] ]
   <= 0
 infeasibility_p_prime_base_50_row_2_leq_0: 5e-06 z_p_prime_base_50[2] + [
   - 0.5 alpha[0] * z_p_prime_base_50[2]
   - 0.5 alpha[1] * z_p_prime_base_50[2]
   + 0.5 alpha[2] * z_p_prime_base_50[2]
   - 0.5 alpha[3] * z_p_prime_base_50[2]
   - 0.5 alpha[4] * z_p_prime_base_50[2] - 0.5 p[0] * z_p_prime_base_50[2]
   + 0.5 p[1] * z_p_prime_base_50[2] - 0.5 p[2] * z_p_prime_base_50[2]
   + 0.5 p[3] * z_p_prime_base_50[2] - 0.5 p[4] * z_p_prime_base_50[2] ]
   <= 0
 infeasibility_p_prime_base_50_row_3_leq_0: 5e-06 z_p_prime_base_50[3] + [
   alpha[4] * z_p_prime_base_50[3] + p[4] * z_p_prime_base_50[3] ] <= 0
 infeasibility_p_prime_base_50_row_4_leq_0: 5e-06 z_p_prime_base_50[4] + [
   - alpha[2] * z_p_prime_base_50[4] + p[2] * z_p_prime_base_50[4] ] <= 0
 infeasibility_p_base_51_row_0_leq_0: 5e-06 z_p_base_51[0] + [
   0.5 p[0] * z_p_base_51[0] + 0.5 p[1] * z_p_base_51[0]
   + 0.5 p[2] * z_p_base_51[0] - 0.5 p[3] * z_p_base_51[0]
   - 0.5 p[4] * z_p_base_51[0] ] <= 0
 infeasibility_p_base_51_row_1_leq_0: 5e-06 z_p_base_51[1] + [
   0.5 p[0] * z_p_base_51[1] - 0.5 p[1] * z_p_base_51[1]
   - 0.5 p[2] * z_p_base_51[1] + 0.5 p[3] * z_p_base_51[1]
   + 0.5 p[4] * z_p_base_51[1] ] <= 0
 infeasibility_p_base_51_row_2_leq_0: 5e-06 z_p_base_51[2] + [
   - 0.5 p[0] * z_p_base_51[2] + 0.5 p[1] * z_p_base_51[2]
   + 0.5 p[2] * z_p_base_51[2] + 0.5 p[3] * z_p_base_51[2]
   - 0.5 p[4] * z_p_base_51[2] ] <= 0
 infeasibility_p_base_51_row_3_leq_0: 5e-06 z_p_base_51[3] + [
   - p[2] * z_p_base_51[3] + p[4] * z_p_base_51[3] ] <= 0
 infeasibility_p_base_51_row_4_leq_0: 5e-06 z_p_base_51[4] + [
   p[2] * z_p_base_51[4] ] <= 0
 infeasibility_p_prime_base_51_row_0_leq_0: 5e-06 z_p_prime_base_51[0] + [
   0.5 alpha[0] * z_p_prime_base_51[0]
   - 0.5 alpha[1] * z_p_prime_base_51[0]
   - 0.5 alpha[2] * z_p_prime_base_51[0]
   + 0.5 alpha[3] * z_p_prime_base_51[0]
   - 0.5 alpha[4] * z_p_prime_base_51[0] + 0.5 p[0] * z_p_prime_base_51[0]
   + 0.5 p[1] * z_p_prime_base_51[0] + 0.5 p[2] * z_p_prime_base_51[0]
   - 0.5 p[3] * z_p_prime_base_51[0] - 0.5 p[4] * z_p_prime_base_51[0] ]
   <= 0
 infeasibility_p_prime_base_51_row_1_leq_0: 5e-06 z_p_prime_base_51[1] + [
   0.5 alpha[0] * z_p_prime_base_51[1]
   + 0.5 alpha[1] * z_p_prime_base_51[1]
   + 0.5 alpha[2] * z_p_prime_base_51[1]
   - 0.5 alpha[3] * z_p_prime_base_51[1]
   + 0.5 alpha[4] * z_p_prime_base_51[1] + 0.5 p[0] * z_p_prime_base_51[1]
   - 0.5 p[1] * z_p_prime_base_51[1] - 0.5 p[2] * z_p_prime_base_51[1]
   + 0.5 p[3] * z_p_prime_base_51[1] + 0.5 p[4] * z_p_prime_base_51[1] ]
   <= 0
 infeasibility_p_prime_base_51_row_2_leq_0: 5e-06 z_p_prime_base_51[2] + [
   - 0.5 alpha[0] * z_p_prime_base_51[2]
   - 0.5 alpha[1] * z_p_prime_base_51[2]
   - 0.5 alpha[2] * z_p_prime_base_51[2]
   - 0.5 alpha[3] * z_p_prime_base_51[2]
   - 0.5 alpha[4] * z_p_prime_base_51[2] - 0.5 p[0] * z_p_prime_base_51[2]
   + 0.5 p[1] * z_p_prime_base_51[2] + 0.5 p[2] * z_p_prime_base_51[2]
   + 0.5 p[3] * z_p_prime_base_51[2] - 0.5 p[4] * z_p_prime_base_51[2] ]
   <= 0
 infeasibility_p_prime_base_51_row_3_leq_0: 5e-06 z_p_prime_base_51[3] + [
   alpha[2] * z_p_prime_base_51[3] + alpha[4] * z_p_prime_base_51[3]
   - p[2] * z_p_prime_base_51[3] + p[4] * z_p_prime_base_51[3] ] <= 0
 infeasibility_p_prime_base_51_row_4_leq_0: 5e-06 z_p_prime_base_51[4] + [
   - alpha[2] * z_p_prime_base_51[4] + p[2] * z_p_prime_base_51[4] ] <= 0
 infeasibility_p_base_52_row_0_leq_0: 5e-06 z_p_base_52[0] + [
   0.5 p[0] * z_p_base_52[0] + 0.5 p[1] * z_p_base_52[0]
   + 0.5 p[2] * z_p_base_52[0] - 0.5 p[3] * z_p_base_52[0]
   - 0.5 p[4] * z_p_base_52[0] ] <= 0
 infeasibility_p_base_52_row_1_leq_0: 5e-06 z_p_base_52[1] + [
   0.5 p[0] * z_p_base_52[1] - 0.5 p[1] * z_p_base_52[1]
   - 0.5 p[2] * z_p_base_52[1] + 0.5 p[3] * z_p_base_52[1]
   + 0.5 p[4] * z_p_base_52[1] ] <= 0
 infeasibility_p_base_52_row_2_leq_0: 5e-06 z_p_base_52[2] + [
   - 0.5 p[0] * z_p_base_52[2] + 0.5 p[1] * z_p_base_52[2]
   - 0.5 p[2] * z_p_base_52[2] + 0.5 p[3] * z_p_base_52[2]
   + 0.5 p[4] * z_p_base_52[2] ] <= 0
 infeasibility_p_base_52_row_3_leq_0: 5e-06 z_p_base_52[3] + [
   p[2] * z_p_base_52[3] - p[4] * z_p_base_52[3] ] <= 0
 infeasibility_p_base_52_row_4_leq_0: 5e-06 z_p_base_52[4] + [
   p[4] * z_p_base_52[4] ] <= 0
 infeasibility_p_prime_base_52_row_0_leq_0: 5e-06 z_p_prime_base_52[0] + [
   0.5 alpha[0] * z_p_prime_base_52[0]
   - 0.5 alpha[1] * z_p_prime_base_52[0]
   - 0.5 alpha[2] * z_p_prime_base_52[0]
   + 0.5 alpha[3] * z_p_prime_base_52[0]
   - 0.5 alpha[4] * z_p_prime_base_52[0] + 0.5 p[0] * z_p_prime_base_52[0]
   + 0.5 p[1] * z_p_prime_base_52[0] + 0.5 p[2] * z_p_prime_base_52[0]
   - 0.5 p[3] * z_p_prime_base_52[0] - 0.5 p[4] * z_p_prime_base_52[0] ]
   <= 0
 infeasibility_p_prime_base_52_row_1_leq_0: 5e-06 z_p_prime_base_52[1] + [
   0.5 alpha[0] * z_p_prime_base_52[1]
   + 0.5 alpha[1] * z_p_prime_base_52[1]
   + 0.5 alpha[2] * z_p_prime_base_52[1]
   - 0.5 alpha[3] * z_p_prime_base_52[1]
   + 0.5 alpha[4] * z_p_prime_base_52[1] + 0.5 p[0] * z_p_prime_base_52[1]
   - 0.5 p[1] * z_p_prime_base_52[1] - 0.5 p[2] * z_p_prime_base_52[1]
   + 0.5 p[3] * z_p_prime_base_52[1] + 0.5 p[4] * z_p_prime_base_52[1] ]
   <= 0
 infeasibility_p_prime_base_52_row_2_leq_0: 5e-06 z_p_prime_base_52[2] + [
   - 0.5 alpha[0] * z_p_prime_base_52[2]
   - 0.5 alpha[1] * z_p_prime_base_52[2]
   + 0.5 alpha[2] * z_p_prime_base_52[2]
   - 0.5 alpha[3] * z_p_prime_base_52[2]
   + 0.5 alpha[4] * z_p_prime_base_52[2] - 0.5 p[0] * z_p_prime_base_52[2]
   + 0.5 p[1] * z_p_prime_base_52[2] - 0.5 p[2] * z_p_prime_base_52[2]
   + 0.5 p[3] * z_p_prime_base_52[2] + 0.5 p[4] * z_p_prime_base_52[2] ]
   <= 0
 infeasibility_p_prime_base_52_row_3_leq_0: 5e-06 z_p_prime_base_52[3] + [
   - alpha[2] * z_p_prime_base_52[3] - alpha[4] * z_p_prime_base_52[3]
   + p[2] * z_p_prime_base_52[3] - p[4] * z_p_prime_base_52[3] ] <= 0
 infeasibility_p_prime_base_52_row_4_leq_0: 5e-06 z_p_prime_base_52[4] + [
   alpha[4] * z_p_prime_base_52[4] + p[4] * z_p_prime_base_52[4] ] <= 0
 infeasibility_p_base_53_row_0_leq_0: 5e-06 z_p_base_53[0] + [
   0.5 p[0] * z_p_base_53[0] + 0.5 p[1] * z_p_base_53[0]
   + 0.5 p[2] * z_p_base_53[0] - 0.5 p[3] * z_p_base_53[0]
   + 0.5 p[4] * z_p_base_53[0] ] <= 0
 infeasibility_p_base_53_row_1_leq_0: 5e-06 z_p_base_53[1] + [
   0.5 p[0] * z_p_base_53[1] - 0.5 p[1] * z_p_base_53[1]
   - 0.5 p[2] * z_p_base_53[1] + 0.5 p[3] * z_p_base_53[1]
   - 0.5 p[4] * z_p_base_53[1] ] <= 0
 infeasibility_p_base_53_row_2_leq_0: 5e-06 z_p_base_53[2] + [
   - 0.5 p[0] * z_p_base_53[2] + 0.5 p[1] * z_p_base_53[2]
   - 0.5 p[2] * z_p_base_53[2] + 0.5 p[3] * z_p_base_53[2]
   - 0.5 p[4] * z_p_base_53[2] ] <= 0
 infeasibility_p_base_53_row_3_leq_0: 5e-06 z_p_base_53[3] + [
   p[2] * z_p_base_53[3] ] <= 0
 infeasibility_p_base_53_row_4_leq_0: 5e-06 z_p_base_53[4] + [
   p[4] * z_p_base_53[4] ] <= 0
 infeasibility_p_prime_base_53_row_0_leq_0: 5e-06 z_p_prime_base_53[0] + [
   0.5 alpha[0] * z_p_prime_base_53[0]
   - 0.5 alpha[1] * z_p_prime_base_53[0]
   - 0.5 alpha[2] * z_p_prime_base_53[0]
   + 0.5 alpha[3] * z_p_prime_base_53[0]
   + 0.5 alpha[4] * z_p_prime_base_53[0] + 0.5 p[0] * z_p_prime_base_53[0]
   + 0.5 p[1] * z_p_prime_base_53[0] + 0.5 p[2] * z_p_prime_base_53[0]
   - 0.5 p[3] * z_p_prime_base_53[0] + 0.5 p[4] * z_p_prime_base_53[0] ]
   <= 0
 infeasibility_p_prime_base_53_row_1_leq_0: 5e-06 z_p_prime_base_53[1] + [
   0.5 alpha[0] * z_p_prime_base_53[1]
   + 0.5 alpha[1] * z_p_prime_base_53[1]
   + 0.5 alpha[2] * z_p_prime_base_53[1]
   - 0.5 alpha[3] * z_p_prime_base_53[1]
   - 0.5 alpha[4] * z_p_prime_base_53[1] + 0.5 p[0] * z_p_prime_base_53[1]
   - 0.5 p[1] * z_p_prime_base_53[1] - 0.5 p[2] * z_p_prime_base_53[1]
   + 0.5 p[3] * z_p_prime_base_53[1] - 0.5 p[4] * z_p_prime_base_53[1] ]
   <= 0
 infeasibility_p_prime_base_53_row_2_leq_0: 5e-06 z_p_prime_base_53[2] + [
   - 0.5 alpha[0] * z_p_prime_base_53[2]
   - 0.5 alpha[1] * z_p_prime_base_53[2]
   + 0.5 alpha[2] * z_p_prime_base_53[2]
   - 0.5 alpha[3] * z_p_prime_base_53[2]
   - 0.5 alpha[4] * z_p_prime_base_53[2] - 0.5 p[0] * z_p_prime_base_53[2]
   + 0.5 p[1] * z_p_prime_base_53[2] - 0.5 p[2] * z_p_prime_base_53[2]
   + 0.5 p[3] * z_p_prime_base_53[2] - 0.5 p[4] * z_p_prime_base_53[2] ]
   <= 0
 infeasibility_p_prime_base_53_row_3_leq_0: 5e-06 z_p_prime_base_53[3] + [
   - alpha[2] * z_p_prime_base_53[3] + p[2] * z_p_prime_base_53[3] ] <= 0
 infeasibility_p_prime_base_53_row_4_leq_0: 5e-06 z_p_prime_base_53[4] + [
   alpha[4] * z_p_prime_base_53[4] + p[4] * z_p_prime_base_53[4] ] <= 0
 infeasibility_p_base_54_row_0_leq_0: 5e-06 z_p_base_54[0] + [
   0.5 p[0] * z_p_base_54[0] + 0.5 p[1] * z_p_base_54[0]
   - 0.5 p[2] * z_p_base_54[0] - 0.5 p[3] * z_p_base_54[0]
   + 0.5 p[4] * z_p_base_54[0] ] <= 0
 infeasibility_p_base_54_row_1_leq_0: 5e-06 z_p_base_54[1] + [
   0.5 p[0] * z_p_base_54[1] - 0.5 p[1] * z_p_base_54[1]
   + 0.5 p[2] * z_p_base_54[1] + 0.5 p[3] * z_p_base_54[1]
   - 0.5 p[4] * z_p_base_54[1] ] <= 0
 infeasibility_p_base_54_row_2_leq_0: 5e-06 z_p_base_54[2] + [
   - 0.5 p[0] * z_p_base_54[2] + 0.5 p[1] * z_p_base_54[2]
   + 0.5 p[2] * z_p_base_54[2] + 0.5 p[3] * z_p_base_54[2]
   - 0.5 p[4] * z_p_base_54[2] ] <= 0
 infeasibility_p_base_54_row_3_leq_0: 5e-06 z_p_base_54[3] + [
   p[2] * z_p_base_54[3] ] <= 0
 infeasibility_p_base_54_row_4_leq_0: 5e-06 z_p_base_54[4] + [
   - p[2] * z_p_base_54[4] + p[4] * z_p_base_54[4] ] <= 0
 infeasibility_p_prime_base_54_row_0_leq_0: 5e-06 z_p_prime_base_54[0] + [
   0.5 alpha[0] * z_p_prime_base_54[0]
   - 0.5 alpha[1] * z_p_prime_base_54[0]
   + 0.5 alpha[2] * z_p_prime_base_54[0]
   + 0.5 alpha[3] * z_p_prime_base_54[0]
   + 0.5 alpha[4] * z_p_prime_base_54[0] + 0.5 p[0] * z_p_prime_base_54[0]
   + 0.5 p[1] * z_p_prime_base_54[0] - 0.5 p[2] * z_p_prime_base_54[0]
   - 0.5 p[3] * z_p_prime_base_54[0] + 0.5 p[4] * z_p_prime_base_54[0] ]
   <= 0
 infeasibility_p_prime_base_54_row_1_leq_0: 5e-06 z_p_prime_base_54[1] + [
   0.5 alpha[0] * z_p_prime_base_54[1]
   + 0.5 alpha[1] * z_p_prime_base_54[1]
   - 0.5 alpha[2] * z_p_prime_base_54[1]
   - 0.5 alpha[3] * z_p_prime_base_54[1]
   - 0.5 alpha[4] * z_p_prime_base_54[1] + 0.5 p[0] * z_p_prime_base_54[1]
   - 0.5 p[1] * z_p_prime_base_54[1] + 0.5 p[2] * z_p_prime_base_54[1]
   + 0.5 p[3] * z_p_prime_base_54[1] - 0.5 p[4] * z_p_prime_base_54[1] ]
   <= 0
 infeasibility_p_prime_base_54_row_2_leq_0: 5e-06 z_p_prime_base_54[2] + [
   - 0.5 alpha[0] * z_p_prime_base_54[2]
   - 0.5 alpha[1] * z_p_prime_base_54[2]
   - 0.5 alpha[2] * z_p_prime_base_54[2]
   - 0.5 alpha[3] * z_p_prime_base_54[2]
   - 0.5 alpha[4] * z_p_prime_base_54[2] - 0.5 p[0] * z_p_prime_base_54[2]
   + 0.5 p[1] * z_p_prime_base_54[2] + 0.5 p[2] * z_p_prime_base_54[2]
   + 0.5 p[3] * z_p_prime_base_54[2] - 0.5 p[4] * z_p_prime_base_54[2] ]
   <= 0
 infeasibility_p_prime_base_54_row_3_leq_0: 5e-06 z_p_prime_base_54[3] + [
   - alpha[2] * z_p_prime_base_54[3] + p[2] * z_p_prime_base_54[3] ] <= 0
 infeasibility_p_prime_base_54_row_4_leq_0: 5e-06 z_p_prime_base_54[4] + [
   alpha[2] * z_p_prime_base_54[4] + alpha[4] * z_p_prime_base_54[4]
   - p[2] * z_p_prime_base_54[4] + p[4] * z_p_prime_base_54[4] ] <= 0
 infeasibility_p_base_55_row_0_leq_0: 5e-06 z_p_base_55[0] + [
   0.5 p[0] * z_p_base_55[0] + 0.5 p[1] * z_p_base_55[0]
   + 0.5 p[2] * z_p_base_55[0] - 0.5 p[3] * z_p_base_55[0]
   - 0.5 p[4] * z_p_base_55[0] ] <= 0
 infeasibility_p_base_55_row_1_leq_0: 5e-06 z_p_base_55[1] + [
   0.5 p[0] * z_p_base_55[1] - 0.5 p[1] * z_p_base_55[1]
   - 0.5 p[2] * z_p_base_55[1] + 0.5 p[3] * z_p_base_55[1]
   + 0.5 p[4] * z_p_base_55[1] ] <= 0
 infeasibility_p_base_55_row_2_leq_0: 5e-06 z_p_base_55[2] + [
   - 0.5 p[0] * z_p_base_55[2] + 0.5 p[1] * z_p_base_55[2]
   - 0.5 p[2] * z_p_base_55[2] + 0.5 p[3] * z_p_base_55[2]
   + 0.5 p[4] * z_p_base_55[2] ] <= 0
 infeasibility_p_base_55_row_3_leq_0: 5e-06 z_p_base_55[3] + [
   - 0.5 p[0] * z_p_base_55[3] + 0.5 p[1] * z_p_base_55[3]
   + 0.5 p[2] * z_p_base_55[3] + 0.5 p[3] * z_p_base_55[3]
   - 0.5 p[4] * z_p_base_55[3] ] <= 0
 infeasibility_p_base_55_row_4_leq_0: 5e-06 z_p_base_55[4] + [
   0.5 p[0] * z_p_base_55[4] - 0.5 p[1] * z_p_base_55[4]
   + 0.5 p[2] * z_p_base_55[4] - 0.5 p[3] * z_p_base_55[4]
   + 0.5 p[4] * z_p_base_55[4] ] <= 0
 infeasibility_p_prime_base_55_row_0_leq_0: 5e-06 z_p_prime_base_55[0] + [
   0.5 alpha[0] * z_p_prime_base_55[0]
   - 0.5 alpha[1] * z_p_prime_base_55[0]
   - 0.5 alpha[2] * z_p_prime_base_55[0]
   + 0.5 alpha[3] * z_p_prime_base_55[0]
   - 0.5 alpha[4] * z_p_prime_base_55[0] + 0.5 p[0] * z_p_prime_base_55[0]
   + 0.5 p[1] * z_p_prime_base_55[0] + 0.5 p[2] * z_p_prime_base_55[0]
   - 0.5 p[3] * z_p_prime_base_55[0] - 0.5 p[4] * z_p_prime_base_55[0] ]
   <= 0
 infeasibility_p_prime_base_55_row_1_leq_0: 5e-06 z_p_prime_base_55[1] + [
   0.5 alpha[0] * z_p_prime_base_55[1]
   + 0.5 alpha[1] * z_p_prime_base_55[1]
   + 0.5 alpha[2] * z_p_prime_base_55[1]
   - 0.5 alpha[3] * z_p_prime_base_55[1]
   + 0.5 alpha[4] * z_p_prime_base_55[1] + 0.5 p[0] * z_p_prime_base_55[1]
   - 0.5 p[1] * z_p_prime_base_55[1] - 0.5 p[2] * z_p_prime_base_55[1]
   + 0.5 p[3] * z_p_prime_base_55[1] + 0.5 p[4] * z_p_prime_base_55[1] ]
   <= 0
 infeasibility_p_prime_base_55_row_2_leq_0: 5e-06 z_p_prime_base_55[2] + [
   - 0.5 alpha[0] * z_p_prime_base_55[2]
   - 0.5 alpha[1] * z_p_prime_base_55[2]
   + 0.5 alpha[2] * z_p_prime_base_55[2]
   - 0.5 alpha[3] * z_p_prime_base_55[2]
   + 0.5 alpha[4] * z_p_prime_base_55[2] - 0.5 p[0] * z_p_prime_base_55[2]
   + 0.5 p[1] * z_p_prime_base_55[2] - 0.5 p[2] * z_p_prime_base_55[2]
   + 0.5 p[3] * z_p_prime_base_55[2] + 0.5 p[4] * z_p_prime_base_55[2] ]
   <= 0
 infeasibility_p_prime_base_55_row_3_leq_0: 5e-06 z_p_prime_base_55[3] + [
   - 0.5 alpha[0] * z_p_prime_base_55[3]
   - 0.5 alpha[1] * z_p_prime_base_55[3]
   - 0.5 alpha[2] * z_p_prime_base_55[3]
   - 0.5 alpha[3] * z_p_prime_base_55[3]
   - 0.5 alpha[4] * z_p_prime_base_55[3] - 0.5 p[0] * z_p_prime_base_55[3]
   + 0.5 p[1] * z_p_prime_base_55[3] + 0.5 p[2] * z_p_prime_base_55[3]
   + 0.5 p[3] * z_p_prime_base_55[3] - 0.5 p[4] * z_p_prime_base_55[3] ]
   <= 0
 infeasibility_p_prime_base_55_row_4_leq_0: 5e-06 z_p_prime_base_55[4] + [
   0.5 alpha[0] * z_p_prime_base_55[4]
   + 0.5 alpha[1] * z_p_prime_base_55[4]
   - 0.5 alpha[2] * z_p_prime_base_55[4]
   + 0.5 alpha[3] * z_p_prime_base_55[4]
   + 0.5 alpha[4] * z_p_prime_base_55[4] + 0.5 p[0] * z_p_prime_base_55[4]
   - 0.5 p[1] * z_p_prime_base_55[4] + 0.5 p[2] * z_p_prime_base_55[4]
   - 0.5 p[3] * z_p_prime_base_55[4] + 0.5 p[4] * z_p_prime_base_55[4] ]
   <= 0
 infeasibility_p_base_56_row_0_leq_0: 5e-06 z_p_base_56[0] + [
   p[1] * z_p_base_56[0] ] <= 0
 infeasibility_p_base_56_row_1_leq_0: 5e-06 z_p_base_56[1] + [
   p[0] * z_p_base_56[1] - p[1] * z_p_base_56[1] ] <= 0
 infeasibility_p_base_56_row_2_leq_0: 5e-06 z_p_base_56[2] + [
   - 0.5 p[0] * z_p_base_56[2] + 0.5 p[1] * z_p_base_56[2]
   + 0.5 p[2] * z_p_base_56[2] + 0.5 p[3] * z_p_base_56[2]
   - 0.5 p[4] * z_p_base_56[2] ] <= 0
 infeasibility_p_base_56_row_3_leq_0: 5e-06 z_p_base_56[3] + [
   0.5 p[0] * z_p_base_56[3] - 0.5 p[1] * z_p_base_56[3]
   + 0.5 p[2] * z_p_base_56[3] - 0.5 p[3] * z_p_base_56[3]
   + 0.5 p[4] * z_p_base_56[3] ] <= 0
 infeasibility_p_base_56_row_4_leq_0: 5e-06 z_p_base_56[4] + [
   - 0.5 p[0] * z_p_base_56[4] + 0.5 p[1] * z_p_base_56[4]
   - 0.5 p[2] * z_p_base_56[4] + 0.5 p[3] * z_p_base_56[4]
   + 0.5 p[4] * z_p_base_56[4] ] <= 0
 infeasibility_p_prime_base_56_row_0_leq_0: 5e-06 z_p_prime_base_56[0] + [
   - alpha[1] * z_p_prime_base_56[0] + p[1] * z_p_prime_base_56[0] ] <= 0
 infeasibility_p_prime_base_56_row_1_leq_0: 5e-06 z_p_prime_base_56[1] + [
   alpha[0] * z_p_prime_base_56[1] + alpha[1] * z_p_prime_base_56[1]
   + p[0] * z_p_prime_base_56[1] - p[1] * z_p_prime_base_56[1] ] <= 0
 infeasibility_p_prime_base_56_row_2_leq_0: 5e-06 z_p_prime_base_56[2] + [
   - 0.5 alpha[0] * z_p_prime_base_56[2]
   - 0.5 alpha[1] * z_p_prime_base_56[2]
   - 0.5 alpha[2] * z_p_prime_base_56[2]
   - 0.5 alpha[3] * z_p_prime_base_56[2]
   - 0.5 alpha[4] * z_p_prime_base_56[2] - 0.5 p[0] * z_p_prime_base_56[2]
   + 0.5 p[1] * z_p_prime_base_56[2] + 0.5 p[2] * z_p_prime_base_56[2]
   + 0.5 p[3] * z_p_prime_base_56[2] - 0.5 p[4] * z_p_prime_base_56[2] ]
   <= 0
 infeasibility_p_prime_base_56_row_3_leq_0: 5e-06 z_p_prime_base_56[3] + [
   0.5 alpha[0] * z_p_prime_base_56[3]
   + 0.5 alpha[1] * z_p_prime_base_56[3]
   - 0.5 alpha[2] * z_p_prime_base_56[3]
   + 0.5 alpha[3] * z_p_prime_base_56[3]
   + 0.5 alpha[4] * z_p_prime_base_56[3] + 0.5 p[0] * z_p_prime_base_56[3]
   - 0.5 p[1] * z_p_prime_base_56[3] + 0.5 p[2] * z_p_prime_base_56[3]
   - 0.5 p[3] * z_p_prime_base_56[3] + 0.5 p[4] * z_p_prime_base_56[3] ]
   <= 0
 infeasibility_p_prime_base_56_row_4_leq_0: 5e-06 z_p_prime_base_56[4] + [
   - 0.5 alpha[0] * z_p_prime_base_56[4]
   - 0.5 alpha[1] * z_p_prime_base_56[4]
   + 0.5 alpha[2] * z_p_prime_base_56[4]
   - 0.5 alpha[3] * z_p_prime_base_56[4]
   + 0.5 alpha[4] * z_p_prime_base_56[4] - 0.5 p[0] * z_p_prime_base_56[4]
   + 0.5 p[1] * z_p_prime_base_56[4] - 0.5 p[2] * z_p_prime_base_56[4]
   + 0.5 p[3] * z_p_prime_base_56[4] + 0.5 p[4] * z_p_prime_base_56[4] ]
   <= 0
 infeasibility_p_base_57_row_0_leq_0: 5e-06 z_p_base_57[0] + [
   0.5 p[0] * z_p_base_57[0] + 0.5 p[1] * z_p_base_57[0]
   - 0.5 p[2] * z_p_base_57[0] - 0.5 p[3] * z_p_base_57[0]
   - 0.5 p[4] * z_p_base_57[0] ] <= 0
 infeasibility_p_base_57_row_1_leq_0: 5e-06 z_p_base_57[1] + [
   0.5 p[0] * z_p_base_57[1] - 0.5 p[1] * z_p_base_57[1]
   + 0.5 p[2] * z_p_base_57[1] + 0.5 p[3] * z_p_base_57[1]
   + 0.5 p[4] * z_p_base_57[1] ] <= 0
 infeasibility_p_base_57_row_2_leq_0: 5e-06 z_p_base_57[2] + [
   p[2] * z_p_base_57[2] ] <= 0
 infeasibility_p_base_57_row_3_leq_0: 5e-06 z_p_base_57[3] + [
   p[3] * z_p_base_57[3] ] <= 0
 infeasibility_p_base_57_row_4_leq_0: 5e-06 z_p_base_57[4] + [
   - 0.5 p[0] * z_p_base_57[4] + 0.5 p[1] * z_p_base_57[4]
   - 0.5 p[2] * z_p_base_57[4] - 0.5 p[3] * z_p_base_57[4]
   + 0.5 p[4] * z_p_base_57[4] ] <= 0
 infeasibility_p_prime_base_57_row_0_leq_0: 5e-06 z_p_prime_base_57[0] + [
   0.5 alpha[0] * z_p_prime_base_57[0]
   - 0.5 alpha[1] * z_p_prime_base_57[0]
   + 0.5 alpha[2] * z_p_prime_base_57[0]
   + 0.5 alpha[3] * z_p_prime_base_57[0]
   - 0.5 alpha[4] * z_p_prime_base_57[0] + 0.5 p[0] * z_p_prime_base_57[0]
   + 0.5 p[1] * z_p_prime_base_57[0] - 0.5 p[2] * z_p_prime_base_57[0]
   - 0.5 p[3] * z_p_prime_base_57[0] - 0.5 p[4] * z_p_prime_base_57[0] ]
   <= 0
 infeasibility_p_prime_base_57_row_1_leq_0: 5e-06 z_p_prime_base_57[1] + [
   0.5 alpha[0] * z_p_prime_base_57[1]
   + 0.5 alpha[1] * z_p_prime_base_57[1]
   - 0.5 alpha[2] * z_p_prime_base_57[1]
   - 0.5 alpha[3] * z_p_prime_base_57[1]
   + 0.5 alpha[4] * z_p_prime_base_57[1] + 0.5 p[0] * z_p_prime_base_57[1]
   - 0.5 p[1] * z_p_prime_base_57[1] + 0.5 p[2] * z_p_prime_base_57[1]
   + 0.5 p[3] * z_p_prime_base_57[1] + 0.5 p[4] * z_p_prime_base_57[1] ]
   <= 0
 infeasibility_p_prime_base_57_row_2_leq_0: 5e-06 z_p_prime_base_57[2] + [
   - alpha[2] * z_p_prime_base_57[2] + p[2] * z_p_prime_base_57[2] ] <= 0
 infeasibility_p_prime_base_57_row_3_leq_0: 5e-06 z_p_prime_base_57[3] + [
   - alpha[3] * z_p_prime_base_57[3] + p[3] * z_p_prime_base_57[3] ] <= 0
 infeasibility_p_prime_base_57_row_4_leq_0: 5e-06 z_p_prime_base_57[4] + [
   - 0.5 alpha[0] * z_p_prime_base_57[4]
   - 0.5 alpha[1] * z_p_prime_base_57[4]
   + 0.5 alpha[2] * z_p_prime_base_57[4]
   + 0.5 alpha[3] * z_p_prime_base_57[4]
   + 0.5 alpha[4] * z_p_prime_base_57[4] - 0.5 p[0] * z_p_prime_base_57[4]
   + 0.5 p[1] * z_p_prime_base_57[4] - 0.5 p[2] * z_p_prime_base_57[4]
   - 0.5 p[3] * z_p_prime_base_57[4] + 0.5 p[4] * z_p_prime_base_57[4] ]
   <= 0
 infeasibility_p_base_58_row_0_leq_0: 5e-06 z_p_base_58[0] + [
   p[0] * z_p_base_58[0] - p[4] * z_p_base_58[0] ] <= 0
 infeasibility_p_base_58_row_1_leq_0: 5e-06 z_p_base_58[1] + [
   p[4] * z_p_base_58[1] ] <= 0
 infeasibility_p_base_58_row_2_leq_0: 5e-06 z_p_base_58[2] + [
   - 0.5 p[0] * z_p_base_58[2] + 0.5 p[1] * z_p_base_58[2]
   + 0.5 p[2] * z_p_base_58[2] - 0.5 p[3] * z_p_base_58[2]
   + 0.5 p[4] * z_p_base_58[2] ] <= 0
 infeasibility_p_base_58_row_3_leq_0: 5e-06 z_p_base_58[3] + [
   - 0.5 p[0] * z_p_base_58[3] + 0.5 p[1] * z_p_base_58[3]
   - 0.5 p[2] * z_p_base_58[3] + 0.5 p[3] * z_p_base_58[3]
   + 0.5 p[4] * z_p_base_58[3] ] <= 0
 infeasibility_p_base_58_row_4_leq_0: 5e-06 z_p_base_58[4] + [
   0.5 p[0] * z_p_base_58[4] - 0.5 p[1] * z_p_base_58[4]
   + 0.5 p[2] * z_p_base_58[4] + 0.5 p[3] * z_p_base_58[4]
   - 0.5 p[4] * z_p_base_58[4] ] <= 0
 infeasibility_p_prime_base_58_row_0_leq_0: 5e-06 z_p_prime_base_58[0] + [
   alpha[0] * z_p_prime_base_58[0] - alpha[4] * z_p_prime_base_58[0]
   + p[0] * z_p_prime_base_58[0] - p[4] * z_p_prime_base_58[0] ] <= 0
 infeasibility_p_prime_base_58_row_1_leq_0: 5e-06 z_p_prime_base_58[1] + [
   alpha[4] * z_p_prime_base_58[1] + p[4] * z_p_prime_base_58[1] ] <= 0
 infeasibility_p_prime_base_58_row_2_leq_0: 5e-06 z_p_prime_base_58[2] + [
   - 0.5 alpha[0] * z_p_prime_base_58[2]
   - 0.5 alpha[1] * z_p_prime_base_58[2]
   - 0.5 alpha[2] * z_p_prime_base_58[2]
   + 0.5 alpha[3] * z_p_prime_base_58[2]
   + 0.5 alpha[4] * z_p_prime_base_58[2] - 0.5 p[0] * z_p_prime_base_58[2]
   + 0.5 p[1] * z_p_prime_base_58[2] + 0.5 p[2] * z_p_prime_base_58[2]
   - 0.5 p[3] * z_p_prime_base_58[2] + 0.5 p[4] * z_p_prime_base_58[2] ]
   <= 0
 infeasibility_p_prime_base_58_row_3_leq_0: 5e-06 z_p_prime_base_58[3] + [
   - 0.5 alpha[0] * z_p_prime_base_58[3]
   - 0.5 alpha[1] * z_p_prime_base_58[3]
   + 0.5 alpha[2] * z_p_prime_base_58[3]
   - 0.5 alpha[3] * z_p_prime_base_58[3]
   + 0.5 alpha[4] * z_p_prime_base_58[3] - 0.5 p[0] * z_p_prime_base_58[3]
   + 0.5 p[1] * z_p_prime_base_58[3] - 0.5 p[2] * z_p_prime_base_58[3]
   + 0.5 p[3] * z_p_prime_base_58[3] + 0.5 p[4] * z_p_prime_base_58[3] ]
   <= 0
 infeasibility_p_prime_base_58_row_4_leq_0: 5e-06 z_p_prime_base_58[4] + [
   0.5 alpha[0] * z_p_prime_base_58[4]
   + 0.5 alpha[1] * z_p_prime_base_58[4]
   - 0.5 alpha[2] * z_p_prime_base_58[4]
   - 0.5 alpha[3] * z_p_prime_base_58[4]
   - 0.5 alpha[4] * z_p_prime_base_58[4] + 0.5 p[0] * z_p_prime_base_58[4]
   - 0.5 p[1] * z_p_prime_base_58[4] + 0.5 p[2] * z_p_prime_base_58[4]
   + 0.5 p[3] * z_p_prime_base_58[4] - 0.5 p[4] * z_p_prime_base_58[4] ]
   <= 0
 infeasibility_p_base_59_row_0_leq_0: 5e-06 z_p_base_59[0] + [
   0.5 p[0] * z_p_base_59[0] + 0.5 p[1] * z_p_base_59[0]
   - 0.5 p[2] * z_p_base_59[0] + 0.5 p[3] * z_p_base_59[0]
   - 0.5 p[4] * z_p_base_59[0] ] <= 0
 infeasibility_p_base_59_row_1_leq_0: 5e-06 z_p_base_59[1] + [
   0.5 p[0] * z_p_base_59[1] - 0.5 p[1] * z_p_base_59[1]
   + 0.5 p[2] * z_p_base_59[1] - 0.5 p[3] * z_p_base_59[1]
   + 0.5 p[4] * z_p_base_59[1] ] <= 0
 infeasibility_p_base_59_row_2_leq_0: 5e-06 z_p_base_59[2] + [
   p[2] * z_p_base_59[2] - p[3] * z_p_base_59[2] ] <= 0
 infeasibility_p_base_59_row_3_leq_0: 5e-06 z_p_base_59[3] + [
   - 0.5 p[0] * z_p_base_59[3] + 0.5 p[1] * z_p_base_59[3]
   - 0.5 p[2] * z_p_base_59[3] + 0.5 p[3] * z_p_base_59[3]
   + 0.5 p[4] * z_p_base_59[3] ] <= 0
 infeasibility_p_base_59_row_4_leq_0: 5e-06 z_p_base_59[4] + [
   p[3] * z_p_base_59[4] ] <= 0
 infeasibility_p_prime_base_59_row_0_leq_0: 5e-06 z_p_prime_base_59[0] + [
   0.5 alpha[0] * z_p_prime_base_59[0]
   - 0.5 alpha[1] * z_p_prime_base_59[0]
   + 0.5 alpha[2] * z_p_prime_base_59[0]
   - 0.5 alpha[3] * z_p_prime_base_59[0]
   - 0.5 alpha[4] * z_p_prime_base_59[0] + 0.5 p[0] * z_p_prime_base_59[0]
   + 0.5 p[1] * z_p_prime_base_59[0] - 0.5 p[2] * z_p_prime_base_59[0]
   + 0.5 p[3] * z_p_prime_base_59[0] - 0.5 p[4] * z_p_prime_base_59[0] ]
   <= 0
 infeasibility_p_prime_base_59_row_1_leq_0: 5e-06 z_p_prime_base_59[1] + [
   0.5 alpha[0] * z_p_prime_base_59[1]
   + 0.5 alpha[1] * z_p_prime_base_59[1]
   - 0.5 alpha[2] * z_p_prime_base_59[1]
   + 0.5 alpha[3] * z_p_prime_base_59[1]
   + 0.5 alpha[4] * z_p_prime_base_59[1] + 0.5 p[0] * z_p_prime_base_59[1]
   - 0.5 p[1] * z_p_prime_base_59[1] + 0.5 p[2] * z_p_prime_base_59[1]
   - 0.5 p[3] * z_p_prime_base_59[1] + 0.5 p[4] * z_p_prime_base_59[1] ]
   <= 0
 infeasibility_p_prime_base_59_row_2_leq_0: 5e-06 z_p_prime_base_59[2] + [
   - alpha[2] * z_p_prime_base_59[2] + alpha[3] * z_p_prime_base_59[2]
   + p[2] * z_p_prime_base_59[2] - p[3] * z_p_prime_base_59[2] ] <= 0
 infeasibility_p_prime_base_59_row_3_leq_0: 5e-06 z_p_prime_base_59[3] + [
   - 0.5 alpha[0] * z_p_prime_base_59[3]
   - 0.5 alpha[1] * z_p_prime_base_59[3]
   + 0.5 alpha[2] * z_p_prime_base_59[3]
   - 0.5 alpha[3] * z_p_prime_base_59[3]
   + 0.5 alpha[4] * z_p_prime_base_59[3] - 0.5 p[0] * z_p_prime_base_59[3]
   + 0.5 p[1] * z_p_prime_base_59[3] - 0.5 p[2] * z_p_prime_base_59[3]
   + 0.5 p[3] * z_p_prime_base_59[3] + 0.5 p[4] * z_p_prime_base_59[3] ]
   <= 0
 infeasibility_p_prime_base_59_row_4_leq_0: 5e-06 z_p_prime_base_59[4] + [
   - alpha[3] * z_p_prime_base_59[4] + p[3] * z_p_prime_base_59[4] ] <= 0
 infeasibility_p_base_60_row_0_leq_0: 5e-06 z_p_base_60[0] + [
   0.5 p[0] * z_p_base_60[0] + 0.5 p[1] * z_p_base_60[0]
   - 0.5 p[2] * z_p_base_60[0] + 0.5 p[3] * z_p_base_60[0]
   - 0.5 p[4] * z_p_base_60[0] ] <= 0
 infeasibility_p_base_60_row_1_leq_0: 5e-06 z_p_base_60[1] + [
   0.5 p[0] * z_p_base_60[1] - 0.5 p[1] * z_p_base_60[1]
   + 0.5 p[2] * z_p_base_60[1] - 0.5 p[3] * z_p_base_60[1]
   + 0.5 p[4] * z_p_base_60[1] ] <= 0
 infeasibility_p_base_60_row_2_leq_0: 5e-06 z_p_base_60[2] + [
   p[2] * z_p_base_60[2] ] <= 0
 infeasibility_p_base_60_row_3_leq_0: 5e-06 z_p_base_60[3] + [
   - 0.5 p[0] * z_p_base_60[3] + 0.5 p[1] * z_p_base_60[3]
   - 0.5 p[2] * z_p_base_60[3] - 0.5 p[3] * z_p_base_60[3]
   + 0.5 p[4] * z_p_base_60[3] ] <= 0
 infeasibility_p_base_60_row_4_leq_0: 5e-06 z_p_base_60[4] + [
   p[3] * z_p_base_60[4] ] <= 0
 infeasibility_p_prime_base_60_row_0_leq_0: 5e-06 z_p_prime_base_60[0] + [
   0.5 alpha[0] * z_p_prime_base_60[0]
   - 0.5 alpha[1] * z_p_prime_base_60[0]
   + 0.5 alpha[2] * z_p_prime_base_60[0]
   - 0.5 alpha[3] * z_p_prime_base_60[0]
   - 0.5 alpha[4] * z_p_prime_base_60[0] + 0.5 p[0] * z_p_prime_base_60[0]
   + 0.5 p[1] * z_p_prime_base_60[0] - 0.5 p[2] * z_p_prime_base_60[0]
   + 0.5 p[3] * z_p_prime_base_60[0] - 0.5 p[4] * z_p_prime_base_60[0] ]
   <= 0
 infeasibility_p_prime_base_60_row_1_leq_0: 5e-06 z_p_prime_base_60[1] + [
   0.5 alpha[0] * z_p_prime_base_60[1]
   + 0.5 alpha[1] * z_p_prime_base_60[1]
   - 0.5 alpha[2] * z_p_prime_base_60[1]
   + 0.5 alpha[3] * z_p_prime_base_60[1]
   + 0.5 alpha[4] * z_p_prime_base_60[1] + 0.5 p[0] * z_p_prime_base_60[1]
   - 0.5 p[1] * z_p_prime_base_60[1] + 0.5 p[2] * z_p_prime_base_60[1]
   - 0.5 p[3] * z_p_prime_base_60[1] + 0.5 p[4] * z_p_prime_base_60[1] ]
   <= 0
 infeasibility_p_prime_base_60_row_2_leq_0: 5e-06 z_p_prime_base_60[2] + [
   - alpha[2] * z_p_prime_base_60[2] + p[2] * z_p_prime_base_60[2] ] <= 0
 infeasibility_p_prime_base_60_row_3_leq_0: 5e-06 z_p_prime_base_60[3] + [
   - 0.5 alpha[0] * z_p_prime_base_60[3]
   - 0.5 alpha[1] * z_p_prime_base_60[3]
   + 0.5 alpha[2] * z_p_prime_base_60[3]
   + 0.5 alpha[3] * z_p_prime_base_60[3]
   + 0.5 alpha[4] * z_p_prime_base_60[3] - 0.5 p[0] * z_p_prime_base_60[3]
   + 0.5 p[1] * z_p_prime_base_60[3] - 0.5 p[2] * z_p_prime_base_60[3]
   - 0.5 p[3] * z_p_prime_base_60[3] + 0.5 p[4] * z_p_prime_base_60[3] ]
   <= 0
 infeasibility_p_prime_base_60_row_4_leq_0: 5e-06 z_p_prime_base_60[4] + [
   - alpha[3] * z_p_prime_base_60[4] + p[3] * z_p_prime_base_60[4] ] <= 0
 infeasibility_p_base_61_row_0_leq_0: 5e-06 z_p_base_61[0] + [
   0.5 p[0] * z_p_base_61[0] + 0.5 p[1] * z_p_base_61[0]
   - 0.5 p[2] * z_p_base_61[0] + 0.5 p[3] * z_p_base_61[0]
   - 0.5 p[4] * z_p_base_61[0] ] <= 0
 infeasibility_p_base_61_row_1_leq_0: 5e-06 z_p_base_61[1] + [
   0.5 p[0] * z_p_base_61[1] - 0.5 p[1] * z_p_base_61[1]
   + 0.5 p[2] * z_p_base_61[1] - 0.5 p[3] * z_p_base_61[1]
   + 0.5 p[4] * z_p_base_61[1] ] <= 0
 infeasibility_p_base_61_row_2_leq_0: 5e-06 z_p_base_61[2] + [
   - 0.5 p[0] * z_p_base_61[2] + 0.5 p[1] * z_p_base_61[2]
   + 0.5 p[2] * z_p_base_61[2] - 0.5 p[3] * z_p_base_61[2]
   + 0.5 p[4] * z_p_base_61[2] ] <= 0
 infeasibility_p_base_61_row_3_leq_0: 5e-06 z_p_base_61[3] + [
   0.5 p[0] * z_p_base_61[3] - 0.5 p[1] * z_p_base_61[3]
   + 0.5 p[2] * z_p_base_61[3] + 0.5 p[3] * z_p_base_61[3]
   - 0.5 p[4] * z_p_base_61[3] ] <= 0
 infeasibility_p_base_61_row_4_leq_0: 5e-06 z_p_base_61[4] + [
   - 0.5 p[0] * z_p_base_61[4] + 0.5 p[1] * z_p_base_61[4]
   - 0.5 p[2] * z_p_base_61[4] + 0.5 p[3] * z_p_base_61[4]
   + 0.5 p[4] * z_p_base_61[4] ] <= 0
 infeasibility_p_prime_base_61_row_0_leq_0: 5e-06 z_p_prime_base_61[0] + [
   0.5 alpha[0] * z_p_prime_base_61[0]
   - 0.5 alpha[1] * z_p_prime_base_61[0]
   + 0.5 alpha[2] * z_p_prime_base_61[0]
   - 0.5 alpha[3] * z_p_prime_base_61[0]
   - 0.5 alpha[4] * z_p_prime_base_61[0] + 0.5 p[0] * z_p_prime_base_61[0]
   + 0.5 p[1] * z_p_prime_base_61[0] - 0.5 p[2] * z_p_prime_base_61[0]
   + 0.5 p[3] * z_p_prime_base_61[0] - 0.5 p[4] * z_p_prime_base_61[0] ]
   <= 0
 infeasibility_p_prime_base_61_row_1_leq_0: 5e-06 z_p_prime_base_61[1] + [
   0.5 alpha[0] * z_p_prime_base_61[1]
   + 0.5 alpha[1] * z_p_prime_base_61[1]
   - 0.5 alpha[2] * z_p_prime_base_61[1]
   + 0.5 alpha[3] * z_p_prime_base_61[1]
   + 0.5 alpha[4] * z_p_prime_base_61[1] + 0.5 p[0] * z_p_prime_base_61[1]
   - 0.5 p[1] * z_p_prime_base_61[1] + 0.5 p[2] * z_p_prime_base_61[1]
   - 0.5 p[3] * z_p_prime_base_61[1] + 0.5 p[4] * z_p_prime_base_61[1] ]
   <= 0
 infeasibility_p_prime_base_61_row_2_leq_0: 5e-06 z_p_prime_base_61[2] + [
   - 0.5 alpha[0] * z_p_prime_base_61[2]
   - 0.5 alpha[1] * z_p_prime_base_61[2]
   - 0.5 alpha[2] * z_p_prime_base_61[2]
   + 0.5 alpha[3] * z_p_prime_base_61[2]
   + 0.5 alpha[4] * z_p_prime_base_61[2] - 0.5 p[0] * z_p_prime_base_61[2]
   + 0.5 p[1] * z_p_prime_base_61[2] + 0.5 p[2] * z_p_prime_base_61[2]
   - 0.5 p[3] * z_p_prime_base_61[2] + 0.5 p[4] * z_p_prime_base_61[2] ]
   <= 0
 infeasibility_p_prime_base_61_row_3_leq_0: 5e-06 z_p_prime_base_61[3] + [
   0.5 alpha[0] * z_p_prime_base_61[3]
   + 0.5 alpha[1] * z_p_prime_base_61[3]
   - 0.5 alpha[2] * z_p_prime_base_61[3]
   - 0.5 alpha[3] * z_p_prime_base_61[3]
   - 0.5 alpha[4] * z_p_prime_base_61[3] + 0.5 p[0] * z_p_prime_base_61[3]
   - 0.5 p[1] * z_p_prime_base_61[3] + 0.5 p[2] * z_p_prime_base_61[3]
   + 0.5 p[3] * z_p_prime_base_61[3] - 0.5 p[4] * z_p_prime_base_61[3] ]
   <= 0
 infeasibility_p_prime_base_61_row_4_leq_0: 5e-06 z_p_prime_base_61[4] + [
   - 0.5 alpha[0] * z_p_prime_base_61[4]
   - 0.5 alpha[1] * z_p_prime_base_61[4]
   + 0.5 alpha[2] * z_p_prime_base_61[4]
   - 0.5 alpha[3] * z_p_prime_base_61[4]
   + 0.5 alpha[4] * z_p_prime_base_61[4] - 0.5 p[0] * z_p_prime_base_61[4]
   + 0.5 p[1] * z_p_prime_base_61[4] - 0.5 p[2] * z_p_prime_base_61[4]
   + 0.5 p[3] * z_p_prime_base_61[4] + 0.5 p[4] * z_p_prime_base_61[4] ]
   <= 0
 infeasibility_p_base_62_row_0_leq_0: 5e-06 z_p_base_62[0] + [
   0.5 p[0] * z_p_base_62[0] + 0.5 p[1] * z_p_base_62[0]
   + 0.5 p[2] * z_p_base_62[0] - 0.5 p[3] * z_p_base_62[0]
   - 0.5 p[4] * z_p_base_62[0] ] <= 0
 infeasibility_p_base_62_row_1_leq_0: 5e-06 z_p_base_62[1] + [
   0.5 p[0] * z_p_base_62[1] - 0.5 p[1] * z_p_base_62[1]
   - 0.5 p[2] * z_p_base_62[1] + 0.5 p[3] * z_p_base_62[1]
   + 0.5 p[4] * z_p_base_62[1] ] <= 0
 infeasibility_p_base_62_row_2_leq_0: 5e-06 z_p_base_62[2] + [
   - p[2] * z_p_base_62[2] + p[3] * z_p_base_62[2] ] <= 0
 infeasibility_p_base_62_row_3_leq_0: 5e-06 z_p_base_62[3] + [
   - 0.5 p[0] * z_p_base_62[3] + 0.5 p[1] * z_p_base_62[3]
   + 0.5 p[2] * z_p_base_62[3] - 0.5 p[3] * z_p_base_62[3]
   + 0.5 p[4] * z_p_base_62[3] ] <= 0
 infeasibility_p_base_62_row_4_leq_0: 5e-06 z_p_base_62[4] + [
   p[2] * z_p_base_62[4] ] <= 0
 infeasibility_p_prime_base_62_row_0_leq_0: 5e-06 z_p_prime_base_62[0] + [
   0.5 alpha[0] * z_p_prime_base_62[0]
   - 0.5 alpha[1] * z_p_prime_base_62[0]
   - 0.5 alpha[2] * z_p_prime_base_62[0]
   + 0.5 alpha[3] * z_p_prime_base_62[0]
   - 0.5 alpha[4] * z_p_prime_base_62[0] + 0.5 p[0] * z_p_prime_base_62[0]
   + 0.5 p[1] * z_p_prime_base_62[0] + 0.5 p[2] * z_p_prime_base_62[0]
   - 0.5 p[3] * z_p_prime_base_62[0] - 0.5 p[4] * z_p_prime_base_62[0] ]
   <= 0
 infeasibility_p_prime_base_62_row_1_leq_0: 5e-06 z_p_prime_base_62[1] + [
   0.5 alpha[0] * z_p_prime_base_62[1]
   + 0.5 alpha[1] * z_p_prime_base_62[1]
   + 0.5 alpha[2] * z_p_prime_base_62[1]
   - 0.5 alpha[3] * z_p_prime_base_62[1]
   + 0.5 alpha[4] * z_p_prime_base_62[1] + 0.5 p[0] * z_p_prime_base_62[1]
   - 0.5 p[1] * z_p_prime_base_62[1] - 0.5 p[2] * z_p_prime_base_62[1]
   + 0.5 p[3] * z_p_prime_base_62[1] + 0.5 p[4] * z_p_prime_base_62[1] ]
   <= 0
 infeasibility_p_prime_base_62_row_2_leq_0: 5e-06 z_p_prime_base_62[2] + [
   alpha[2] * z_p_prime_base_62[2] - alpha[3] * z_p_prime_base_62[2]
   - p[2] * z_p_prime_base_62[2] + p[3] * z_p_prime_base_62[2] ] <= 0
 infeasibility_p_prime_base_62_row_3_leq_0: 5e-06 z_p_prime_base_62[3] + [
   - 0.5 alpha[0] * z_p_prime_base_62[3]
   - 0.5 alpha[1] * z_p_prime_base_62[3]
   - 0.5 alpha[2] * z_p_prime_base_62[3]
   + 0.5 alpha[3] * z_p_prime_base_62[3]
   + 0.5 alpha[4] * z_p_prime_base_62[3] - 0.5 p[0] * z_p_prime_base_62[3]
   + 0.5 p[1] * z_p_prime_base_62[3] + 0.5 p[2] * z_p_prime_base_62[3]
   - 0.5 p[3] * z_p_prime_base_62[3] + 0.5 p[4] * z_p_prime_base_62[3] ]
   <= 0
 infeasibility_p_prime_base_62_row_4_leq_0: 5e-06 z_p_prime_base_62[4] + [
   - alpha[2] * z_p_prime_base_62[4] + p[2] * z_p_prime_base_62[4] ] <= 0
 infeasibility_p_base_63_row_0_leq_0: 5e-06 z_p_base_63[0] + [
   0.5 p[0] * z_p_base_63[0] + 0.5 p[1] * z_p_base_63[0]
   + 0.5 p[2] * z_p_base_63[0] - 0.5 p[3] * z_p_base_63[0]
   - 0.5 p[4] * z_p_base_63[0] ] <= 0
 infeasibility_p_base_63_row_1_leq_0: 5e-06 z_p_base_63[1] + [
   0.5 p[0] * z_p_base_63[1] - 0.5 p[1] * z_p_base_63[1]
   - 0.5 p[2] * z_p_base_63[1] + 0.5 p[3] * z_p_base_63[1]
   + 0.5 p[4] * z_p_base_63[1] ] <= 0
 infeasibility_p_base_63_row_2_leq_0: 5e-06 z_p_base_63[2] + [
   p[3] * z_p_base_63[2] ] <= 0
 infeasibility_p_base_63_row_3_leq_0: 5e-06 z_p_base_63[3] + [
   - 0.5 p[0] * z_p_base_63[3] + 0.5 p[1] * z_p_base_63[3]
   - 0.5 p[2] * z_p_base_63[3] - 0.5 p[3] * z_p_base_63[3]
   + 0.5 p[4] * z_p_base_63[3] ] <= 0
 infeasibility_p_base_63_row_4_leq_0: 5e-06 z_p_base_63[4] + [
   p[2] * z_p_base_63[4] ] <= 0
 infeasibility_p_prime_base_63_row_0_leq_0: 5e-06 z_p_prime_base_63[0] + [
   0.5 alpha[0] * z_p_prime_base_63[0]
   - 0.5 alpha[1] * z_p_prime_base_63[0]
   - 0.5 alpha[2] * z_p_prime_base_63[0]
   + 0.5 alpha[3] * z_p_prime_base_63[0]
   - 0.5 alpha[4] * z_p_prime_base_63[0] + 0.5 p[0] * z_p_prime_base_63[0]
   + 0.5 p[1] * z_p_prime_base_63[0] + 0.5 p[2] * z_p_prime_base_63[0]
   - 0.5 p[3] * z_p_prime_base_63[0] - 0.5 p[4] * z_p_prime_base_63[0] ]
   <= 0
 infeasibility_p_prime_base_63_row_1_leq_0: 5e-06 z_p_prime_base_63[1] + [
   0.5 alpha[0] * z_p_prime_base_63[1]
   + 0.5 alpha[1] * z_p_prime_base_63[1]
   + 0.5 alpha[2] * z_p_prime_base_63[1]
   - 0.5 alpha[3] * z_p_prime_base_63[1]
   + 0.5 alpha[4] * z_p_prime_base_63[1] + 0.5 p[0] * z_p_prime_base_63[1]
   - 0.5 p[1] * z_p_prime_base_63[1] - 0.5 p[2] * z_p_prime_base_63[1]
   + 0.5 p[3] * z_p_prime_base_63[1] + 0.5 p[4] * z_p_prime_base_63[1] ]
   <= 0
 infeasibility_p_prime_base_63_row_2_leq_0: 5e-06 z_p_prime_base_63[2] + [
   - alpha[3] * z_p_prime_base_63[2] + p[3] * z_p_prime_base_63[2] ] <= 0
 infeasibility_p_prime_base_63_row_3_leq_0: 5e-06 z_p_prime_base_63[3] + [
   - 0.5 alpha[0] * z_p_prime_base_63[3]
   - 0.5 alpha[1] * z_p_prime_base_63[3]
   + 0.5 alpha[2] * z_p_prime_base_63[3]
   + 0.5 alpha[3] * z_p_prime_base_63[3]
   + 0.5 alpha[4] * z_p_prime_base_63[3] - 0.5 p[0] * z_p_prime_base_63[3]
   + 0.5 p[1] * z_p_prime_base_63[3] - 0.5 p[2] * z_p_prime_base_63[3]
   - 0.5 p[3] * z_p_prime_base_63[3] + 0.5 p[4] * z_p_prime_base_63[3] ]
   <= 0
 infeasibility_p_prime_base_63_row_4_leq_0: 5e-06 z_p_prime_base_63[4] + [
   - alpha[2] * z_p_prime_base_63[4] + p[2] * z_p_prime_base_63[4] ] <= 0
 infeasibility_p_base_64_row_0_leq_0: 5e-06 z_p_base_64[0] + [
   0.5 p[0] * z_p_base_64[0] + 0.5 p[1] * z_p_base_64[0]
   + 0.5 p[2] * z_p_base_64[0] - 0.5 p[3] * z_p_base_64[0]
   - 0.5 p[4] * z_p_base_64[0] ] <= 0
 infeasibility_p_base_64_row_1_leq_0: 5e-06 z_p_base_64[1] + [
   0.5 p[0] * z_p_base_64[1] - 0.5 p[1] * z_p_base_64[1]
   - 0.5 p[2] * z_p_base_64[1] + 0.5 p[3] * z_p_base_64[1]
   + 0.5 p[4] * z_p_base_64[1] ] <= 0
 infeasibility_p_base_64_row_2_leq_0: 5e-06 z_p_base_64[2] + [
   - 0.5 p[0] * z_p_base_64[2] + 0.5 p[1] * z_p_base_64[2]
   - 0.5 p[2] * z_p_base_64[2] + 0.5 p[3] * z_p_base_64[2]
   + 0.5 p[4] * z_p_base_64[2] ] <= 0
 infeasibility_p_base_64_row_3_leq_0: 5e-06 z_p_base_64[3] + [
   0.5 p[0] * z_p_base_64[3] - 0.5 p[1] * z_p_base_64[3]
   + 0.5 p[2] * z_p_base_64[3] + 0.5 p[3] * z_p_base_64[3]
   - 0.5 p[4] * z_p_base_64[3] ] <= 0
 infeasibility_p_base_64_row_4_leq_0: 5e-06 z_p_base_64[4] + [
   - 0.5 p[0] * z_p_base_64[4] + 0.5 p[1] * z_p_base_64[4]
   + 0.5 p[2] * z_p_base_64[4] - 0.5 p[3] * z_p_base_64[4]
   + 0.5 p[4] * z_p_base_64[4] ] <= 0
 infeasibility_p_prime_base_64_row_0_leq_0: 5e-06 z_p_prime_base_64[0] + [
   0.5 alpha[0] * z_p_prime_base_64[0]
   - 0.5 alpha[1] * z_p_prime_base_64[0]
   - 0.5 alpha[2] * z_p_prime_base_64[0]
   + 0.5 alpha[3] * z_p_prime_base_64[0]
   - 0.5 alpha[4] * z_p_prime_base_64[0] + 0.5 p[0] * z_p_prime_base_64[0]
   + 0.5 p[1] * z_p_prime_base_64[0] + 0.5 p[2] * z_p_prime_base_64[0]
   - 0.5 p[3] * z_p_prime_base_64[0] - 0.5 p[4] * z_p_prime_base_64[0] ]
   <= 0
 infeasibility_p_prime_base_64_row_1_leq_0: 5e-06 z_p_prime_base_64[1] + [
   0.5 alpha[0] * z_p_prime_base_64[1]
   + 0.5 alpha[1] * z_p_prime_base_64[1]
   + 0.5 alpha[2] * z_p_prime_base_64[1]
   - 0.5 alpha[3] * z_p_prime_base_64[1]
   + 0.5 alpha[4] * z_p_prime_base_64[1] + 0.5 p[0] * z_p_prime_base_64[1]
   - 0.5 p[1] * z_p_prime_base_64[1] - 0.5 p[2] * z_p_prime_base_64[1]
   + 0.5 p[3] * z_p_prime_base_64[1] + 0.5 p[4] * z_p_prime_base_64[1] ]
   <= 0
 infeasibility_p_prime_base_64_row_2_leq_0: 5e-06 z_p_prime_base_64[2] + [
   - 0.5 alpha[0] * z_p_prime_base_64[2]
   - 0.5 alpha[1] * z_p_prime_base_64[2]
   + 0.5 alpha[2] * z_p_prime_base_64[2]
   - 0.5 alpha[3] * z_p_prime_base_64[2]
   + 0.5 alpha[4] * z_p_prime_base_64[2] - 0.5 p[0] * z_p_prime_base_64[2]
   + 0.5 p[1] * z_p_prime_base_64[2] - 0.5 p[2] * z_p_prime_base_64[2]
   + 0.5 p[3] * z_p_prime_base_64[2] + 0.5 p[4] * z_p_prime_base_64[2] ]
   <= 0
 infeasibility_p_prime_base_64_row_3_leq_0: 5e-06 z_p_prime_base_64[3] + [
   0.5 alpha[0] * z_p_prime_base_64[3]
   + 0.5 alpha[1] * z_p_prime_base_64[3]
   - 0.5 alpha[2] * z_p_prime_base_64[3]
   - 0.5 alpha[3] * z_p_prime_base_64[3]
   - 0.5 alpha[4] * z_p_prime_base_64[3] + 0.5 p[0] * z_p_prime_base_64[3]
   - 0.5 p[1] * z_p_prime_base_64[3] + 0.5 p[2] * z_p_prime_base_64[3]
   + 0.5 p[3] * z_p_prime_base_64[3] - 0.5 p[4] * z_p_prime_base_64[3] ]
   <= 0
 infeasibility_p_prime_base_64_row_4_leq_0: 5e-06 z_p_prime_base_64[4] + [
   - 0.5 alpha[0] * z_p_prime_base_64[4]
   - 0.5 alpha[1] * z_p_prime_base_64[4]
   - 0.5 alpha[2] * z_p_prime_base_64[4]
   + 0.5 alpha[3] * z_p_prime_base_64[4]
   + 0.5 alpha[4] * z_p_prime_base_64[4] - 0.5 p[0] * z_p_prime_base_64[4]
   + 0.5 p[1] * z_p_prime_base_64[4] + 0.5 p[2] * z_p_prime_base_64[4]
   - 0.5 p[3] * z_p_prime_base_64[4] + 0.5 p[4] * z_p_prime_base_64[4] ]
   <= 0
 infeasibility_p_base_66_row_0_leq_0: 5e-06 z_p_base_66[0] + [
   0.5 p[0] * z_p_base_66[0] + 0.5 p[1] * z_p_base_66[0]
   - 0.5 p[2] * z_p_base_66[0] + 0.5 p[3] * z_p_base_66[0]
   - 0.5 p[4] * z_p_base_66[0] ] <= 0
 infeasibility_p_base_66_row_1_leq_0: 5e-06 z_p_base_66[1] + [
   0.5 p[0] * z_p_base_66[1] - 0.5 p[1] * z_p_base_66[1]
   + 0.5 p[2] * z_p_base_66[1] - 0.5 p[3] * z_p_base_66[1]
   + 0.5 p[4] * z_p_base_66[1] ] <= 0
 infeasibility_p_base_66_row_2_leq_0: 5e-06 z_p_base_66[2] + [
   - 0.5 p[0] * z_p_base_66[2] + 0.5 p[1] * z_p_base_66[2]
   + 0.5 p[2] * z_p_base_66[2] - 0.5 p[3] * z_p_base_66[2]
   + 0.5 p[4] * z_p_base_66[2] ] <= 0
 infeasibility_p_base_66_row_3_leq_0: 5e-06 z_p_base_66[3] + [
   p[2] * z_p_base_66[3] ] <= 0
 infeasibility_p_base_66_row_4_leq_0: 5e-06 z_p_base_66[4] + [
   - p[2] * z_p_base_66[4] + p[3] * z_p_base_66[4] ] <= 0
 infeasibility_p_prime_base_66_row_0_leq_0: 5e-06 z_p_prime_base_66[0] + [
   0.5 alpha[0] * z_p_prime_base_66[0]
   - 0.5 alpha[1] * z_p_prime_base_66[0]
   + 0.5 alpha[2] * z_p_prime_base_66[0]
   - 0.5 alpha[3] * z_p_prime_base_66[0]
   - 0.5 alpha[4] * z_p_prime_base_66[0] + 0.5 p[0] * z_p_prime_base_66[0]
   + 0.5 p[1] * z_p_prime_base_66[0] - 0.5 p[2] * z_p_prime_base_66[0]
   + 0.5 p[3] * z_p_prime_base_66[0] - 0.5 p[4] * z_p_prime_base_66[0] ]
   <= 0
 infeasibility_p_prime_base_66_row_1_leq_0: 5e-06 z_p_prime_base_66[1] + [
   0.5 alpha[0] * z_p_prime_base_66[1]
   + 0.5 alpha[1] * z_p_prime_base_66[1]
   - 0.5 alpha[2] * z_p_prime_base_66[1]
   + 0.5 alpha[3] * z_p_prime_base_66[1]
   + 0.5 alpha[4] * z_p_prime_base_66[1] + 0.5 p[0] * z_p_prime_base_66[1]
   - 0.5 p[1] * z_p_prime_base_66[1] + 0.5 p[2] * z_p_prime_base_66[1]
   - 0.5 p[3] * z_p_prime_base_66[1] + 0.5 p[4] * z_p_prime_base_66[1] ]
   <= 0
 infeasibility_p_prime_base_66_row_2_leq_0: 5e-06 z_p_prime_base_66[2] + [
   - 0.5 alpha[0] * z_p_prime_base_66[2]
   - 0.5 alpha[1] * z_p_prime_base_66[2]
   - 0.5 alpha[2] * z_p_prime_base_66[2]
   + 0.5 alpha[3] * z_p_prime_base_66[2]
   + 0.5 alpha[4] * z_p_prime_base_66[2] - 0.5 p[0] * z_p_prime_base_66[2]
   + 0.5 p[1] * z_p_prime_base_66[2] + 0.5 p[2] * z_p_prime_base_66[2]
   - 0.5 p[3] * z_p_prime_base_66[2] + 0.5 p[4] * z_p_prime_base_66[2] ]
   <= 0
 infeasibility_p_prime_base_66_row_3_leq_0: 5e-06 z_p_prime_base_66[3] + [
   - alpha[2] * z_p_prime_base_66[3] + p[2] * z_p_prime_base_66[3] ] <= 0
 infeasibility_p_prime_base_66_row_4_leq_0: 5e-06 z_p_prime_base_66[4] + [
   alpha[2] * z_p_prime_base_66[4] - alpha[3] * z_p_prime_base_66[4]
   - p[2] * z_p_prime_base_66[4] + p[3] * z_p_prime_base_66[4] ] <= 0
 infeasibility_p_base_67_row_0_leq_0: 5e-06 z_p_base_67[0] + [
   0.5 p[0] * z_p_base_67[0] + 0.5 p[1] * z_p_base_67[0]
   + 0.5 p[2] * z_p_base_67[0] + 0.5 p[3] * z_p_base_67[0]
   - 0.5 p[4] * z_p_base_67[0] ] <= 0
 infeasibility_p_base_67_row_1_leq_0: 5e-06 z_p_base_67[1] + [
   0.5 p[0] * z_p_base_67[1] - 0.5 p[1] * z_p_base_67[1]
   - 0.5 p[2] * z_p_base_67[1] - 0.5 p[3] * z_p_base_67[1]
   + 0.5 p[4] * z_p_base_67[1] ] <= 0
 infeasibility_p_base_67_row_2_leq_0: 5e-06 z_p_base_67[2] + [
   - 0.5 p[0] * z_p_base_67[2] + 0.5 p[1] * z_p_base_67[2]
   - 0.5 p[2] * z_p_base_67[2] - 0.5 p[3] * z_p_base_67[2]
   + 0.5 p[4] * z_p_base_67[2] ] <= 0
 infeasibility_p_base_67_row_3_leq_0: 5e-06 z_p_base_67[3] + [
   p[2] * z_p_base_67[3] ] <= 0
 infeasibility_p_base_67_row_4_leq_0: 5e-06 z_p_base_67[4] + [
   p[3] * z_p_base_67[4] ] <= 0
 infeasibility_p_prime_base_67_row_0_leq_0: 5e-06 z_p_prime_base_67[0] + [
   0.5 alpha[0] * z_p_prime_base_67[0]
   - 0.5 alpha[1] * z_p_prime_base_67[0]
   - 0.5 alpha[2] * z_p_prime_base_67[0]
   - 0.5 alpha[3] * z_p_prime_base_67[0]
   - 0.5 alpha[4] * z_p_prime_base_67[0] + 0.5 p[0] * z_p_prime_base_67[0]
   + 0.5 p[1] * z_p_prime_base_67[0] + 0.5 p[2] * z_p_prime_base_67[0]
   + 0.5 p[3] * z_p_prime_base_67[0] - 0.5 p[4] * z_p_prime_base_67[0] ]
   <= 0
 infeasibility_p_prime_base_67_row_1_leq_0: 5e-06 z_p_prime_base_67[1] + [
   0.5 alpha[0] * z_p_prime_base_67[1]
   + 0.5 alpha[1] * z_p_prime_base_67[1]
   + 0.5 alpha[2] * z_p_prime_base_67[1]
   + 0.5 alpha[3] * z_p_prime_base_67[1]
   + 0.5 alpha[4] * z_p_prime_base_67[1] + 0.5 p[0] * z_p_prime_base_67[1]
   - 0.5 p[1] * z_p_prime_base_67[1] - 0.5 p[2] * z_p_prime_base_67[1]
   - 0.5 p[3] * z_p_prime_base_67[1] + 0.5 p[4] * z_p_prime_base_67[1] ]
   <= 0
 infeasibility_p_prime_base_67_row_2_leq_0: 5e-06 z_p_prime_base_67[2] + [
   - 0.5 alpha[0] * z_p_prime_base_67[2]
   - 0.5 alpha[1] * z_p_prime_base_67[2]
   + 0.5 alpha[2] * z_p_prime_base_67[2]
   + 0.5 alpha[3] * z_p_prime_base_67[2]
   + 0.5 alpha[4] * z_p_prime_base_67[2] - 0.5 p[0] * z_p_prime_base_67[2]
   + 0.5 p[1] * z_p_prime_base_67[2] - 0.5 p[2] * z_p_prime_base_67[2]
   - 0.5 p[3] * z_p_prime_base_67[2] + 0.5 p[4] * z_p_prime_base_67[2] ]
   <= 0
 infeasibility_p_prime_base_67_row_3_leq_0: 5e-06 z_p_prime_base_67[3] + [
   - alpha[2] * z_p_prime_base_67[3] + p[2] * z_p_prime_base_67[3] ] <= 0
 infeasibility_p_prime_base_67_row_4_leq_0: 5e-06 z_p_prime_base_67[4] + [
   - alpha[3] * z_p_prime_base_67[4] + p[3] * z_p_prime_base_67[4] ] <= 0
 infeasibility_p_base_68_row_0_leq_0: 5e-06 z_p_base_68[0] + [
   p[1] * z_p_base_68[0] ] <= 0
 infeasibility_p_base_68_row_1_leq_0: 5e-06 z_p_base_68[1] + [
   p[0] * z_p_base_68[1] - p[1] * z_p_base_68[1] ] <= 0
 infeasibility_p_base_68_row_2_leq_0: 5e-06 z_p_base_68[2] + [
   0.5 p[0] * z_p_base_68[2] - 0.5 p[1] * z_p_base_68[2]
   + 0.5 p[2] * z_p_base_68[2] + 0.5 p[3] * z_p_base_68[2]
   - 0.5 p[4] * z_p_base_68[2] ] <= 0
 infeasibility_p_base_68_row_3_leq_0: 5e-06 z_p_base_68[3] + [
   - 0.5 p[0] * z_p_base_68[3] + 0.5 p[1] * z_p_base_68[3]
   + 0.5 p[2] * z_p_base_68[3] - 0.5 p[3] * z_p_base_68[3]
   + 0.5 p[4] * z_p_base_68[3] ] <= 0
 infeasibility_p_base_68_row_4_leq_0: 5e-06 z_p_base_68[4] + [
   - 0.5 p[0] * z_p_base_68[4] + 0.5 p[1] * z_p_base_68[4]
   - 0.5 p[2] * z_p_base_68[4] + 0.5 p[3] * z_p_base_68[4]
   + 0.5 p[4] * z_p_base_68[4] ] <= 0
 infeasibility_p_prime_base_68_row_0_leq_0: 5e-06 z_p_prime_base_68[0] + [
   - alpha[1] * z_p_prime_base_68[0] + p[1] * z_p_prime_base_68[0] ] <= 0
 infeasibility_p_prime_base_68_row_1_leq_0: 5e-06 z_p_prime_base_68[1] + [
   alpha[0] * z_p_prime_base_68[1] + alpha[1] * z_p_prime_base_68[1]
   + p[0] * z_p_prime_base_68[1] - p[1] * z_p_prime_base_68[1] ] <= 0
 infeasibility_p_prime_base_68_row_2_leq_0: 5e-06 z_p_prime_base_68[2] + [
   0.5 alpha[0] * z_p_prime_base_68[2]
   + 0.5 alpha[1] * z_p_prime_base_68[2]
   - 0.5 alpha[2] * z_p_prime_base_68[2]
   - 0.5 alpha[3] * z_p_prime_base_68[2]
   - 0.5 alpha[4] * z_p_prime_base_68[2] + 0.5 p[0] * z_p_prime_base_68[2]
   - 0.5 p[1] * z_p_prime_base_68[2] + 0.5 p[2] * z_p_prime_base_68[2]
   + 0.5 p[3] * z_p_prime_base_68[2] - 0.5 p[4] * z_p_prime_base_68[2] ]
   <= 0
 infeasibility_p_prime_base_68_row_3_leq_0: 5e-06 z_p_prime_base_68[3] + [
   - 0.5 alpha[0] * z_p_prime_base_68[3]
   - 0.5 alpha[1] * z_p_prime_base_68[3]
   - 0.5 alpha[2] * z_p_prime_base_68[3]
   + 0.5 alpha[3] * z_p_prime_base_68[3]
   + 0.5 alpha[4] * z_p_prime_base_68[3] - 0.5 p[0] * z_p_prime_base_68[3]
   + 0.5 p[1] * z_p_prime_base_68[3] + 0.5 p[2] * z_p_prime_base_68[3]
   - 0.5 p[3] * z_p_prime_base_68[3] + 0.5 p[4] * z_p_prime_base_68[3] ]
   <= 0
 infeasibility_p_prime_base_68_row_4_leq_0: 5e-06 z_p_prime_base_68[4] + [
   - 0.5 alpha[0] * z_p_prime_base_68[4]
   - 0.5 alpha[1] * z_p_prime_base_68[4]
   + 0.5 alpha[2] * z_p_prime_base_68[4]
   - 0.5 alpha[3] * z_p_prime_base_68[4]
   + 0.5 alpha[4] * z_p_prime_base_68[4] - 0.5 p[0] * z_p_prime_base_68[4]
   + 0.5 p[1] * z_p_prime_base_68[4] - 0.5 p[2] * z_p_prime_base_68[4]
   + 0.5 p[3] * z_p_prime_base_68[4] + 0.5 p[4] * z_p_prime_base_68[4] ]
   <= 0
 infeasibility_p_base_69_row_0_leq_0: 5e-06 z_p_base_69[0] + [
   p[0] * z_p_base_69[0] ] <= 0
 infeasibility_p_base_69_row_1_leq_0: 5e-06 z_p_base_69[1] + [
   - 0.5 p[0] * z_p_base_69[1] + 0.5 p[1] * z_p_base_69[1]
   + 0.5 p[2] * z_p_base_69[1] - 0.5 p[3] * z_p_base_69[1]
   - 0.5 p[4] * z_p_base_69[1] ] <= 0
 infeasibility_p_base_69_row_2_leq_0: 5e-06 z_p_base_69[2] + [
   - 0.5 p[0] * z_p_base_69[2] + 0.5 p[1] * z_p_base_69[2]
   - 0.5 p[2] * z_p_base_69[2] + 0.5 p[3] * z_p_base_69[2]
   - 0.5 p[4] * z_p_base_69[2] ] <= 0
 infeasibility_p_base_69_row_3_leq_0: 5e-06 z_p_base_69[3] + [
   p[4] * z_p_base_69[3] ] <= 0
 infeasibility_p_base_69_row_4_leq_0: 5e-06 z_p_base_69[4] + [
   0.5 p[0] * z_p_base_69[4] - 0.5 p[1] * z_p_base_69[4]
   + 0.5 p[2] * z_p_base_69[4] + 0.5 p[3] * z_p_base_69[4]
   + 0.5 p[4] * z_p_base_69[4] ] <= 0
 infeasibility_p_prime_base_69_row_0_leq_0: 5e-06 z_p_prime_base_69[0] + [
   alpha[0] * z_p_prime_base_69[0] + p[0] * z_p_prime_base_69[0] ] <= 0
 infeasibility_p_prime_base_69_row_1_leq_0: 5e-06 z_p_prime_base_69[1] + [
   - 0.5 alpha[0] * z_p_prime_base_69[1]
   - 0.5 alpha[1] * z_p_prime_base_69[1]
   - 0.5 alpha[2] * z_p_prime_base_69[1]
   + 0.5 alpha[3] * z_p_prime_base_69[1]
   - 0.5 alpha[4] * z_p_prime_base_69[1] - 0.5 p[0] * z_p_prime_base_69[1]
   + 0.5 p[1] * z_p_prime_base_69[1] + 0.5 p[2] * z_p_prime_base_69[1]
   - 0.5 p[3] * z_p_prime_base_69[1] - 0.5 p[4] * z_p_prime_base_69[1] ]
   <= 0
 infeasibility_p_prime_base_69_row_2_leq_0: 5e-06 z_p_prime_base_69[2] + [
   - 0.5 alpha[0] * z_p_prime_base_69[2]
   - 0.5 alpha[1] * z_p_prime_base_69[2]
   + 0.5 alpha[2] * z_p_prime_base_69[2]
   - 0.5 alpha[3] * z_p_prime_base_69[2]
   - 0.5 alpha[4] * z_p_prime_base_69[2] - 0.5 p[0] * z_p_prime_base_69[2]
   + 0.5 p[1] * z_p_prime_base_69[2] - 0.5 p[2] * z_p_prime_base_69[2]
   + 0.5 p[3] * z_p_prime_base_69[2] - 0.5 p[4] * z_p_prime_base_69[2] ]
   <= 0
 infeasibility_p_prime_base_69_row_3_leq_0: 5e-06 z_p_prime_base_69[3] + [
   alpha[4] * z_p_prime_base_69[3] + p[4] * z_p_prime_base_69[3] ] <= 0
 infeasibility_p_prime_base_69_row_4_leq_0: 5e-06 z_p_prime_base_69[4] + [
   0.5 alpha[0] * z_p_prime_base_69[4]
   + 0.5 alpha[1] * z_p_prime_base_69[4]
   - 0.5 alpha[2] * z_p_prime_base_69[4]
   - 0.5 alpha[3] * z_p_prime_base_69[4]
   + 0.5 alpha[4] * z_p_prime_base_69[4] + 0.5 p[0] * z_p_prime_base_69[4]
   - 0.5 p[1] * z_p_prime_base_69[4] + 0.5 p[2] * z_p_prime_base_69[4]
   + 0.5 p[3] * z_p_prime_base_69[4] + 0.5 p[4] * z_p_prime_base_69[4] ]
   <= 0
 infeasibility_p_base_70_row_0_leq_0: 5e-06 z_p_base_70[0] + [
   p[0] * z_p_base_70[0] ] <= 0
 infeasibility_p_base_70_row_1_leq_0: 5e-06 z_p_base_70[1] + [
   - 0.5 p[0] * z_p_base_70[1] + 0.5 p[1] * z_p_base_70[1]
   + 0.5 p[2] * z_p_base_70[1] - 0.5 p[3] * z_p_base_70[1]
   - 0.5 p[4] * z_p_base_70[1] ] <= 0
 infeasibility_p_base_70_row_2_leq_0: 5e-06 z_p_base_70[2] + [
   p[3] * z_p_base_70[2] ] <= 0
 infeasibility_p_base_70_row_3_leq_0: 5e-06 z_p_base_70[3] + [
   - 0.5 p[0] * z_p_base_70[3] + 0.5 p[1] * z_p_base_70[3]
   - 0.5 p[2] * z_p_base_70[3] - 0.5 p[3] * z_p_base_70[3]
   + 0.5 p[4] * z_p_base_70[3] ] <= 0
 infeasibility_p_base_70_row_4_leq_0: 5e-06 z_p_base_70[4] + [
   0.5 p[0] * z_p_base_70[4] - 0.5 p[1] * z_p_base_70[4]
   + 0.5 p[2] * z_p_base_70[4] + 0.5 p[3] * z_p_base_70[4]
   + 0.5 p[4] * z_p_base_70[4] ] <= 0
 infeasibility_p_prime_base_70_row_0_leq_0: 5e-06 z_p_prime_base_70[0] + [
   alpha[0] * z_p_prime_base_70[0] + p[0] * z_p_prime_base_70[0] ] <= 0
 infeasibility_p_prime_base_70_row_1_leq_0: 5e-06 z_p_prime_base_70[1] + [
   - 0.5 alpha[0] * z_p_prime_base_70[1]
   - 0.5 alpha[1] * z_p_prime_base_70[1]
   - 0.5 alpha[2] * z_p_prime_base_70[1]
   + 0.5 alpha[3] * z_p_prime_base_70[1]
   - 0.5 alpha[4] * z_p_prime_base_70[1] - 0.5 p[0] * z_p_prime_base_70[1]
   + 0.5 p[1] * z_p_prime_base_70[1] + 0.5 p[2] * z_p_prime_base_70[1]
   - 0.5 p[3] * z_p_prime_base_70[1] - 0.5 p[4] * z_p_prime_base_70[1] ]
   <= 0
 infeasibility_p_prime_base_70_row_2_leq_0: 5e-06 z_p_prime_base_70[2] + [
   - alpha[3] * z_p_prime_base_70[2] + p[3] * z_p_prime_base_70[2] ] <= 0
 infeasibility_p_prime_base_70_row_3_leq_0: 5e-06 z_p_prime_base_70[3] + [
   - 0.5 alpha[0] * z_p_prime_base_70[3]
   - 0.5 alpha[1] * z_p_prime_base_70[3]
   + 0.5 alpha[2] * z_p_prime_base_70[3]
   + 0.5 alpha[3] * z_p_prime_base_70[3]
   + 0.5 alpha[4] * z_p_prime_base_70[3] - 0.5 p[0] * z_p_prime_base_70[3]
   + 0.5 p[1] * z_p_prime_base_70[3] - 0.5 p[2] * z_p_prime_base_70[3]
   - 0.5 p[3] * z_p_prime_base_70[3] + 0.5 p[4] * z_p_prime_base_70[3] ]
   <= 0
 infeasibility_p_prime_base_70_row_4_leq_0: 5e-06 z_p_prime_base_70[4] + [
   0.5 alpha[0] * z_p_prime_base_70[4]
   + 0.5 alpha[1] * z_p_prime_base_70[4]
   - 0.5 alpha[2] * z_p_prime_base_70[4]
   - 0.5 alpha[3] * z_p_prime_base_70[4]
   + 0.5 alpha[4] * z_p_prime_base_70[4] + 0.5 p[0] * z_p_prime_base_70[4]
   - 0.5 p[1] * z_p_prime_base_70[4] + 0.5 p[2] * z_p_prime_base_70[4]
   + 0.5 p[3] * z_p_prime_base_70[4] + 0.5 p[4] * z_p_prime_base_70[4] ]
   <= 0
 infeasibility_p_base_71_row_0_leq_0: 5e-06 z_p_base_71[0] + [
   p[0] * z_p_base_71[0] ] <= 0
 infeasibility_p_base_71_row_1_leq_0: 5e-06 z_p_base_71[1] + [
   p[2] * z_p_base_71[1] ] <= 0
 infeasibility_p_base_71_row_2_leq_0: 5e-06 z_p_base_71[2] + [
   - 0.5 p[0] * z_p_base_71[2] + 0.5 p[1] * z_p_base_71[2]
   - 0.5 p[2] * z_p_base_71[2] + 0.5 p[3] * z_p_base_71[2]
   - 0.5 p[4] * z_p_base_71[2] ] <= 0
 infeasibility_p_base_71_row_3_leq_0: 5e-06 z_p_base_71[3] + [
   - 0.5 p[0] * z_p_base_71[3] + 0.5 p[1] * z_p_base_71[3]
   - 0.5 p[2] * z_p_base_71[3] - 0.5 p[3] * z_p_base_71[3]
   + 0.5 p[4] * z_p_base_71[3] ] <= 0
 infeasibility_p_base_71_row_4_leq_0: 5e-06 z_p_base_71[4] + [
   0.5 p[0] * z_p_base_71[4] - 0.5 p[1] * z_p_base_71[4]
   + 0.5 p[2] * z_p_base_71[4] + 0.5 p[3] * z_p_base_71[4]
   + 0.5 p[4] * z_p_base_71[4] ] <= 0
 infeasibility_p_prime_base_71_row_0_leq_0: 5e-06 z_p_prime_base_71[0] + [
   alpha[0] * z_p_prime_base_71[0] + p[0] * z_p_prime_base_71[0] ] <= 0
 infeasibility_p_prime_base_71_row_1_leq_0: 5e-06 z_p_prime_base_71[1] + [
   - alpha[2] * z_p_prime_base_71[1] + p[2] * z_p_prime_base_71[1] ] <= 0
 infeasibility_p_prime_base_71_row_2_leq_0: 5e-06 z_p_prime_base_71[2] + [
   - 0.5 alpha[0] * z_p_prime_base_71[2]
   - 0.5 alpha[1] * z_p_prime_base_71[2]
   + 0.5 alpha[2] * z_p_prime_base_71[2]
   - 0.5 alpha[3] * z_p_prime_base_71[2]
   - 0.5 alpha[4] * z_p_prime_base_71[2] - 0.5 p[0] * z_p_prime_base_71[2]
   + 0.5 p[1] * z_p_prime_base_71[2] - 0.5 p[2] * z_p_prime_base_71[2]
   + 0.5 p[3] * z_p_prime_base_71[2] - 0.5 p[4] * z_p_prime_base_71[2] ]
   <= 0
 infeasibility_p_prime_base_71_row_3_leq_0: 5e-06 z_p_prime_base_71[3] + [
   - 0.5 alpha[0] * z_p_prime_base_71[3]
   - 0.5 alpha[1] * z_p_prime_base_71[3]
   + 0.5 alpha[2] * z_p_prime_base_71[3]
   + 0.5 alpha[3] * z_p_prime_base_71[3]
   + 0.5 alpha[4] * z_p_prime_base_71[3] - 0.5 p[0] * z_p_prime_base_71[3]
   + 0.5 p[1] * z_p_prime_base_71[3] - 0.5 p[2] * z_p_prime_base_71[3]
   - 0.5 p[3] * z_p_prime_base_71[3] + 0.5 p[4] * z_p_prime_base_71[3] ]
   <= 0
 infeasibility_p_prime_base_71_row_4_leq_0: 5e-06 z_p_prime_base_71[4] + [
   0.5 alpha[0] * z_p_prime_base_71[4]
   + 0.5 alpha[1] * z_p_prime_base_71[4]
   - 0.5 alpha[2] * z_p_prime_base_71[4]
   - 0.5 alpha[3] * z_p_prime_base_71[4]
   + 0.5 alpha[4] * z_p_prime_base_71[4] + 0.5 p[0] * z_p_prime_base_71[4]
   - 0.5 p[1] * z_p_prime_base_71[4] + 0.5 p[2] * z_p_prime_base_71[4]
   + 0.5 p[3] * z_p_prime_base_71[4] + 0.5 p[4] * z_p_prime_base_71[4] ]
   <= 0
 infeasibility_p_base_72_row_0_leq_0: 5e-06 z_p_base_72[0] + [
   p[0] * z_p_base_72[0] ] <= 0
 infeasibility_p_base_72_row_1_leq_0: 5e-06 z_p_base_72[1] + [
   - 0.5 p[0] * z_p_base_72[1] + 0.5 p[1] * z_p_base_72[1]
   + 0.5 p[2] * z_p_base_72[1] - 0.5 p[3] * z_p_base_72[1]
   - 0.5 p[4] * z_p_base_72[1] ] <= 0
 infeasibility_p_base_72_row_2_leq_0: 5e-06 z_p_base_72[2] + [
   - 0.5 p[0] * z_p_base_72[2] + 0.5 p[1] * z_p_base_72[2]
   - 0.5 p[2] * z_p_base_72[2] + 0.5 p[3] * z_p_base_72[2]
   + 0.5 p[4] * z_p_base_72[2] ] <= 0
 infeasibility_p_base_72_row_3_leq_0: 5e-06 z_p_base_72[3] + [
   0.5 p[0] * z_p_base_72[3] - 0.5 p[1] * z_p_base_72[3]
   + 0.5 p[2] * z_p_base_72[3] + 0.5 p[3] * z_p_base_72[3]
   - 0.5 p[4] * z_p_base_72[3] ] <= 0
 infeasibility_p_base_72_row_4_leq_0: 5e-06 z_p_base_72[4] + [
   p[4] * z_p_base_72[4] ] <= 0
 infeasibility_p_prime_base_72_row_0_leq_0: 5e-06 z_p_prime_base_72[0] + [
   alpha[0] * z_p_prime_base_72[0] + p[0] * z_p_prime_base_72[0] ] <= 0
 infeasibility_p_prime_base_72_row_1_leq_0: 5e-06 z_p_prime_base_72[1] + [
   - 0.5 alpha[0] * z_p_prime_base_72[1]
   - 0.5 alpha[1] * z_p_prime_base_72[1]
   - 0.5 alpha[2] * z_p_prime_base_72[1]
   + 0.5 alpha[3] * z_p_prime_base_72[1]
   - 0.5 alpha[4] * z_p_prime_base_72[1] - 0.5 p[0] * z_p_prime_base_72[1]
   + 0.5 p[1] * z_p_prime_base_72[1] + 0.5 p[2] * z_p_prime_base_72[1]
   - 0.5 p[3] * z_p_prime_base_72[1] - 0.5 p[4] * z_p_prime_base_72[1] ]
   <= 0
 infeasibility_p_prime_base_72_row_2_leq_0: 5e-06 z_p_prime_base_72[2] + [
   - 0.5 alpha[0] * z_p_prime_base_72[2]
   - 0.5 alpha[1] * z_p_prime_base_72[2]
   + 0.5 alpha[2] * z_p_prime_base_72[2]
   - 0.5 alpha[3] * z_p_prime_base_72[2]
   + 0.5 alpha[4] * z_p_prime_base_72[2] - 0.5 p[0] * z_p_prime_base_72[2]
   + 0.5 p[1] * z_p_prime_base_72[2] - 0.5 p[2] * z_p_prime_base_72[2]
   + 0.5 p[3] * z_p_prime_base_72[2] + 0.5 p[4] * z_p_prime_base_72[2] ]
   <= 0
 infeasibility_p_prime_base_72_row_3_leq_0: 5e-06 z_p_prime_base_72[3] + [
   0.5 alpha[0] * z_p_prime_base_72[3]
   + 0.5 alpha[1] * z_p_prime_base_72[3]
   - 0.5 alpha[2] * z_p_prime_base_72[3]
   - 0.5 alpha[3] * z_p_prime_base_72[3]
   - 0.5 alpha[4] * z_p_prime_base_72[3] + 0.5 p[0] * z_p_prime_base_72[3]
   - 0.5 p[1] * z_p_prime_base_72[3] + 0.5 p[2] * z_p_prime_base_72[3]
   + 0.5 p[3] * z_p_prime_base_72[3] - 0.5 p[4] * z_p_prime_base_72[3] ]
   <= 0
 infeasibility_p_prime_base_72_row_4_leq_0: 5e-06 z_p_prime_base_72[4] + [
   alpha[4] * z_p_prime_base_72[4] + p[4] * z_p_prime_base_72[4] ] <= 0
 infeasibility_p_base_73_row_0_leq_0: 5e-06 z_p_base_73[0] + [
   p[0] * z_p_base_73[0] ] <= 0
 infeasibility_p_base_73_row_1_leq_0: 5e-06 z_p_base_73[1] + [
   - 0.5 p[0] * z_p_base_73[1] + 0.5 p[1] * z_p_base_73[1]
   + 0.5 p[2] * z_p_base_73[1] - 0.5 p[3] * z_p_base_73[1]
   + 0.5 p[4] * z_p_base_73[1] ] <= 0
 infeasibility_p_base_73_row_2_leq_0: 5e-06 z_p_base_73[2] + [
   - 0.5 p[0] * z_p_base_73[2] + 0.5 p[1] * z_p_base_73[2]
   - 0.5 p[2] * z_p_base_73[2] + 0.5 p[3] * z_p_base_73[2]
   - 0.5 p[4] * z_p_base_73[2] ] <= 0
 infeasibility_p_base_73_row_3_leq_0: 5e-06 z_p_base_73[3] + [
   0.5 p[0] * z_p_base_73[3] - 0.5 p[1] * z_p_base_73[3]
   + 0.5 p[2] * z_p_base_73[3] + 0.5 p[3] * z_p_base_73[3]
   - 0.5 p[4] * z_p_base_73[3] ] <= 0
 infeasibility_p_base_73_row_4_leq_0: 5e-06 z_p_base_73[4] + [
   p[4] * z_p_base_73[4] ] <= 0
 infeasibility_p_prime_base_73_row_0_leq_0: 5e-06 z_p_prime_base_73[0] + [
   alpha[0] * z_p_prime_base_73[0] + p[0] * z_p_prime_base_73[0] ] <= 0
 infeasibility_p_prime_base_73_row_1_leq_0: 5e-06 z_p_prime_base_73[1] + [
   - 0.5 alpha[0] * z_p_prime_base_73[1]
   - 0.5 alpha[1] * z_p_prime_base_73[1]
   - 0.5 alpha[2] * z_p_prime_base_73[1]
   + 0.5 alpha[3] * z_p_prime_base_73[1]
   + 0.5 alpha[4] * z_p_prime_base_73[1] - 0.5 p[0] * z_p_prime_base_73[1]
   + 0.5 p[1] * z_p_prime_base_73[1] + 0.5 p[2] * z_p_prime_base_73[1]
   - 0.5 p[3] * z_p_prime_base_73[1] + 0.5 p[4] * z_p_prime_base_73[1] ]
   <= 0
 infeasibility_p_prime_base_73_row_2_leq_0: 5e-06 z_p_prime_base_73[2] + [
   - 0.5 alpha[0] * z_p_prime_base_73[2]
   - 0.5 alpha[1] * z_p_prime_base_73[2]
   + 0.5 alpha[2] * z_p_prime_base_73[2]
   - 0.5 alpha[3] * z_p_prime_base_73[2]
   - 0.5 alpha[4] * z_p_prime_base_73[2] - 0.5 p[0] * z_p_prime_base_73[2]
   + 0.5 p[1] * z_p_prime_base_73[2] - 0.5 p[2] * z_p_prime_base_73[2]
   + 0.5 p[3] * z_p_prime_base_73[2] - 0.5 p[4] * z_p_prime_base_73[2] ]
   <= 0
 infeasibility_p_prime_base_73_row_3_leq_0: 5e-06 z_p_prime_base_73[3] + [
   0.5 alpha[0] * z_p_prime_base_73[3]
   + 0.5 alpha[1] * z_p_prime_base_73[3]
   - 0.5 alpha[2] * z_p_prime_base_73[3]
   - 0.5 alpha[3] * z_p_prime_base_73[3]
   - 0.5 alpha[4] * z_p_prime_base_73[3] + 0.5 p[0] * z_p_prime_base_73[3]
   - 0.5 p[1] * z_p_prime_base_73[3] + 0.5 p[2] * z_p_prime_base_73[3]
   + 0.5 p[3] * z_p_prime_base_73[3] - 0.5 p[4] * z_p_prime_base_73[3] ]
   <= 0
 infeasibility_p_prime_base_73_row_4_leq_0: 5e-06 z_p_prime_base_73[4] + [
   alpha[4] * z_p_prime_base_73[4] + p[4] * z_p_prime_base_73[4] ] <= 0
 infeasibility_p_base_74_row_0_leq_0: 5e-06 z_p_base_74[0] + [
   p[0] * z_p_base_74[0] ] <= 0
 infeasibility_p_base_74_row_1_leq_0: 5e-06 z_p_base_74[1] + [
   - 0.5 p[0] * z_p_base_74[1] + 0.5 p[1] * z_p_base_74[1]
   + 0.5 p[2] * z_p_base_74[1] - 0.5 p[3] * z_p_base_74[1]
   - 0.5 p[4] * z_p_base_74[1] ] <= 0
 infeasibility_p_base_74_row_2_leq_0: 5e-06 z_p_base_74[2] + [
   - 0.5 p[0] * z_p_base_74[2] + 0.5 p[1] * z_p_base_74[2]
   - 0.5 p[2] * z_p_base_74[2] + 0.5 p[3] * z_p_base_74[2]
   + 0.5 p[4] * z_p_base_74[2] ] <= 0
 infeasibility_p_base_74_row_3_leq_0: 5e-06 z_p_base_74[3] + [
   p[3] * z_p_base_74[3] ] <= 0
 infeasibility_p_base_74_row_4_leq_0: 5e-06 z_p_base_74[4] + [
   0.5 p[0] * z_p_base_74[4] - 0.5 p[1] * z_p_base_74[4]
   + 0.5 p[2] * z_p_base_74[4] - 0.5 p[3] * z_p_base_74[4]
   + 0.5 p[4] * z_p_base_74[4] ] <= 0
 infeasibility_p_prime_base_74_row_0_leq_0: 5e-06 z_p_prime_base_74[0] + [
   alpha[0] * z_p_prime_base_74[0] + p[0] * z_p_prime_base_74[0] ] <= 0
 infeasibility_p_prime_base_74_row_1_leq_0: 5e-06 z_p_prime_base_74[1] + [
   - 0.5 alpha[0] * z_p_prime_base_74[1]
   - 0.5 alpha[1] * z_p_prime_base_74[1]
   - 0.5 alpha[2] * z_p_prime_base_74[1]
   + 0.5 alpha[3] * z_p_prime_base_74[1]
   - 0.5 alpha[4] * z_p_prime_base_74[1] - 0.5 p[0] * z_p_prime_base_74[1]
   + 0.5 p[1] * z_p_prime_base_74[1] + 0.5 p[2] * z_p_prime_base_74[1]
   - 0.5 p[3] * z_p_prime_base_74[1] - 0.5 p[4] * z_p_prime_base_74[1] ]
   <= 0
 infeasibility_p_prime_base_74_row_2_leq_0: 5e-06 z_p_prime_base_74[2] + [
   - 0.5 alpha[0] * z_p_prime_base_74[2]
   - 0.5 alpha[1] * z_p_prime_base_74[2]
   + 0.5 alpha[2] * z_p_prime_base_74[2]
   - 0.5 alpha[3] * z_p_prime_base_74[2]
   + 0.5 alpha[4] * z_p_prime_base_74[2] - 0.5 p[0] * z_p_prime_base_74[2]
   + 0.5 p[1] * z_p_prime_base_74[2] - 0.5 p[2] * z_p_prime_base_74[2]
   + 0.5 p[3] * z_p_prime_base_74[2] + 0.5 p[4] * z_p_prime_base_74[2] ]
   <= 0
 infeasibility_p_prime_base_74_row_3_leq_0: 5e-06 z_p_prime_base_74[3] + [
   - alpha[3] * z_p_prime_base_74[3] + p[3] * z_p_prime_base_74[3] ] <= 0
 infeasibility_p_prime_base_74_row_4_leq_0: 5e-06 z_p_prime_base_74[4] + [
   0.5 alpha[0] * z_p_prime_base_74[4]
   + 0.5 alpha[1] * z_p_prime_base_74[4]
   - 0.5 alpha[2] * z_p_prime_base_74[4]
   + 0.5 alpha[3] * z_p_prime_base_74[4]
   + 0.5 alpha[4] * z_p_prime_base_74[4] + 0.5 p[0] * z_p_prime_base_74[4]
   - 0.5 p[1] * z_p_prime_base_74[4] + 0.5 p[2] * z_p_prime_base_74[4]
   - 0.5 p[3] * z_p_prime_base_74[4] + 0.5 p[4] * z_p_prime_base_74[4] ]
   <= 0
 infeasibility_p_base_75_row_0_leq_0: 5e-06 z_p_base_75[0] + [
   p[0] * z_p_base_75[0] ] <= 0
 infeasibility_p_base_75_row_1_leq_0: 5e-06 z_p_base_75[1] + [
   - 0.5 p[0] * z_p_base_75[1] + 0.5 p[1] * z_p_base_75[1]
   + 0.5 p[2] * z_p_base_75[1] + 0.5 p[3] * z_p_base_75[1]
   - 0.5 p[4] * z_p_base_75[1] ] <= 0
 infeasibility_p_base_75_row_2_leq_0: 5e-06 z_p_base_75[2] + [
   - 0.5 p[0] * z_p_base_75[2] + 0.5 p[1] * z_p_base_75[2]
   - 0.5 p[2] * z_p_base_75[2] - 0.5 p[3] * z_p_base_75[2]
   + 0.5 p[4] * z_p_base_75[2] ] <= 0
 infeasibility_p_base_75_row_3_leq_0: 5e-06 z_p_base_75[3] + [
   0.5 p[0] * z_p_base_75[3] - 0.5 p[1] * z_p_base_75[3]
   + 0.5 p[2] * z_p_base_75[3] - 0.5 p[3] * z_p_base_75[3]
   + 0.5 p[4] * z_p_base_75[3] ] <= 0
 infeasibility_p_base_75_row_4_leq_0: 5e-06 z_p_base_75[4] + [
   p[3] * z_p_base_75[4] ] <= 0
 infeasibility_p_prime_base_75_row_0_leq_0: 5e-06 z_p_prime_base_75[0] + [
   alpha[0] * z_p_prime_base_75[0] + p[0] * z_p_prime_base_75[0] ] <= 0
 infeasibility_p_prime_base_75_row_1_leq_0: 5e-06 z_p_prime_base_75[1] + [
   - 0.5 alpha[0] * z_p_prime_base_75[1]
   - 0.5 alpha[1] * z_p_prime_base_75[1]
   - 0.5 alpha[2] * z_p_prime_base_75[1]
   - 0.5 alpha[3] * z_p_prime_base_75[1]
   - 0.5 alpha[4] * z_p_prime_base_75[1] - 0.5 p[0] * z_p_prime_base_75[1]
   + 0.5 p[1] * z_p_prime_base_75[1] + 0.5 p[2] * z_p_prime_base_75[1]
   + 0.5 p[3] * z_p_prime_base_75[1] - 0.5 p[4] * z_p_prime_base_75[1] ]
   <= 0
 infeasibility_p_prime_base_75_row_2_leq_0: 5e-06 z_p_prime_base_75[2] + [
   - 0.5 alpha[0] * z_p_prime_base_75[2]
   - 0.5 alpha[1] * z_p_prime_base_75[2]
   + 0.5 alpha[2] * z_p_prime_base_75[2]
   + 0.5 alpha[3] * z_p_prime_base_75[2]
   + 0.5 alpha[4] * z_p_prime_base_75[2] - 0.5 p[0] * z_p_prime_base_75[2]
   + 0.5 p[1] * z_p_prime_base_75[2] - 0.5 p[2] * z_p_prime_base_75[2]
   - 0.5 p[3] * z_p_prime_base_75[2] + 0.5 p[4] * z_p_prime_base_75[2] ]
   <= 0
 infeasibility_p_prime_base_75_row_3_leq_0: 5e-06 z_p_prime_base_75[3] + [
   0.5 alpha[0] * z_p_prime_base_75[3]
   + 0.5 alpha[1] * z_p_prime_base_75[3]
   - 0.5 alpha[2] * z_p_prime_base_75[3]
   + 0.5 alpha[3] * z_p_prime_base_75[3]
   + 0.5 alpha[4] * z_p_prime_base_75[3] + 0.5 p[0] * z_p_prime_base_75[3]
   - 0.5 p[1] * z_p_prime_base_75[3] + 0.5 p[2] * z_p_prime_base_75[3]
   - 0.5 p[3] * z_p_prime_base_75[3] + 0.5 p[4] * z_p_prime_base_75[3] ]
   <= 0
 infeasibility_p_prime_base_75_row_4_leq_0: 5e-06 z_p_prime_base_75[4] + [
   - alpha[3] * z_p_prime_base_75[4] + p[3] * z_p_prime_base_75[4] ] <= 0
 infeasibility_p_base_76_row_0_leq_0: 5e-06 z_p_base_76[0] + [
   p[0] * z_p_base_76[0] ] <= 0
 infeasibility_p_base_76_row_1_leq_0: 5e-06 z_p_base_76[1] + [
   - p[0] * z_p_base_76[1] + p[1] * z_p_base_76[1] ] <= 0
 infeasibility_p_base_76_row_2_leq_0: 5e-06 z_p_base_76[2] + [
   0.5 p[0] * z_p_base_76[2] - 0.5 p[1] * z_p_base_76[2]
   + 0.5 p[2] * z_p_base_76[2] + 0.5 p[3] * z_p_base_76[2]
   - 0.5 p[4] * z_p_base_76[2] ] <= 0
 infeasibility_p_base_76_row_3_leq_0: 5e-06 z_p_base_76[3] + [
   0.5 p[0] * z_p_base_76[3] - 0.5 p[1] * z_p_base_76[3]
   + 0.5 p[2] * z_p_base_76[3] - 0.5 p[3] * z_p_base_76[3]
   + 0.5 p[4] * z_p_base_76[3] ] <= 0
 infeasibility_p_base_76_row_4_leq_0: 5e-06 z_p_base_76[4] + [
   - 0.5 p[0] * z_p_base_76[4] + 0.5 p[1] * z_p_base_76[4]
   - 0.5 p[2] * z_p_base_76[4] + 0.5 p[3] * z_p_base_76[4]
   + 0.5 p[4] * z_p_base_76[4] ] <= 0
 infeasibility_p_prime_base_76_row_0_leq_0: 5e-06 z_p_prime_base_76[0] + [
   alpha[0] * z_p_prime_base_76[0] + p[0] * z_p_prime_base_76[0] ] <= 0
 infeasibility_p_prime_base_76_row_1_leq_0: 5e-06 z_p_prime_base_76[1] + [
   - alpha[0] * z_p_prime_base_76[1] - alpha[1] * z_p_prime_base_76[1]
   - p[0] * z_p_prime_base_76[1] + p[1] * z_p_prime_base_76[1] ] <= 0
 infeasibility_p_prime_base_76_row_2_leq_0: 5e-06 z_p_prime_base_76[2] + [
   0.5 alpha[0] * z_p_prime_base_76[2]
   + 0.5 alpha[1] * z_p_prime_base_76[2]
   - 0.5 alpha[2] * z_p_prime_base_76[2]
   - 0.5 alpha[3] * z_p_prime_base_76[2]
   - 0.5 alpha[4] * z_p_prime_base_76[2] + 0.5 p[0] * z_p_prime_base_76[2]
   - 0.5 p[1] * z_p_prime_base_76[2] + 0.5 p[2] * z_p_prime_base_76[2]
   + 0.5 p[3] * z_p_prime_base_76[2] - 0.5 p[4] * z_p_prime_base_76[2] ]
   <= 0
 infeasibility_p_prime_base_76_row_3_leq_0: 5e-06 z_p_prime_base_76[3] + [
   0.5 alpha[0] * z_p_prime_base_76[3]
   + 0.5 alpha[1] * z_p_prime_base_76[3]
   - 0.5 alpha[2] * z_p_prime_base_76[3]
   + 0.5 alpha[3] * z_p_prime_base_76[3]
   + 0.5 alpha[4] * z_p_prime_base_76[3] + 0.5 p[0] * z_p_prime_base_76[3]
   - 0.5 p[1] * z_p_prime_base_76[3] + 0.5 p[2] * z_p_prime_base_76[3]
   - 0.5 p[3] * z_p_prime_base_76[3] + 0.5 p[4] * z_p_prime_base_76[3] ]
   <= 0
 infeasibility_p_prime_base_76_row_4_leq_0: 5e-06 z_p_prime_base_76[4] + [
   - 0.5 alpha[0] * z_p_prime_base_76[4]
   - 0.5 alpha[1] * z_p_prime_base_76[4]
   + 0.5 alpha[2] * z_p_prime_base_76[4]
   - 0.5 alpha[3] * z_p_prime_base_76[4]
   + 0.5 alpha[4] * z_p_prime_base_76[4] - 0.5 p[0] * z_p_prime_base_76[4]
   + 0.5 p[1] * z_p_prime_base_76[4] - 0.5 p[2] * z_p_prime_base_76[4]
   + 0.5 p[3] * z_p_prime_base_76[4] + 0.5 p[4] * z_p_prime_base_76[4] ]
   <= 0
 infeasibility_p_base_77_row_0_leq_0: 5e-06 z_p_base_77[0] + [
   p[0] * z_p_base_77[0] ] <= 0
 infeasibility_p_base_77_row_1_leq_0: 5e-06 z_p_base_77[1] + [
   - 0.5 p[0] * z_p_base_77[1] + 0.5 p[1] * z_p_base_77[1]
   - 0.5 p[2] * z_p_base_77[1] + 0.5 p[3] * z_p_base_77[1]
   - 0.5 p[4] * z_p_base_77[1] ] <= 0
 infeasibility_p_base_77_row_2_leq_0: 5e-06 z_p_base_77[2] + [
   - 0.5 p[0] * z_p_base_77[2] + 0.5 p[1] * z_p_base_77[2]
   + 0.5 p[2] * z_p_base_77[2] - 0.5 p[3] * z_p_base_77[2]
   + 0.5 p[4] * z_p_base_77[2] ] <= 0
 infeasibility_p_base_77_row_3_leq_0: 5e-06 z_p_base_77[3] + [
   p[2] * z_p_base_77[3] ] <= 0
 infeasibility_p_base_77_row_4_leq_0: 5e-06 z_p_base_77[4] + [
   0.5 p[0] * z_p_base_77[4] - 0.5 p[1] * z_p_base_77[4]
   - 0.5 p[2] * z_p_base_77[4] + 0.5 p[3] * z_p_base_77[4]
   + 0.5 p[4] * z_p_base_77[4] ] <= 0
 infeasibility_p_prime_base_77_row_0_leq_0: 5e-06 z_p_prime_base_77[0] + [
   alpha[0] * z_p_prime_base_77[0] + p[0] * z_p_prime_base_77[0] ] <= 0
 infeasibility_p_prime_base_77_row_1_leq_0: 5e-06 z_p_prime_base_77[1] + [
   - 0.5 alpha[0] * z_p_prime_base_77[1]
   - 0.5 alpha[1] * z_p_prime_base_77[1]
   + 0.5 alpha[2] * z_p_prime_base_77[1]
   - 0.5 alpha[3] * z_p_prime_base_77[1]
   - 0.5 alpha[4] * z_p_prime_base_77[1] - 0.5 p[0] * z_p_prime_base_77[1]
   + 0.5 p[1] * z_p_prime_base_77[1] - 0.5 p[2] * z_p_prime_base_77[1]
   + 0.5 p[3] * z_p_prime_base_77[1] - 0.5 p[4] * z_p_prime_base_77[1] ]
   <= 0
 infeasibility_p_prime_base_77_row_2_leq_0: 5e-06 z_p_prime_base_77[2] + [
   - 0.5 alpha[0] * z_p_prime_base_77[2]
   - 0.5 alpha[1] * z_p_prime_base_77[2]
   - 0.5 alpha[2] * z_p_prime_base_77[2]
   + 0.5 alpha[3] * z_p_prime_base_77[2]
   + 0.5 alpha[4] * z_p_prime_base_77[2] - 0.5 p[0] * z_p_prime_base_77[2]
   + 0.5 p[1] * z_p_prime_base_77[2] + 0.5 p[2] * z_p_prime_base_77[2]
   - 0.5 p[3] * z_p_prime_base_77[2] + 0.5 p[4] * z_p_prime_base_77[2] ]
   <= 0
 infeasibility_p_prime_base_77_row_3_leq_0: 5e-06 z_p_prime_base_77[3] + [
   - alpha[2] * z_p_prime_base_77[3] + p[2] * z_p_prime_base_77[3] ] <= 0
 infeasibility_p_prime_base_77_row_4_leq_0: 5e-06 z_p_prime_base_77[4] + [
   0.5 alpha[0] * z_p_prime_base_77[4]
   + 0.5 alpha[1] * z_p_prime_base_77[4]
   + 0.5 alpha[2] * z_p_prime_base_77[4]
   - 0.5 alpha[3] * z_p_prime_base_77[4]
   + 0.5 alpha[4] * z_p_prime_base_77[4] + 0.5 p[0] * z_p_prime_base_77[4]
   - 0.5 p[1] * z_p_prime_base_77[4] - 0.5 p[2] * z_p_prime_base_77[4]
   + 0.5 p[3] * z_p_prime_base_77[4] + 0.5 p[4] * z_p_prime_base_77[4] ]
   <= 0
 infeasibility_p_base_78_row_0_leq_0: 5e-06 z_p_base_78[0] + [
   p[0] * z_p_base_78[0] ] <= 0
 infeasibility_p_base_78_row_1_leq_0: 5e-06 z_p_base_78[1] + [
   - 0.5 p[0] * z_p_base_78[1] + 0.5 p[1] * z_p_base_78[1]
   + 0.5 p[2] * z_p_base_78[1] + 0.5 p[3] * z_p_base_78[1]
   - 0.5 p[4] * z_p_base_78[1] ] <= 0
 infeasibility_p_base_78_row_2_leq_0: 5e-06 z_p_base_78[2] + [
   - 0.5 p[0] * z_p_base_78[2] + 0.5 p[1] * z_p_base_78[2]
   - 0.5 p[2] * z_p_base_78[2] - 0.5 p[3] * z_p_base_78[2]
   + 0.5 p[4] * z_p_base_78[2] ] <= 0
 infeasibility_p_base_78_row_3_leq_0: 5e-06 z_p_base_78[3] + [
   p[2] * z_p_base_78[3] ] <= 0
 infeasibility_p_base_78_row_4_leq_0: 5e-06 z_p_base_78[4] + [
   0.5 p[0] * z_p_base_78[4] - 0.5 p[1] * z_p_base_78[4]
   - 0.5 p[2] * z_p_base_78[4] + 0.5 p[3] * z_p_base_78[4]
   + 0.5 p[4] * z_p_base_78[4] ] <= 0
 infeasibility_p_prime_base_78_row_0_leq_0: 5e-06 z_p_prime_base_78[0] + [
   alpha[0] * z_p_prime_base_78[0] + p[0] * z_p_prime_base_78[0] ] <= 0
 infeasibility_p_prime_base_78_row_1_leq_0: 5e-06 z_p_prime_base_78[1] + [
   - 0.5 alpha[0] * z_p_prime_base_78[1]
   - 0.5 alpha[1] * z_p_prime_base_78[1]
   - 0.5 alpha[2] * z_p_prime_base_78[1]
   - 0.5 alpha[3] * z_p_prime_base_78[1]
   - 0.5 alpha[4] * z_p_prime_base_78[1] - 0.5 p[0] * z_p_prime_base_78[1]
   + 0.5 p[1] * z_p_prime_base_78[1] + 0.5 p[2] * z_p_prime_base_78[1]
   + 0.5 p[3] * z_p_prime_base_78[1] - 0.5 p[4] * z_p_prime_base_78[1] ]
   <= 0
 infeasibility_p_prime_base_78_row_2_leq_0: 5e-06 z_p_prime_base_78[2] + [
   - 0.5 alpha[0] * z_p_prime_base_78[2]
   - 0.5 alpha[1] * z_p_prime_base_78[2]
   + 0.5 alpha[2] * z_p_prime_base_78[2]
   + 0.5 alpha[3] * z_p_prime_base_78[2]
   + 0.5 alpha[4] * z_p_prime_base_78[2] - 0.5 p[0] * z_p_prime_base_78[2]
   + 0.5 p[1] * z_p_prime_base_78[2] - 0.5 p[2] * z_p_prime_base_78[2]
   - 0.5 p[3] * z_p_prime_base_78[2] + 0.5 p[4] * z_p_prime_base_78[2] ]
   <= 0
 infeasibility_p_prime_base_78_row_3_leq_0: 5e-06 z_p_prime_base_78[3] + [
   - alpha[2] * z_p_prime_base_78[3] + p[2] * z_p_prime_base_78[3] ] <= 0
 infeasibility_p_prime_base_78_row_4_leq_0: 5e-06 z_p_prime_base_78[4] + [
   0.5 alpha[0] * z_p_prime_base_78[4]
   + 0.5 alpha[1] * z_p_prime_base_78[4]
   + 0.5 alpha[2] * z_p_prime_base_78[4]
   - 0.5 alpha[3] * z_p_prime_base_78[4]
   + 0.5 alpha[4] * z_p_prime_base_78[4] + 0.5 p[0] * z_p_prime_base_78[4]
   - 0.5 p[1] * z_p_prime_base_78[4] - 0.5 p[2] * z_p_prime_base_78[4]
   + 0.5 p[3] * z_p_prime_base_78[4] + 0.5 p[4] * z_p_prime_base_78[4] ]
   <= 0
 infeasibility_p_base_79_row_0_leq_0: 5e-06 z_p_base_79[0] + [
   p[0] * z_p_base_79[0] ] <= 0
 infeasibility_p_base_79_row_1_leq_0: 5e-06 z_p_base_79[1] + [
   - p[0] * z_p_base_79[1] + p[1] * z_p_base_79[1] ] <= 0
 infeasibility_p_base_79_row_2_leq_0: 5e-06 z_p_base_79[2] + [
   0.5 p[0] * z_p_base_79[2] - 0.5 p[1] * z_p_base_79[2]
   + 0.5 p[2] * z_p_base_79[2] + 0.5 p[3] * z_p_base_79[2]
   - 0.5 p[4] * z_p_base_79[2] ] <= 0
 infeasibility_p_base_79_row_3_leq_0: 5e-06 z_p_base_79[3] + [
   - 0.5 p[0] * z_p_base_79[3] + 0.5 p[1] * z_p_base_79[3]
   + 0.5 p[2] * z_p_base_79[3] - 0.5 p[3] * z_p_base_79[3]
   + 0.5 p[4] * z_p_base_79[3] ] <= 0
 infeasibility_p_base_79_row_4_leq_0: 5e-06 z_p_base_79[4] + [
   0.5 p[0] * z_p_base_79[4] - 0.5 p[1] * z_p_base_79[4]
   - 0.5 p[2] * z_p_base_79[4] + 0.5 p[3] * z_p_base_79[4]
   + 0.5 p[4] * z_p_base_79[4] ] <= 0
 infeasibility_p_prime_base_79_row_0_leq_0: 5e-06 z_p_prime_base_79[0] + [
   alpha[0] * z_p_prime_base_79[0] + p[0] * z_p_prime_base_79[0] ] <= 0
 infeasibility_p_prime_base_79_row_1_leq_0: 5e-06 z_p_prime_base_79[1] + [
   - alpha[0] * z_p_prime_base_79[1] - alpha[1] * z_p_prime_base_79[1]
   - p[0] * z_p_prime_base_79[1] + p[1] * z_p_prime_base_79[1] ] <= 0
 infeasibility_p_prime_base_79_row_2_leq_0: 5e-06 z_p_prime_base_79[2] + [
   0.5 alpha[0] * z_p_prime_base_79[2]
   + 0.5 alpha[1] * z_p_prime_base_79[2]
   - 0.5 alpha[2] * z_p_prime_base_79[2]
   - 0.5 alpha[3] * z_p_prime_base_79[2]
   - 0.5 alpha[4] * z_p_prime_base_79[2] + 0.5 p[0] * z_p_prime_base_79[2]
   - 0.5 p[1] * z_p_prime_base_79[2] + 0.5 p[2] * z_p_prime_base_79[2]
   + 0.5 p[3] * z_p_prime_base_79[2] - 0.5 p[4] * z_p_prime_base_79[2] ]
   <= 0
 infeasibility_p_prime_base_79_row_3_leq_0: 5e-06 z_p_prime_base_79[3] + [
   - 0.5 alpha[0] * z_p_prime_base_79[3]
   - 0.5 alpha[1] * z_p_prime_base_79[3]
   - 0.5 alpha[2] * z_p_prime_base_79[3]
   + 0.5 alpha[3] * z_p_prime_base_79[3]
   + 0.5 alpha[4] * z_p_prime_base_79[3] - 0.5 p[0] * z_p_prime_base_79[3]
   + 0.5 p[1] * z_p_prime_base_79[3] + 0.5 p[2] * z_p_prime_base_79[3]
   - 0.5 p[3] * z_p_prime_base_79[3] + 0.5 p[4] * z_p_prime_base_79[3] ]
   <= 0
 infeasibility_p_prime_base_79_row_4_leq_0: 5e-06 z_p_prime_base_79[4] + [
   0.5 alpha[0] * z_p_prime_base_79[4]
   + 0.5 alpha[1] * z_p_prime_base_79[4]
   + 0.5 alpha[2] * z_p_prime_base_79[4]
   - 0.5 alpha[3] * z_p_prime_base_79[4]
   + 0.5 alpha[4] * z_p_prime_base_79[4] + 0.5 p[0] * z_p_prime_base_79[4]
   - 0.5 p[1] * z_p_prime_base_79[4] - 0.5 p[2] * z_p_prime_base_79[4]
   + 0.5 p[3] * z_p_prime_base_79[4] + 0.5 p[4] * z_p_prime_base_79[4] ]
   <= 0
 infeasibility_p_base_80_row_0_leq_0: 5e-06 z_p_base_80[0] + [
   p[0] * z_p_base_80[0] ] <= 0
 infeasibility_p_base_80_row_1_leq_0: 5e-06 z_p_base_80[1] + [
   - p[0] * z_p_base_80[1] + p[1] * z_p_base_80[1] ] <= 0
 infeasibility_p_base_80_row_2_leq_0: 5e-06 z_p_base_80[2] + [
   - 0.5 p[0] * z_p_base_80[2] + 0.5 p[1] * z_p_base_80[2]
   + 0.5 p[2] * z_p_base_80[2] + 0.5 p[3] * z_p_base_80[2]
   - 0.5 p[4] * z_p_base_80[2] ] <= 0
 infeasibility_p_base_80_row_3_leq_0: 5e-06 z_p_base_80[3] + [
   0.5 p[0] * z_p_base_80[3] - 0.5 p[1] * z_p_base_80[3]
   + 0.5 p[2] * z_p_base_80[3] - 0.5 p[3] * z_p_base_80[3]
   + 0.5 p[4] * z_p_base_80[3] ] <= 0
 infeasibility_p_base_80_row_4_leq_0: 5e-06 z_p_base_80[4] + [
   0.5 p[0] * z_p_base_80[4] - 0.5 p[1] * z_p_base_80[4]
   - 0.5 p[2] * z_p_base_80[4] + 0.5 p[3] * z_p_base_80[4]
   + 0.5 p[4] * z_p_base_80[4] ] <= 0
 infeasibility_p_prime_base_80_row_0_leq_0: 5e-06 z_p_prime_base_80[0] + [
   alpha[0] * z_p_prime_base_80[0] + p[0] * z_p_prime_base_80[0] ] <= 0
 infeasibility_p_prime_base_80_row_1_leq_0: 5e-06 z_p_prime_base_80[1] + [
   - alpha[0] * z_p_prime_base_80[1] - alpha[1] * z_p_prime_base_80[1]
   - p[0] * z_p_prime_base_80[1] + p[1] * z_p_prime_base_80[1] ] <= 0
 infeasibility_p_prime_base_80_row_2_leq_0: 5e-06 z_p_prime_base_80[2] + [
   - 0.5 alpha[0] * z_p_prime_base_80[2]
   - 0.5 alpha[1] * z_p_prime_base_80[2]
   - 0.5 alpha[2] * z_p_prime_base_80[2]
   - 0.5 alpha[3] * z_p_prime_base_80[2]
   - 0.5 alpha[4] * z_p_prime_base_80[2] - 0.5 p[0] * z_p_prime_base_80[2]
   + 0.5 p[1] * z_p_prime_base_80[2] + 0.5 p[2] * z_p_prime_base_80[2]
   + 0.5 p[3] * z_p_prime_base_80[2] - 0.5 p[4] * z_p_prime_base_80[2] ]
   <= 0
 infeasibility_p_prime_base_80_row_3_leq_0: 5e-06 z_p_prime_base_80[3] + [
   0.5 alpha[0] * z_p_prime_base_80[3]
   + 0.5 alpha[1] * z_p_prime_base_80[3]
   - 0.5 alpha[2] * z_p_prime_base_80[3]
   + 0.5 alpha[3] * z_p_prime_base_80[3]
   + 0.5 alpha[4] * z_p_prime_base_80[3] + 0.5 p[0] * z_p_prime_base_80[3]
   - 0.5 p[1] * z_p_prime_base_80[3] + 0.5 p[2] * z_p_prime_base_80[3]
   - 0.5 p[3] * z_p_prime_base_80[3] + 0.5 p[4] * z_p_prime_base_80[3] ]
   <= 0
 infeasibility_p_prime_base_80_row_4_leq_0: 5e-06 z_p_prime_base_80[4] + [
   0.5 alpha[0] * z_p_prime_base_80[4]
   + 0.5 alpha[1] * z_p_prime_base_80[4]
   + 0.5 alpha[2] * z_p_prime_base_80[4]
   - 0.5 alpha[3] * z_p_prime_base_80[4]
   + 0.5 alpha[4] * z_p_prime_base_80[4] + 0.5 p[0] * z_p_prime_base_80[4]
   - 0.5 p[1] * z_p_prime_base_80[4] - 0.5 p[2] * z_p_prime_base_80[4]
   + 0.5 p[3] * z_p_prime_base_80[4] + 0.5 p[4] * z_p_prime_base_80[4] ]
   <= 0
 infeasibility_p_base_81_row_0_leq_0: 5e-06 z_p_base_81[0] + [
   0.5 p[0] * z_p_base_81[0] - 0.5 p[1] * z_p_base_81[0]
   + 0.5 p[2] * z_p_base_81[0] - 0.5 p[3] * z_p_base_81[0]
   - 0.5 p[4] * z_p_base_81[0] ] <= 0
 infeasibility_p_base_81_row_1_leq_0: 5e-06 z_p_base_81[1] + [
   0.5 p[0] * z_p_base_81[1] + 0.5 p[1] * z_p_base_81[1]
   - 0.5 p[2] * z_p_base_81[1] + 0.5 p[3] * z_p_base_81[1]
   - 0.5 p[4] * z_p_base_81[1] ] <= 0
 infeasibility_p_base_81_row_2_leq_0: 5e-06 z_p_base_81[2] + [
   p[4] * z_p_base_81[2] ] <= 0
 infeasibility_p_base_81_row_3_leq_0: 5e-06 z_p_base_81[3] + [
   p[1] * z_p_base_81[3] ] <= 0
 infeasibility_p_base_81_row_4_leq_0: 5e-06 z_p_base_81[4] + [
   - 0.5 p[0] * z_p_base_81[4] - 0.5 p[1] * z_p_base_81[4]
   + 0.5 p[2] * z_p_base_81[4] + 0.5 p[3] * z_p_base_81[4]
   + 0.5 p[4] * z_p_base_81[4] ] <= 0
 infeasibility_p_prime_base_81_row_0_leq_0: 5e-06 z_p_prime_base_81[0] + [
   0.5 alpha[0] * z_p_prime_base_81[0]
   + 0.5 alpha[1] * z_p_prime_base_81[0]
   - 0.5 alpha[2] * z_p_prime_base_81[0]
   + 0.5 alpha[3] * z_p_prime_base_81[0]
   - 0.5 alpha[4] * z_p_prime_base_81[0] + 0.5 p[0] * z_p_prime_base_81[0]
   - 0.5 p[1] * z_p_prime_base_81[0] + 0.5 p[2] * z_p_prime_base_81[0]
   - 0.5 p[3] * z_p_prime_base_81[0] - 0.5 p[4] * z_p_prime_base_81[0] ]
   <= 0
 infeasibility_p_prime_base_81_row_1_leq_0: 5e-06 z_p_prime_base_81[1] + [
   0.5 alpha[0] * z_p_prime_base_81[1]
   - 0.5 alpha[1] * z_p_prime_base_81[1]
   + 0.5 alpha[2] * z_p_prime_base_81[1]
   - 0.5 alpha[3] * z_p_prime_base_81[1]
   - 0.5 alpha[4] * z_p_prime_base_81[1] + 0.5 p[0] * z_p_prime_base_81[1]
   + 0.5 p[1] * z_p_prime_base_81[1] - 0.5 p[2] * z_p_prime_base_81[1]
   + 0.5 p[3] * z_p_prime_base_81[1] - 0.5 p[4] * z_p_prime_base_81[1] ]
   <= 0
 infeasibility_p_prime_base_81_row_2_leq_0: 5e-06 z_p_prime_base_81[2] + [
   alpha[4] * z_p_prime_base_81[2] + p[4] * z_p_prime_base_81[2] ] <= 0
 infeasibility_p_prime_base_81_row_3_leq_0: 5e-06 z_p_prime_base_81[3] + [
   - alpha[1] * z_p_prime_base_81[3] + p[1] * z_p_prime_base_81[3] ] <= 0
 infeasibility_p_prime_base_81_row_4_leq_0: 5e-06 z_p_prime_base_81[4] + [
   - 0.5 alpha[0] * z_p_prime_base_81[4]
   + 0.5 alpha[1] * z_p_prime_base_81[4]
   - 0.5 alpha[2] * z_p_prime_base_81[4]
   - 0.5 alpha[3] * z_p_prime_base_81[4]
   + 0.5 alpha[4] * z_p_prime_base_81[4] - 0.5 p[0] * z_p_prime_base_81[4]
   - 0.5 p[1] * z_p_prime_base_81[4] + 0.5 p[2] * z_p_prime_base_81[4]
   + 0.5 p[3] * z_p_prime_base_81[4] + 0.5 p[4] * z_p_prime_base_81[4] ]
   <= 0
 infeasibility_p_base_82_row_0_leq_0: 5e-06 z_p_base_82[0] + [
   0.5 p[0] * z_p_base_82[0] - 0.5 p[1] * z_p_base_82[0]
   + 0.5 p[2] * z_p_base_82[0] - 0.5 p[3] * z_p_base_82[0]
   - 0.5 p[4] * z_p_base_82[0] ] <= 0
 infeasibility_p_base_82_row_1_leq_0: 5e-06 z_p_base_82[1] + [
   p[3] * z_p_base_82[1] ] <= 0
 infeasibility_p_base_82_row_2_leq_0: 5e-06 z_p_base_82[2] + [
   0.5 p[0] * z_p_base_82[2] + 0.5 p[1] * z_p_base_82[2]
   - 0.5 p[2] * z_p_base_82[2] - 0.5 p[3] * z_p_base_82[2]
   + 0.5 p[4] * z_p_base_82[2] ] <= 0
 infeasibility_p_base_82_row_3_leq_0: 5e-06 z_p_base_82[3] + [
   p[1] * z_p_base_82[3] ] <= 0
 infeasibility_p_base_82_row_4_leq_0: 5e-06 z_p_base_82[4] + [
   - 0.5 p[0] * z_p_base_82[4] - 0.5 p[1] * z_p_base_82[4]
   + 0.5 p[2] * z_p_base_82[4] + 0.5 p[3] * z_p_base_82[4]
   + 0.5 p[4] * z_p_base_82[4] ] <= 0
 infeasibility_p_prime_base_82_row_0_leq_0: 5e-06 z_p_prime_base_82[0] + [
   0.5 alpha[0] * z_p_prime_base_82[0]
   + 0.5 alpha[1] * z_p_prime_base_82[0]
   - 0.5 alpha[2] * z_p_prime_base_82[0]
   + 0.5 alpha[3] * z_p_prime_base_82[0]
   - 0.5 alpha[4] * z_p_prime_base_82[0] + 0.5 p[0] * z_p_prime_base_82[0]
   - 0.5 p[1] * z_p_prime_base_82[0] + 0.5 p[2] * z_p_prime_base_82[0]
   - 0.5 p[3] * z_p_prime_base_82[0] - 0.5 p[4] * z_p_prime_base_82[0] ]
   <= 0
 infeasibility_p_prime_base_82_row_1_leq_0: 5e-06 z_p_prime_base_82[1] + [
   - alpha[3] * z_p_prime_base_82[1] + p[3] * z_p_prime_base_82[1] ] <= 0
 infeasibility_p_prime_base_82_row_2_leq_0: 5e-06 z_p_prime_base_82[2] + [
   0.5 alpha[0] * z_p_prime_base_82[2]
   - 0.5 alpha[1] * z_p_prime_base_82[2]
   + 0.5 alpha[2] * z_p_prime_base_82[2]
   + 0.5 alpha[3] * z_p_prime_base_82[2]
   + 0.5 alpha[4] * z_p_prime_base_82[2] + 0.5 p[0] * z_p_prime_base_82[2]
   + 0.5 p[1] * z_p_prime_base_82[2] - 0.5 p[2] * z_p_prime_base_82[2]
   - 0.5 p[3] * z_p_prime_base_82[2] + 0.5 p[4] * z_p_prime_base_82[2] ]
   <= 0
 infeasibility_p_prime_base_82_row_3_leq_0: 5e-06 z_p_prime_base_82[3] + [
   - alpha[1] * z_p_prime_base_82[3] + p[1] * z_p_prime_base_82[3] ] <= 0
 infeasibility_p_prime_base_82_row_4_leq_0: 5e-06 z_p_prime_base_82[4] + [
   - 0.5 alpha[0] * z_p_prime_base_82[4]
   + 0.5 alpha[1] * z_p_prime_base_82[4]
   - 0.5 alpha[2] * z_p_prime_base_82[4]
   - 0.5 alpha[3] * z_p_prime_base_82[4]
   + 0.5 alpha[4] * z_p_prime_base_82[4] - 0.5 p[0] * z_p_prime_base_82[4]
   - 0.5 p[1] * z_p_prime_base_82[4] + 0.5 p[2] * z_p_prime_base_82[4]
   + 0.5 p[3] * z_p_prime_base_82[4] + 0.5 p[4] * z_p_prime_base_82[4] ]
   <= 0
 infeasibility_p_base_83_row_0_leq_0: 5e-06 z_p_base_83[0] + [
   - p[1] * z_p_base_83[0] + p[2] * z_p_base_83[0] ] <= 0
 infeasibility_p_base_83_row_1_leq_0: 5e-06 z_p_base_83[1] + [
   0.5 p[0] * z_p_base_83[1] + 0.5 p[1] * z_p_base_83[1]
   - 0.5 p[2] * z_p_base_83[1] + 0.5 p[3] * z_p_base_83[1]
   - 0.5 p[4] * z_p_base_83[1] ] <= 0
 infeasibility_p_base_83_row_2_leq_0: 5e-06 z_p_base_83[2] + [
   0.5 p[0] * z_p_base_83[2] + 0.5 p[1] * z_p_base_83[2]
   - 0.5 p[2] * z_p_base_83[2] - 0.5 p[3] * z_p_base_83[2]
   + 0.5 p[4] * z_p_base_83[2] ] <= 0
 infeasibility_p_base_83_row_3_leq_0: 5e-06 z_p_base_83[3] + [
   p[1] * z_p_base_83[3] ] <= 0
 infeasibility_p_base_83_row_4_leq_0: 5e-06 z_p_base_83[4] + [
   - 0.5 p[0] * z_p_base_83[4] - 0.5 p[1] * z_p_base_83[4]
   + 0.5 p[2] * z_p_base_83[4] + 0.5 p[3] * z_p_base_83[4]
   + 0.5 p[4] * z_p_base_83[4] ] <= 0
 infeasibility_p_prime_base_83_row_0_leq_0: 5e-06 z_p_prime_base_83[0] + [
   alpha[1] * z_p_prime_base_83[0] - alpha[2] * z_p_prime_base_83[0]
   - p[1] * z_p_prime_base_83[0] + p[2] * z_p_prime_base_83[0] ] <= 0
 infeasibility_p_prime_base_83_row_1_leq_0: 5e-06 z_p_prime_base_83[1] + [
   0.5 alpha[0] * z_p_prime_base_83[1]
   - 0.5 alpha[1] * z_p_prime_base_83[1]
   + 0.5 alpha[2] * z_p_prime_base_83[1]
   - 0.5 alpha[3] * z_p_prime_base_83[1]
   - 0.5 alpha[4] * z_p_prime_base_83[1] + 0.5 p[0] * z_p_prime_base_83[1]
   + 0.5 p[1] * z_p_prime_base_83[1] - 0.5 p[2] * z_p_prime_base_83[1]
   + 0.5 p[3] * z_p_prime_base_83[1] - 0.5 p[4] * z_p_prime_base_83[1] ]
   <= 0
 infeasibility_p_prime_base_83_row_2_leq_0: 5e-06 z_p_prime_base_83[2] + [
   0.5 alpha[0] * z_p_prime_base_83[2]
   - 0.5 alpha[1] * z_p_prime_base_83[2]
   + 0.5 alpha[2] * z_p_prime_base_83[2]
   + 0.5 alpha[3] * z_p_prime_base_83[2]
   + 0.5 alpha[4] * z_p_prime_base_83[2] + 0.5 p[0] * z_p_prime_base_83[2]
   + 0.5 p[1] * z_p_prime_base_83[2] - 0.5 p[2] * z_p_prime_base_83[2]
   - 0.5 p[3] * z_p_prime_base_83[2] + 0.5 p[4] * z_p_prime_base_83[2] ]
   <= 0
 infeasibility_p_prime_base_83_row_3_leq_0: 5e-06 z_p_prime_base_83[3] + [
   - alpha[1] * z_p_prime_base_83[3] + p[1] * z_p_prime_base_83[3] ] <= 0
 infeasibility_p_prime_base_83_row_4_leq_0: 5e-06 z_p_prime_base_83[4] + [
   - 0.5 alpha[0] * z_p_prime_base_83[4]
   + 0.5 alpha[1] * z_p_prime_base_83[4]
   - 0.5 alpha[2] * z_p_prime_base_83[4]
   - 0.5 alpha[3] * z_p_prime_base_83[4]
   + 0.5 alpha[4] * z_p_prime_base_83[4] - 0.5 p[0] * z_p_prime_base_83[4]
   - 0.5 p[1] * z_p_prime_base_83[4] + 0.5 p[2] * z_p_prime_base_83[4]
   + 0.5 p[3] * z_p_prime_base_83[4] + 0.5 p[4] * z_p_prime_base_83[4] ]
   <= 0
 infeasibility_p_base_84_row_0_leq_0: 5e-06 z_p_base_84[0] + [
   0.5 p[0] * z_p_base_84[0] + 0.5 p[1] * z_p_base_84[0]
   + 0.5 p[2] * z_p_base_84[0] - 0.5 p[3] * z_p_base_84[0]
   - 0.5 p[4] * z_p_base_84[0] ] <= 0
 infeasibility_p_base_84_row_1_leq_0: 5e-06 z_p_base_84[1] + [
   0.5 p[0] * z_p_base_84[1] - 0.5 p[1] * z_p_base_84[1]
   - 0.5 p[2] * z_p_base_84[1] + 0.5 p[3] * z_p_base_84[1]
   - 0.5 p[4] * z_p_base_84[1] ] <= 0
 infeasibility_p_base_84_row_2_leq_0: 5e-06 z_p_base_84[2] + [
   p[4] * z_p_base_84[2] ] <= 0
 infeasibility_p_base_84_row_3_leq_0: 5e-06 z_p_base_84[3] + [
   p[1] * z_p_base_84[3] ] <= 0
 infeasibility_p_base_84_row_4_leq_0: 5e-06 z_p_base_84[4] + [
   - 0.5 p[0] * z_p_base_84[4] - 0.5 p[1] * z_p_base_84[4]
   + 0.5 p[2] * z_p_base_84[4] + 0.5 p[3] * z_p_base_84[4]
   + 0.5 p[4] * z_p_base_84[4] ] <= 0
 infeasibility_p_prime_base_84_row_0_leq_0: 5e-06 z_p_prime_base_84[0] + [
   0.5 alpha[0] * z_p_prime_base_84[0]
   - 0.5 alpha[1] * z_p_prime_base_84[0]
   - 0.5 alpha[2] * z_p_prime_base_84[0]
   + 0.5 alpha[3] * z_p_prime_base_84[0]
   - 0.5 alpha[4] * z_p_prime_base_84[0] + 0.5 p[0] * z_p_prime_base_84[0]
   + 0.5 p[1] * z_p_prime_base_84[0] + 0.5 p[2] * z_p_prime_base_84[0]
   - 0.5 p[3] * z_p_prime_base_84[0] - 0.5 p[4] * z_p_prime_base_84[0] ]
   <= 0
 infeasibility_p_prime_base_84_row_1_leq_0: 5e-06 z_p_prime_base_84[1] + [
   0.5 alpha[0] * z_p_prime_base_84[1]
   + 0.5 alpha[1] * z_p_prime_base_84[1]
   + 0.5 alpha[2] * z_p_prime_base_84[1]
   - 0.5 alpha[3] * z_p_prime_base_84[1]
   - 0.5 alpha[4] * z_p_prime_base_84[1] + 0.5 p[0] * z_p_prime_base_84[1]
   - 0.5 p[1] * z_p_prime_base_84[1] - 0.5 p[2] * z_p_prime_base_84[1]
   + 0.5 p[3] * z_p_prime_base_84[1] - 0.5 p[4] * z_p_prime_base_84[1] ]
   <= 0
 infeasibility_p_prime_base_84_row_2_leq_0: 5e-06 z_p_prime_base_84[2] + [
   alpha[4] * z_p_prime_base_84[2] + p[4] * z_p_prime_base_84[2] ] <= 0
 infeasibility_p_prime_base_84_row_3_leq_0: 5e-06 z_p_prime_base_84[3] + [
   - alpha[1] * z_p_prime_base_84[3] + p[1] * z_p_prime_base_84[3] ] <= 0
 infeasibility_p_prime_base_84_row_4_leq_0: 5e-06 z_p_prime_base_84[4] + [
   - 0.5 alpha[0] * z_p_prime_base_84[4]
   + 0.5 alpha[1] * z_p_prime_base_84[4]
   - 0.5 alpha[2] * z_p_prime_base_84[4]
   - 0.5 alpha[3] * z_p_prime_base_84[4]
   + 0.5 alpha[4] * z_p_prime_base_84[4] - 0.5 p[0] * z_p_prime_base_84[4]
   - 0.5 p[1] * z_p_prime_base_84[4] + 0.5 p[2] * z_p_prime_base_84[4]
   + 0.5 p[3] * z_p_prime_base_84[4] + 0.5 p[4] * z_p_prime_base_84[4] ]
   <= 0
 infeasibility_p_base_85_row_0_leq_0: 5e-06 z_p_base_85[0] + [
   0.5 p[0] * z_p_base_85[0] + 0.5 p[1] * z_p_base_85[0]
   + 0.5 p[2] * z_p_base_85[0] - 0.5 p[3] * z_p_base_85[0]
   - 0.5 p[4] * z_p_base_85[0] ] <= 0
 infeasibility_p_base_85_row_1_leq_0: 5e-06 z_p_base_85[1] + [
   - p[1] * z_p_base_85[1] + p[3] * z_p_base_85[1] ] <= 0
 infeasibility_p_base_85_row_2_leq_0: 5e-06 z_p_base_85[2] + [
   0.5 p[0] * z_p_base_85[2] + 0.5 p[1] * z_p_base_85[2]
   - 0.5 p[2] * z_p_base_85[2] - 0.5 p[3] * z_p_base_85[2]
   + 0.5 p[4] * z_p_base_85[2] ] <= 0
 infeasibility_p_base_85_row_3_leq_0: 5e-06 z_p_base_85[3] + [
   p[1] * z_p_base_85[3] ] <= 0
 infeasibility_p_base_85_row_4_leq_0: 5e-06 z_p_base_85[4] + [
   - 0.5 p[0] * z_p_base_85[4] - 0.5 p[1] * z_p_base_85[4]
   + 0.5 p[2] * z_p_base_85[4] + 0.5 p[3] * z_p_base_85[4]
   + 0.5 p[4] * z_p_base_85[4] ] <= 0
 infeasibility_p_prime_base_85_row_0_leq_0: 5e-06 z_p_prime_base_85[0] + [
   0.5 alpha[0] * z_p_prime_base_85[0]
   - 0.5 alpha[1] * z_p_prime_base_85[0]
   - 0.5 alpha[2] * z_p_prime_base_85[0]
   + 0.5 alpha[3] * z_p_prime_base_85[0]
   - 0.5 alpha[4] * z_p_prime_base_85[0] + 0.5 p[0] * z_p_prime_base_85[0]
   + 0.5 p[1] * z_p_prime_base_85[0] + 0.5 p[2] * z_p_prime_base_85[0]
   - 0.5 p[3] * z_p_prime_base_85[0] - 0.5 p[4] * z_p_prime_base_85[0] ]
   <= 0
 infeasibility_p_prime_base_85_row_1_leq_0: 5e-06 z_p_prime_base_85[1] + [
   alpha[1] * z_p_prime_base_85[1] - alpha[3] * z_p_prime_base_85[1]
   - p[1] * z_p_prime_base_85[1] + p[3] * z_p_prime_base_85[1] ] <= 0
 infeasibility_p_prime_base_85_row_2_leq_0: 5e-06 z_p_prime_base_85[2] + [
   0.5 alpha[0] * z_p_prime_base_85[2]
   - 0.5 alpha[1] * z_p_prime_base_85[2]
   + 0.5 alpha[2] * z_p_prime_base_85[2]
   + 0.5 alpha[3] * z_p_prime_base_85[2]
   + 0.5 alpha[4] * z_p_prime_base_85[2] + 0.5 p[0] * z_p_prime_base_85[2]
   + 0.5 p[1] * z_p_prime_base_85[2] - 0.5 p[2] * z_p_prime_base_85[2]
   - 0.5 p[3] * z_p_prime_base_85[2] + 0.5 p[4] * z_p_prime_base_85[2] ]
   <= 0
 infeasibility_p_prime_base_85_row_3_leq_0: 5e-06 z_p_prime_base_85[3] + [
   - alpha[1] * z_p_prime_base_85[3] + p[1] * z_p_prime_base_85[3] ] <= 0
 infeasibility_p_prime_base_85_row_4_leq_0: 5e-06 z_p_prime_base_85[4] + [
   - 0.5 alpha[0] * z_p_prime_base_85[4]
   + 0.5 alpha[1] * z_p_prime_base_85[4]
   - 0.5 alpha[2] * z_p_prime_base_85[4]
   - 0.5 alpha[3] * z_p_prime_base_85[4]
   + 0.5 alpha[4] * z_p_prime_base_85[4] - 0.5 p[0] * z_p_prime_base_85[4]
   - 0.5 p[1] * z_p_prime_base_85[4] + 0.5 p[2] * z_p_prime_base_85[4]
   + 0.5 p[3] * z_p_prime_base_85[4] + 0.5 p[4] * z_p_prime_base_85[4] ]
   <= 0
 infeasibility_p_base_86_row_0_leq_0: 5e-06 z_p_base_86[0] + [
   p[2] * z_p_base_86[0] ] <= 0
 infeasibility_p_base_86_row_1_leq_0: 5e-06 z_p_base_86[1] + [
   0.5 p[0] * z_p_base_86[1] - 0.5 p[1] * z_p_base_86[1]
   - 0.5 p[2] * z_p_base_86[1] + 0.5 p[3] * z_p_base_86[1]
   - 0.5 p[4] * z_p_base_86[1] ] <= 0
 infeasibility_p_base_86_row_2_leq_0: 5e-06 z_p_base_86[2] + [
   0.5 p[0] * z_p_base_86[2] + 0.5 p[1] * z_p_base_86[2]
   - 0.5 p[2] * z_p_base_86[2] - 0.5 p[3] * z_p_base_86[2]
   + 0.5 p[4] * z_p_base_86[2] ] <= 0
 infeasibility_p_base_86_row_3_leq_0: 5e-06 z_p_base_86[3] + [
   p[1] * z_p_base_86[3] ] <= 0
 infeasibility_p_base_86_row_4_leq_0: 5e-06 z_p_base_86[4] + [
   - 0.5 p[0] * z_p_base_86[4] - 0.5 p[1] * z_p_base_86[4]
   + 0.5 p[2] * z_p_base_86[4] + 0.5 p[3] * z_p_base_86[4]
   + 0.5 p[4] * z_p_base_86[4] ] <= 0
 infeasibility_p_prime_base_86_row_0_leq_0: 5e-06 z_p_prime_base_86[0] + [
   - alpha[2] * z_p_prime_base_86[0] + p[2] * z_p_prime_base_86[0] ] <= 0
 infeasibility_p_prime_base_86_row_1_leq_0: 5e-06 z_p_prime_base_86[1] + [
   0.5 alpha[0] * z_p_prime_base_86[1]
   + 0.5 alpha[1] * z_p_prime_base_86[1]
   + 0.5 alpha[2] * z_p_prime_base_86[1]
   - 0.5 alpha[3] * z_p_prime_base_86[1]
   - 0.5 alpha[4] * z_p_prime_base_86[1] + 0.5 p[0] * z_p_prime_base_86[1]
   - 0.5 p[1] * z_p_prime_base_86[1] - 0.5 p[2] * z_p_prime_base_86[1]
   + 0.5 p[3] * z_p_prime_base_86[1] - 0.5 p[4] * z_p_prime_base_86[1] ]
   <= 0
 infeasibility_p_prime_base_86_row_2_leq_0: 5e-06 z_p_prime_base_86[2] + [
   0.5 alpha[0] * z_p_prime_base_86[2]
   - 0.5 alpha[1] * z_p_prime_base_86[2]
   + 0.5 alpha[2] * z_p_prime_base_86[2]
   + 0.5 alpha[3] * z_p_prime_base_86[2]
   + 0.5 alpha[4] * z_p_prime_base_86[2] + 0.5 p[0] * z_p_prime_base_86[2]
   + 0.5 p[1] * z_p_prime_base_86[2] - 0.5 p[2] * z_p_prime_base_86[2]
   - 0.5 p[3] * z_p_prime_base_86[2] + 0.5 p[4] * z_p_prime_base_86[2] ]
   <= 0
 infeasibility_p_prime_base_86_row_3_leq_0: 5e-06 z_p_prime_base_86[3] + [
   - alpha[1] * z_p_prime_base_86[3] + p[1] * z_p_prime_base_86[3] ] <= 0
 infeasibility_p_prime_base_86_row_4_leq_0: 5e-06 z_p_prime_base_86[4] + [
   - 0.5 alpha[0] * z_p_prime_base_86[4]
   + 0.5 alpha[1] * z_p_prime_base_86[4]
   - 0.5 alpha[2] * z_p_prime_base_86[4]
   - 0.5 alpha[3] * z_p_prime_base_86[4]
   + 0.5 alpha[4] * z_p_prime_base_86[4] - 0.5 p[0] * z_p_prime_base_86[4]
   - 0.5 p[1] * z_p_prime_base_86[4] + 0.5 p[2] * z_p_prime_base_86[4]
   + 0.5 p[3] * z_p_prime_base_86[4] + 0.5 p[4] * z_p_prime_base_86[4] ]
   <= 0
 infeasibility_p_base_87_row_0_leq_0: 5e-06 z_p_base_87[0] + [
   0.5 p[0] * z_p_base_87[0] + 0.5 p[1] * z_p_base_87[0]
   + 0.5 p[2] * z_p_base_87[0] - 0.5 p[3] * z_p_base_87[0]
   - 0.5 p[4] * z_p_base_87[0] ] <= 0
 infeasibility_p_base_87_row_1_leq_0: 5e-06 z_p_base_87[1] + [
   0.5 p[0] * z_p_base_87[1] + 0.5 p[1] * z_p_base_87[1]
   - 0.5 p[2] * z_p_base_87[1] + 0.5 p[3] * z_p_base_87[1]
   - 0.5 p[4] * z_p_base_87[1] ] <= 0
 infeasibility_p_base_87_row_2_leq_0: 5e-06 z_p_base_87[2] + [
   - p[1] * z_p_base_87[2] + p[4] * z_p_base_87[2] ] <= 0
 infeasibility_p_base_87_row_3_leq_0: 5e-06 z_p_base_87[3] + [
   p[1] * z_p_base_87[3] ] <= 0
 infeasibility_p_base_87_row_4_leq_0: 5e-06 z_p_base_87[4] + [
   - 0.5 p[0] * z_p_base_87[4] - 0.5 p[1] * z_p_base_87[4]
   + 0.5 p[2] * z_p_base_87[4] + 0.5 p[3] * z_p_base_87[4]
   + 0.5 p[4] * z_p_base_87[4] ] <= 0
 infeasibility_p_prime_base_87_row_0_leq_0: 5e-06 z_p_prime_base_87[0] + [
   0.5 alpha[0] * z_p_prime_base_87[0]
   - 0.5 alpha[1] * z_p_prime_base_87[0]
   - 0.5 alpha[2] * z_p_prime_base_87[0]
   + 0.5 alpha[3] * z_p_prime_base_87[0]
   - 0.5 alpha[4] * z_p_prime_base_87[0] + 0.5 p[0] * z_p_prime_base_87[0]
   + 0.5 p[1] * z_p_prime_base_87[0] + 0.5 p[2] * z_p_prime_base_87[0]
   - 0.5 p[3] * z_p_prime_base_87[0] - 0.5 p[4] * z_p_prime_base_87[0] ]
   <= 0
 infeasibility_p_prime_base_87_row_1_leq_0: 5e-06 z_p_prime_base_87[1] + [
   0.5 alpha[0] * z_p_prime_base_87[1]
   - 0.5 alpha[1] * z_p_prime_base_87[1]
   + 0.5 alpha[2] * z_p_prime_base_87[1]
   - 0.5 alpha[3] * z_p_prime_base_87[1]
   - 0.5 alpha[4] * z_p_prime_base_87[1] + 0.5 p[0] * z_p_prime_base_87[1]
   + 0.5 p[1] * z_p_prime_base_87[1] - 0.5 p[2] * z_p_prime_base_87[1]
   + 0.5 p[3] * z_p_prime_base_87[1] - 0.5 p[4] * z_p_prime_base_87[1] ]
   <= 0
 infeasibility_p_prime_base_87_row_2_leq_0: 5e-06 z_p_prime_base_87[2] + [
   alpha[1] * z_p_prime_base_87[2] + alpha[4] * z_p_prime_base_87[2]
   - p[1] * z_p_prime_base_87[2] + p[4] * z_p_prime_base_87[2] ] <= 0
 infeasibility_p_prime_base_87_row_3_leq_0: 5e-06 z_p_prime_base_87[3] + [
   - alpha[1] * z_p_prime_base_87[3] + p[1] * z_p_prime_base_87[3] ] <= 0
 infeasibility_p_prime_base_87_row_4_leq_0: 5e-06 z_p_prime_base_87[4] + [
   - 0.5 alpha[0] * z_p_prime_base_87[4]
   + 0.5 alpha[1] * z_p_prime_base_87[4]
   - 0.5 alpha[2] * z_p_prime_base_87[4]
   - 0.5 alpha[3] * z_p_prime_base_87[4]
   + 0.5 alpha[4] * z_p_prime_base_87[4] - 0.5 p[0] * z_p_prime_base_87[4]
   - 0.5 p[1] * z_p_prime_base_87[4] + 0.5 p[2] * z_p_prime_base_87[4]
   + 0.5 p[3] * z_p_prime_base_87[4] + 0.5 p[4] * z_p_prime_base_87[4] ]
   <= 0
 infeasibility_p_base_88_row_0_leq_0: 5e-06 z_p_base_88[0] + [
   0.5 p[0] * z_p_base_88[0] + 0.5 p[1] * z_p_base_88[0]
   + 0.5 p[2] * z_p_base_88[0] - 0.5 p[3] * z_p_base_88[0]
   - 0.5 p[4] * z_p_base_88[0] ] <= 0
 infeasibility_p_base_88_row_1_leq_0: 5e-06 z_p_base_88[1] + [
   p[3] * z_p_base_88[1] ] <= 0
 infeasibility_p_base_88_row_2_leq_0: 5e-06 z_p_base_88[2] + [
   0.5 p[0] * z_p_base_88[2] - 0.5 p[1] * z_p_base_88[2]
   - 0.5 p[2] * z_p_base_88[2] - 0.5 p[3] * z_p_base_88[2]
   + 0.5 p[4] * z_p_base_88[2] ] <= 0
 infeasibility_p_base_88_row_3_leq_0: 5e-06 z_p_base_88[3] + [
   p[1] * z_p_base_88[3] ] <= 0
 infeasibility_p_base_88_row_4_leq_0: 5e-06 z_p_base_88[4] + [
   - 0.5 p[0] * z_p_base_88[4] - 0.5 p[1] * z_p_base_88[4]
   + 0.5 p[2] * z_p_base_88[4] + 0.5 p[3] * z_p_base_88[4]
   + 0.5 p[4] * z_p_base_88[4] ] <= 0
 infeasibility_p_prime_base_88_row_0_leq_0: 5e-06 z_p_prime_base_88[0] + [
   0.5 alpha[0] * z_p_prime_base_88[0]
   - 0.5 alpha[1] * z_p_prime_base_88[0]
   - 0.5 alpha[2] * z_p_prime_base_88[0]
   + 0.5 alpha[3] * z_p_prime_base_88[0]
   - 0.5 alpha[4] * z_p_prime_base_88[0] + 0.5 p[0] * z_p_prime_base_88[0]
   + 0.5 p[1] * z_p_prime_base_88[0] + 0.5 p[2] * z_p_prime_base_88[0]
   - 0.5 p[3] * z_p_prime_base_88[0] - 0.5 p[4] * z_p_prime_base_88[0] ]
   <= 0
 infeasibility_p_prime_base_88_row_1_leq_0: 5e-06 z_p_prime_base_88[1] + [
   - alpha[3] * z_p_prime_base_88[1] + p[3] * z_p_prime_base_88[1] ] <= 0
 infeasibility_p_prime_base_88_row_2_leq_0: 5e-06 z_p_prime_base_88[2] + [
   0.5 alpha[0] * z_p_prime_base_88[2]
   + 0.5 alpha[1] * z_p_prime_base_88[2]
   + 0.5 alpha[2] * z_p_prime_base_88[2]
   + 0.5 alpha[3] * z_p_prime_base_88[2]
   + 0.5 alpha[4] * z_p_prime_base_88[2] + 0.5 p[0] * z_p_prime_base_88[2]
   - 0.5 p[1] * z_p_prime_base_88[2] - 0.5 p[2] * z_p_prime_base_88[2]
   - 0.5 p[3] * z_p_prime_base_88[2] + 0.5 p[4] * z_p_prime_base_88[2] ]
   <= 0
 infeasibility_p_prime_base_88_row_3_leq_0: 5e-06 z_p_prime_base_88[3] + [
   - alpha[1] * z_p_prime_base_88[3] + p[1] * z_p_prime_base_88[3] ] <= 0
 infeasibility_p_prime_base_88_row_4_leq_0: 5e-06 z_p_prime_base_88[4] + [
   - 0.5 alpha[0] * z_p_prime_base_88[4]
   + 0.5 alpha[1] * z_p_prime_base_88[4]
   - 0.5 alpha[2] * z_p_prime_base_88[4]
   - 0.5 alpha[3] * z_p_prime_base_88[4]
   + 0.5 alpha[4] * z_p_prime_base_88[4] - 0.5 p[0] * z_p_prime_base_88[4]
   - 0.5 p[1] * z_p_prime_base_88[4] + 0.5 p[2] * z_p_prime_base_88[4]
   + 0.5 p[3] * z_p_prime_base_88[4] + 0.5 p[4] * z_p_prime_base_88[4] ]
   <= 0
 infeasibility_p_base_89_row_0_leq_0: 5e-06 z_p_base_89[0] + [
   p[2] * z_p_base_89[0] ] <= 0
 infeasibility_p_base_89_row_1_leq_0: 5e-06 z_p_base_89[1] + [
   0.5 p[0] * z_p_base_89[1] + 0.5 p[1] * z_p_base_89[1]
   - 0.5 p[2] * z_p_base_89[1] + 0.5 p[3] * z_p_base_89[1]
   - 0.5 p[4] * z_p_base_89[1] ] <= 0
 infeasibility_p_base_89_row_2_leq_0: 5e-06 z_p_base_89[2] + [
   0.5 p[0] * z_p_base_89[2] - 0.5 p[1] * z_p_base_89[2]
   - 0.5 p[2] * z_p_base_89[2] - 0.5 p[3] * z_p_base_89[2]
   + 0.5 p[4] * z_p_base_89[2] ] <= 0
 infeasibility_p_base_89_row_3_leq_0: 5e-06 z_p_base_89[3] + [
   p[1] * z_p_base_89[3] ] <= 0
 infeasibility_p_base_89_row_4_leq_0: 5e-06 z_p_base_89[4] + [
   - 0.5 p[0] * z_p_base_89[4] - 0.5 p[1] * z_p_base_89[4]
   + 0.5 p[2] * z_p_base_89[4] + 0.5 p[3] * z_p_base_89[4]
   + 0.5 p[4] * z_p_base_89[4] ] <= 0
 infeasibility_p_prime_base_89_row_0_leq_0: 5e-06 z_p_prime_base_89[0] + [
   - alpha[2] * z_p_prime_base_89[0] + p[2] * z_p_prime_base_89[0] ] <= 0
 infeasibility_p_prime_base_89_row_1_leq_0: 5e-06 z_p_prime_base_89[1] + [
   0.5 alpha[0] * z_p_prime_base_89[1]
   - 0.5 alpha[1] * z_p_prime_base_89[1]
   + 0.5 alpha[2] * z_p_prime_base_89[1]
   - 0.5 alpha[3] * z_p_prime_base_89[1]
   - 0.5 alpha[4] * z_p_prime_base_89[1] + 0.5 p[0] * z_p_prime_base_89[1]
   + 0.5 p[1] * z_p_prime_base_89[1] - 0.5 p[2] * z_p_prime_base_89[1]
   + 0.5 p[3] * z_p_prime_base_89[1] - 0.5 p[4] * z_p_prime_base_89[1] ]
   <= 0
 infeasibility_p_prime_base_89_row_2_leq_0: 5e-06 z_p_prime_base_89[2] + [
   0.5 alpha[0] * z_p_prime_base_89[2]
   + 0.5 alpha[1] * z_p_prime_base_89[2]
   + 0.5 alpha[2] * z_p_prime_base_89[2]
   + 0.5 alpha[3] * z_p_prime_base_89[2]
   + 0.5 alpha[4] * z_p_prime_base_89[2] + 0.5 p[0] * z_p_prime_base_89[2]
   - 0.5 p[1] * z_p_prime_base_89[2] - 0.5 p[2] * z_p_prime_base_89[2]
   - 0.5 p[3] * z_p_prime_base_89[2] + 0.5 p[4] * z_p_prime_base_89[2] ]
   <= 0
 infeasibility_p_prime_base_89_row_3_leq_0: 5e-06 z_p_prime_base_89[3] + [
   - alpha[1] * z_p_prime_base_89[3] + p[1] * z_p_prime_base_89[3] ] <= 0
 infeasibility_p_prime_base_89_row_4_leq_0: 5e-06 z_p_prime_base_89[4] + [
   - 0.5 alpha[0] * z_p_prime_base_89[4]
   + 0.5 alpha[1] * z_p_prime_base_89[4]
   - 0.5 alpha[2] * z_p_prime_base_89[4]
   - 0.5 alpha[3] * z_p_prime_base_89[4]
   + 0.5 alpha[4] * z_p_prime_base_89[4] - 0.5 p[0] * z_p_prime_base_89[4]
   - 0.5 p[1] * z_p_prime_base_89[4] + 0.5 p[2] * z_p_prime_base_89[4]
   + 0.5 p[3] * z_p_prime_base_89[4] + 0.5 p[4] * z_p_prime_base_89[4] ]
   <= 0
 infeasibility_p_base_90_row_0_leq_0: 5e-06 z_p_base_90[0] + [
   0.5 p[0] * z_p_base_90[0] - 0.5 p[1] * z_p_base_90[0]
   + 0.5 p[2] * z_p_base_90[0] - 0.5 p[3] * z_p_base_90[0]
   + 0.5 p[4] * z_p_base_90[0] ] <= 0
 infeasibility_p_base_90_row_1_leq_0: 5e-06 z_p_base_90[1] + [
   0.5 p[0] * z_p_base_90[1] + 0.5 p[1] * z_p_base_90[1]
   - 0.5 p[2] * z_p_base_90[1] + 0.5 p[3] * z_p_base_90[1]
   - 0.5 p[4] * z_p_base_90[1] ] <= 0
 infeasibility_p_base_90_row_2_leq_0: 5e-06 z_p_base_90[2] + [
   p[1] * z_p_base_90[2] - p[4] * z_p_base_90[2] ] <= 0
 infeasibility_p_base_90_row_3_leq_0: 5e-06 z_p_base_90[3] + [
   p[4] * z_p_base_90[3] ] <= 0
 infeasibility_p_base_90_row_4_leq_0: 5e-06 z_p_base_90[4] + [
   - 0.5 p[0] * z_p_base_90[4] - 0.5 p[1] * z_p_base_90[4]
   + 0.5 p[2] * z_p_base_90[4] + 0.5 p[3] * z_p_base_90[4]
   + 0.5 p[4] * z_p_base_90[4] ] <= 0
 infeasibility_p_prime_base_90_row_0_leq_0: 5e-06 z_p_prime_base_90[0] + [
   0.5 alpha[0] * z_p_prime_base_90[0]
   + 0.5 alpha[1] * z_p_prime_base_90[0]
   - 0.5 alpha[2] * z_p_prime_base_90[0]
   + 0.5 alpha[3] * z_p_prime_base_90[0]
   + 0.5 alpha[4] * z_p_prime_base_90[0] + 0.5 p[0] * z_p_prime_base_90[0]
   - 0.5 p[1] * z_p_prime_base_90[0] + 0.5 p[2] * z_p_prime_base_90[0]
   - 0.5 p[3] * z_p_prime_base_90[0] + 0.5 p[4] * z_p_prime_base_90[0] ]
   <= 0
 infeasibility_p_prime_base_90_row_1_leq_0: 5e-06 z_p_prime_base_90[1] + [
   0.5 alpha[0] * z_p_prime_base_90[1]
   - 0.5 alpha[1] * z_p_prime_base_90[1]
   + 0.5 alpha[2] * z_p_prime_base_90[1]
   - 0.5 alpha[3] * z_p_prime_base_90[1]
   - 0.5 alpha[4] * z_p_prime_base_90[1] + 0.5 p[0] * z_p_prime_base_90[1]
   + 0.5 p[1] * z_p_prime_base_90[1] - 0.5 p[2] * z_p_prime_base_90[1]
   + 0.5 p[3] * z_p_prime_base_90[1] - 0.5 p[4] * z_p_prime_base_90[1] ]
   <= 0
 infeasibility_p_prime_base_90_row_2_leq_0: 5e-06 z_p_prime_base_90[2] + [
   - alpha[1] * z_p_prime_base_90[2] - alpha[4] * z_p_prime_base_90[2]
   + p[1] * z_p_prime_base_90[2] - p[4] * z_p_prime_base_90[2] ] <= 0
 infeasibility_p_prime_base_90_row_3_leq_0: 5e-06 z_p_prime_base_90[3] + [
   alpha[4] * z_p_prime_base_90[3] + p[4] * z_p_prime_base_90[3] ] <= 0
 infeasibility_p_prime_base_90_row_4_leq_0: 5e-06 z_p_prime_base_90[4] + [
   - 0.5 alpha[0] * z_p_prime_base_90[4]
   + 0.5 alpha[1] * z_p_prime_base_90[4]
   - 0.5 alpha[2] * z_p_prime_base_90[4]
   - 0.5 alpha[3] * z_p_prime_base_90[4]
   + 0.5 alpha[4] * z_p_prime_base_90[4] - 0.5 p[0] * z_p_prime_base_90[4]
   - 0.5 p[1] * z_p_prime_base_90[4] + 0.5 p[2] * z_p_prime_base_90[4]
   + 0.5 p[3] * z_p_prime_base_90[4] + 0.5 p[4] * z_p_prime_base_90[4] ]
   <= 0
 infeasibility_p_base_91_row_0_leq_0: 5e-06 z_p_base_91[0] + [
   p[0] * z_p_base_91[0] - p[3] * z_p_base_91[0] ] <= 0
 infeasibility_p_base_91_row_1_leq_0: 5e-06 z_p_base_91[1] + [
   p[3] * z_p_base_91[1] ] <= 0
 infeasibility_p_base_91_row_2_leq_0: 5e-06 z_p_base_91[2] + [
   - 0.5 p[0] * z_p_base_91[2] + 0.5 p[1] * z_p_base_91[2]
   + 0.5 p[2] * z_p_base_91[2] + 0.5 p[3] * z_p_base_91[2]
   - 0.5 p[4] * z_p_base_91[2] ] <= 0
 infeasibility_p_base_91_row_3_leq_0: 5e-06 z_p_base_91[3] + [
   0.5 p[0] * z_p_base_91[3] + 0.5 p[1] * z_p_base_91[3]
   - 0.5 p[2] * z_p_base_91[3] - 0.5 p[3] * z_p_base_91[3]
   + 0.5 p[4] * z_p_base_91[3] ] <= 0
 infeasibility_p_base_91_row_4_leq_0: 5e-06 z_p_base_91[4] + [
   - 0.5 p[0] * z_p_base_91[4] - 0.5 p[1] * z_p_base_91[4]
   + 0.5 p[2] * z_p_base_91[4] + 0.5 p[3] * z_p_base_91[4]
   + 0.5 p[4] * z_p_base_91[4] ] <= 0
 infeasibility_p_prime_base_91_row_0_leq_0: 5e-06 z_p_prime_base_91[0] + [
   alpha[0] * z_p_prime_base_91[0] + alpha[3] * z_p_prime_base_91[0]
   + p[0] * z_p_prime_base_91[0] - p[3] * z_p_prime_base_91[0] ] <= 0
 infeasibility_p_prime_base_91_row_1_leq_0: 5e-06 z_p_prime_base_91[1] + [
   - alpha[3] * z_p_prime_base_91[1] + p[3] * z_p_prime_base_91[1] ] <= 0
 infeasibility_p_prime_base_91_row_2_leq_0: 5e-06 z_p_prime_base_91[2] + [
   - 0.5 alpha[0] * z_p_prime_base_91[2]
   - 0.5 alpha[1] * z_p_prime_base_91[2]
   - 0.5 alpha[2] * z_p_prime_base_91[2]
   - 0.5 alpha[3] * z_p_prime_base_91[2]
   - 0.5 alpha[4] * z_p_prime_base_91[2] - 0.5 p[0] * z_p_prime_base_91[2]
   + 0.5 p[1] * z_p_prime_base_91[2] + 0.5 p[2] * z_p_prime_base_91[2]
   + 0.5 p[3] * z_p_prime_base_91[2] - 0.5 p[4] * z_p_prime_base_91[2] ]
   <= 0
 infeasibility_p_prime_base_91_row_3_leq_0: 5e-06 z_p_prime_base_91[3] + [
   0.5 alpha[0] * z_p_prime_base_91[3]
   - 0.5 alpha[1] * z_p_prime_base_91[3]
   + 0.5 alpha[2] * z_p_prime_base_91[3]
   + 0.5 alpha[3] * z_p_prime_base_91[3]
   + 0.5 alpha[4] * z_p_prime_base_91[3] + 0.5 p[0] * z_p_prime_base_91[3]
   + 0.5 p[1] * z_p_prime_base_91[3] - 0.5 p[2] * z_p_prime_base_91[3]
   - 0.5 p[3] * z_p_prime_base_91[3] + 0.5 p[4] * z_p_prime_base_91[3] ]
   <= 0
 infeasibility_p_prime_base_91_row_4_leq_0: 5e-06 z_p_prime_base_91[4] + [
   - 0.5 alpha[0] * z_p_prime_base_91[4]
   + 0.5 alpha[1] * z_p_prime_base_91[4]
   - 0.5 alpha[2] * z_p_prime_base_91[4]
   - 0.5 alpha[3] * z_p_prime_base_91[4]
   + 0.5 alpha[4] * z_p_prime_base_91[4] - 0.5 p[0] * z_p_prime_base_91[4]
   - 0.5 p[1] * z_p_prime_base_91[4] + 0.5 p[2] * z_p_prime_base_91[4]
   + 0.5 p[3] * z_p_prime_base_91[4] + 0.5 p[4] * z_p_prime_base_91[4] ]
   <= 0
 infeasibility_p_base_92_row_0_leq_0: 5e-06 z_p_base_92[0] + [
   0.5 p[0] * z_p_base_92[0] - 0.5 p[1] * z_p_base_92[0]
   + 0.5 p[2] * z_p_base_92[0] - 0.5 p[3] * z_p_base_92[0]
   + 0.5 p[4] * z_p_base_92[0] ] <= 0
 infeasibility_p_base_92_row_1_leq_0: 5e-06 z_p_base_92[1] + [
   0.5 p[0] * z_p_base_92[1] + 0.5 p[1] * z_p_base_92[1]
   - 0.5 p[2] * z_p_base_92[1] + 0.5 p[3] * z_p_base_92[1]
   - 0.5 p[4] * z_p_base_92[1] ] <= 0
 infeasibility_p_base_92_row_2_leq_0: 5e-06 z_p_base_92[2] + [
   - 0.5 p[0] * z_p_base_92[2] + 0.5 p[1] * z_p_base_92[2]
   + 0.5 p[2] * z_p_base_92[2] + 0.5 p[3] * z_p_base_92[2]
   - 0.5 p[4] * z_p_base_92[2] ] <= 0
 infeasibility_p_base_92_row_3_leq_0: 5e-06 z_p_base_92[3] + [
   0.5 p[0] * z_p_base_92[3] + 0.5 p[1] * z_p_base_92[3]
   - 0.5 p[2] * z_p_base_92[3] - 0.5 p[3] * z_p_base_92[3]
   + 0.5 p[4] * z_p_base_92[3] ] <= 0
 infeasibility_p_base_92_row_4_leq_0: 5e-06 z_p_base_92[4] + [
   - 0.5 p[0] * z_p_base_92[4] - 0.5 p[1] * z_p_base_92[4]
   + 0.5 p[2] * z_p_base_92[4] + 0.5 p[3] * z_p_base_92[4]
   + 0.5 p[4] * z_p_base_92[4] ] <= 0
 infeasibility_p_prime_base_92_row_0_leq_0: 5e-06 z_p_prime_base_92[0] + [
   0.5 alpha[0] * z_p_prime_base_92[0]
   + 0.5 alpha[1] * z_p_prime_base_92[0]
   - 0.5 alpha[2] * z_p_prime_base_92[0]
   + 0.5 alpha[3] * z_p_prime_base_92[0]
   + 0.5 alpha[4] * z_p_prime_base_92[0] + 0.5 p[0] * z_p_prime_base_92[0]
   - 0.5 p[1] * z_p_prime_base_92[0] + 0.5 p[2] * z_p_prime_base_92[0]
   - 0.5 p[3] * z_p_prime_base_92[0] + 0.5 p[4] * z_p_prime_base_92[0] ]
   <= 0
 infeasibility_p_prime_base_92_row_1_leq_0: 5e-06 z_p_prime_base_92[1] + [
   0.5 alpha[0] * z_p_prime_base_92[1]
   - 0.5 alpha[1] * z_p_prime_base_92[1]
   + 0.5 alpha[2] * z_p_prime_base_92[1]
   - 0.5 alpha[3] * z_p_prime_base_92[1]
   - 0.5 alpha[4] * z_p_prime_base_92[1] + 0.5 p[0] * z_p_prime_base_92[1]
   + 0.5 p[1] * z_p_prime_base_92[1] - 0.5 p[2] * z_p_prime_base_92[1]
   + 0.5 p[3] * z_p_prime_base_92[1] - 0.5 p[4] * z_p_prime_base_92[1] ]
   <= 0
 infeasibility_p_prime_base_92_row_2_leq_0: 5e-06 z_p_prime_base_92[2] + [
   - 0.5 alpha[0] * z_p_prime_base_92[2]
   - 0.5 alpha[1] * z_p_prime_base_92[2]
   - 0.5 alpha[2] * z_p_prime_base_92[2]
   - 0.5 alpha[3] * z_p_prime_base_92[2]
   - 0.5 alpha[4] * z_p_prime_base_92[2] - 0.5 p[0] * z_p_prime_base_92[2]
   + 0.5 p[1] * z_p_prime_base_92[2] + 0.5 p[2] * z_p_prime_base_92[2]
   + 0.5 p[3] * z_p_prime_base_92[2] - 0.5 p[4] * z_p_prime_base_92[2] ]
   <= 0
 infeasibility_p_prime_base_92_row_3_leq_0: 5e-06 z_p_prime_base_92[3] + [
   0.5 alpha[0] * z_p_prime_base_92[3]
   - 0.5 alpha[1] * z_p_prime_base_92[3]
   + 0.5 alpha[2] * z_p_prime_base_92[3]
   + 0.5 alpha[3] * z_p_prime_base_92[3]
   + 0.5 alpha[4] * z_p_prime_base_92[3] + 0.5 p[0] * z_p_prime_base_92[3]
   + 0.5 p[1] * z_p_prime_base_92[3] - 0.5 p[2] * z_p_prime_base_92[3]
   - 0.5 p[3] * z_p_prime_base_92[3] + 0.5 p[4] * z_p_prime_base_92[3] ]
   <= 0
 infeasibility_p_prime_base_92_row_4_leq_0: 5e-06 z_p_prime_base_92[4] + [
   - 0.5 alpha[0] * z_p_prime_base_92[4]
   + 0.5 alpha[1] * z_p_prime_base_92[4]
   - 0.5 alpha[2] * z_p_prime_base_92[4]
   - 0.5 alpha[3] * z_p_prime_base_92[4]
   + 0.5 alpha[4] * z_p_prime_base_92[4] - 0.5 p[0] * z_p_prime_base_92[4]
   - 0.5 p[1] * z_p_prime_base_92[4] + 0.5 p[2] * z_p_prime_base_92[4]
   + 0.5 p[3] * z_p_prime_base_92[4] + 0.5 p[4] * z_p_prime_base_92[4] ]
   <= 0
 infeasibility_p_base_93_row_0_leq_0: 5e-06 z_p_base_93[0] + [
   0.5 p[0] * z_p_base_93[0] - 0.5 p[1] * z_p_base_93[0]
   + 0.5 p[2] * z_p_base_93[0] - 0.5 p[3] * z_p_base_93[0]
   - 0.5 p[4] * z_p_base_93[0] ] <= 0
 infeasibility_p_base_93_row_1_leq_0: 5e-06 z_p_base_93[1] + [
   0.5 p[0] * z_p_base_93[1] + 0.5 p[1] * z_p_base_93[1]
   - 0.5 p[2] * z_p_base_93[1] + 0.5 p[3] * z_p_base_93[1]
   + 0.5 p[4] * z_p_base_93[1] ] <= 0
 infeasibility_p_base_93_row_2_leq_0: 5e-06 z_p_base_93[2] + [
   p[1] * z_p_base_93[2] ] <= 0
 infeasibility_p_base_93_row_3_leq_0: 5e-06 z_p_base_93[3] + [
   - 0.5 p[0] * z_p_base_93[3] - 0.5 p[1] * z_p_base_93[3]
   + 0.5 p[2] * z_p_base_93[3] + 0.5 p[3] * z_p_base_93[3]
   - 0.5 p[4] * z_p_base_93[3] ] <= 0
 infeasibility_p_base_93_row_4_leq_0: 5e-06 z_p_base_93[4] + [
   p[4] * z_p_base_93[4] ] <= 0
 infeasibility_p_prime_base_93_row_0_leq_0: 5e-06 z_p_prime_base_93[0] + [
   0.5 alpha[0] * z_p_prime_base_93[0]
   + 0.5 alpha[1] * z_p_prime_base_93[0]
   - 0.5 alpha[2] * z_p_prime_base_93[0]
   + 0.5 alpha[3] * z_p_prime_base_93[0]
   - 0.5 alpha[4] * z_p_prime_base_93[0] + 0.5 p[0] * z_p_prime_base_93[0]
   - 0.5 p[1] * z_p_prime_base_93[0] + 0.5 p[2] * z_p_prime_base_93[0]
   - 0.5 p[3] * z_p_prime_base_93[0] - 0.5 p[4] * z_p_prime_base_93[0] ]
   <= 0
 infeasibility_p_prime_base_93_row_1_leq_0: 5e-06 z_p_prime_base_93[1] + [
   0.5 alpha[0] * z_p_prime_base_93[1]
   - 0.5 alpha[1] * z_p_prime_base_93[1]
   + 0.5 alpha[2] * z_p_prime_base_93[1]
   - 0.5 alpha[3] * z_p_prime_base_93[1]
   + 0.5 alpha[4] * z_p_prime_base_93[1] + 0.5 p[0] * z_p_prime_base_93[1]
   + 0.5 p[1] * z_p_prime_base_93[1] - 0.5 p[2] * z_p_prime_base_93[1]
   + 0.5 p[3] * z_p_prime_base_93[1] + 0.5 p[4] * z_p_prime_base_93[1] ]
   <= 0
 infeasibility_p_prime_base_93_row_2_leq_0: 5e-06 z_p_prime_base_93[2] + [
   - alpha[1] * z_p_prime_base_93[2] + p[1] * z_p_prime_base_93[2] ] <= 0
 infeasibility_p_prime_base_93_row_3_leq_0: 5e-06 z_p_prime_base_93[3] + [
   - 0.5 alpha[0] * z_p_prime_base_93[3]
   + 0.5 alpha[1] * z_p_prime_base_93[3]
   - 0.5 alpha[2] * z_p_prime_base_93[3]
   - 0.5 alpha[3] * z_p_prime_base_93[3]
   - 0.5 alpha[4] * z_p_prime_base_93[3] - 0.5 p[0] * z_p_prime_base_93[3]
   - 0.5 p[1] * z_p_prime_base_93[3] + 0.5 p[2] * z_p_prime_base_93[3]
   + 0.5 p[3] * z_p_prime_base_93[3] - 0.5 p[4] * z_p_prime_base_93[3] ]
   <= 0
 infeasibility_p_prime_base_93_row_4_leq_0: 5e-06 z_p_prime_base_93[4] + [
   alpha[4] * z_p_prime_base_93[4] + p[4] * z_p_prime_base_93[4] ] <= 0
 infeasibility_p_base_94_row_0_leq_0: 5e-06 z_p_base_94[0] + [
   0.5 p[0] * z_p_base_94[0] - 0.5 p[1] * z_p_base_94[0]
   + 0.5 p[2] * z_p_base_94[0] - 0.5 p[3] * z_p_base_94[0]
   + 0.5 p[4] * z_p_base_94[0] ] <= 0
 infeasibility_p_base_94_row_1_leq_0: 5e-06 z_p_base_94[1] + [
   0.5 p[0] * z_p_base_94[1] + 0.5 p[1] * z_p_base_94[1]
   - 0.5 p[2] * z_p_base_94[1] + 0.5 p[3] * z_p_base_94[1]
   - 0.5 p[4] * z_p_base_94[1] ] <= 0
 infeasibility_p_base_94_row_2_leq_0: 5e-06 z_p_base_94[2] + [
   p[1] * z_p_base_94[2] ] <= 0
 infeasibility_p_base_94_row_3_leq_0: 5e-06 z_p_base_94[3] + [
   - 0.5 p[0] * z_p_base_94[3] - 0.5 p[1] * z_p_base_94[3]
   + 0.5 p[2] * z_p_base_94[3] + 0.5 p[3] * z_p_base_94[3]
   - 0.5 p[4] * z_p_base_94[3] ] <= 0
 infeasibility_p_base_94_row_4_leq_0: 5e-06 z_p_base_94[4] + [
   p[4] * z_p_base_94[4] ] <= 0
 infeasibility_p_prime_base_94_row_0_leq_0: 5e-06 z_p_prime_base_94[0] + [
   0.5 alpha[0] * z_p_prime_base_94[0]
   + 0.5 alpha[1] * z_p_prime_base_94[0]
   - 0.5 alpha[2] * z_p_prime_base_94[0]
   + 0.5 alpha[3] * z_p_prime_base_94[0]
   + 0.5 alpha[4] * z_p_prime_base_94[0] + 0.5 p[0] * z_p_prime_base_94[0]
   - 0.5 p[1] * z_p_prime_base_94[0] + 0.5 p[2] * z_p_prime_base_94[0]
   - 0.5 p[3] * z_p_prime_base_94[0] + 0.5 p[4] * z_p_prime_base_94[0] ]
   <= 0
 infeasibility_p_prime_base_94_row_1_leq_0: 5e-06 z_p_prime_base_94[1] + [
   0.5 alpha[0] * z_p_prime_base_94[1]
   - 0.5 alpha[1] * z_p_prime_base_94[1]
   + 0.5 alpha[2] * z_p_prime_base_94[1]
   - 0.5 alpha[3] * z_p_prime_base_94[1]
   - 0.5 alpha[4] * z_p_prime_base_94[1] + 0.5 p[0] * z_p_prime_base_94[1]
   + 0.5 p[1] * z_p_prime_base_94[1] - 0.5 p[2] * z_p_prime_base_94[1]
   + 0.5 p[3] * z_p_prime_base_94[1] - 0.5 p[4] * z_p_prime_base_94[1] ]
   <= 0
 infeasibility_p_prime_base_94_row_2_leq_0: 5e-06 z_p_prime_base_94[2] + [
   - alpha[1] * z_p_prime_base_94[2] + p[1] * z_p_prime_base_94[2] ] <= 0
 infeasibility_p_prime_base_94_row_3_leq_0: 5e-06 z_p_prime_base_94[3] + [
   - 0.5 alpha[0] * z_p_prime_base_94[3]
   + 0.5 alpha[1] * z_p_prime_base_94[3]
   - 0.5 alpha[2] * z_p_prime_base_94[3]
   - 0.5 alpha[3] * z_p_prime_base_94[3]
   - 0.5 alpha[4] * z_p_prime_base_94[3] - 0.5 p[0] * z_p_prime_base_94[3]
   - 0.5 p[1] * z_p_prime_base_94[3] + 0.5 p[2] * z_p_prime_base_94[3]
   + 0.5 p[3] * z_p_prime_base_94[3] - 0.5 p[4] * z_p_prime_base_94[3] ]
   <= 0
 infeasibility_p_prime_base_94_row_4_leq_0: 5e-06 z_p_prime_base_94[4] + [
   alpha[4] * z_p_prime_base_94[4] + p[4] * z_p_prime_base_94[4] ] <= 0
 infeasibility_p_base_95_row_0_leq_0: 5e-06 z_p_base_95[0] + [
   0.5 p[0] * z_p_base_95[0] + 0.5 p[1] * z_p_base_95[0]
   + 0.5 p[2] * z_p_base_95[0] - 0.5 p[3] * z_p_base_95[0]
   - 0.5 p[4] * z_p_base_95[0] ] <= 0
 infeasibility_p_base_95_row_1_leq_0: 5e-06 z_p_base_95[1] + [
   0.5 p[0] * z_p_base_95[1] - 0.5 p[1] * z_p_base_95[1]
   - 0.5 p[2] * z_p_base_95[1] + 0.5 p[3] * z_p_base_95[1]
   + 0.5 p[4] * z_p_base_95[1] ] <= 0
 infeasibility_p_base_95_row_2_leq_0: 5e-06 z_p_base_95[2] + [
   p[1] * z_p_base_95[2] - p[4] * z_p_base_95[2] ] <= 0
 infeasibility_p_base_95_row_3_leq_0: 5e-06 z_p_base_95[3] + [
   p[4] * z_p_base_95[3] ] <= 0
 infeasibility_p_base_95_row_4_leq_0: 5e-06 z_p_base_95[4] + [
   - 0.5 p[0] * z_p_base_95[4] - 0.5 p[1] * z_p_base_95[4]
   + 0.5 p[2] * z_p_base_95[4] + 0.5 p[3] * z_p_base_95[4]
   + 0.5 p[4] * z_p_base_95[4] ] <= 0
 infeasibility_p_prime_base_95_row_0_leq_0: 5e-06 z_p_prime_base_95[0] + [
   0.5 alpha[0] * z_p_prime_base_95[0]
   - 0.5 alpha[1] * z_p_prime_base_95[0]
   - 0.5 alpha[2] * z_p_prime_base_95[0]
   + 0.5 alpha[3] * z_p_prime_base_95[0]
   - 0.5 alpha[4] * z_p_prime_base_95[0] + 0.5 p[0] * z_p_prime_base_95[0]
   + 0.5 p[1] * z_p_prime_base_95[0] + 0.5 p[2] * z_p_prime_base_95[0]
   - 0.5 p[3] * z_p_prime_base_95[0] - 0.5 p[4] * z_p_prime_base_95[0] ]
   <= 0
 infeasibility_p_prime_base_95_row_1_leq_0: 5e-06 z_p_prime_base_95[1] + [
   0.5 alpha[0] * z_p_prime_base_95[1]
   + 0.5 alpha[1] * z_p_prime_base_95[1]
   + 0.5 alpha[2] * z_p_prime_base_95[1]
   - 0.5 alpha[3] * z_p_prime_base_95[1]
   + 0.5 alpha[4] * z_p_prime_base_95[1] + 0.5 p[0] * z_p_prime_base_95[1]
   - 0.5 p[1] * z_p_prime_base_95[1] - 0.5 p[2] * z_p_prime_base_95[1]
   + 0.5 p[3] * z_p_prime_base_95[1] + 0.5 p[4] * z_p_prime_base_95[1] ]
   <= 0
 infeasibility_p_prime_base_95_row_2_leq_0: 5e-06 z_p_prime_base_95[2] + [
   - alpha[1] * z_p_prime_base_95[2] - alpha[4] * z_p_prime_base_95[2]
   + p[1] * z_p_prime_base_95[2] - p[4] * z_p_prime_base_95[2] ] <= 0
 infeasibility_p_prime_base_95_row_3_leq_0: 5e-06 z_p_prime_base_95[3] + [
   alpha[4] * z_p_prime_base_95[3] + p[4] * z_p_prime_base_95[3] ] <= 0
 infeasibility_p_prime_base_95_row_4_leq_0: 5e-06 z_p_prime_base_95[4] + [
   - 0.5 alpha[0] * z_p_prime_base_95[4]
   + 0.5 alpha[1] * z_p_prime_base_95[4]
   - 0.5 alpha[2] * z_p_prime_base_95[4]
   - 0.5 alpha[3] * z_p_prime_base_95[4]
   + 0.5 alpha[4] * z_p_prime_base_95[4] - 0.5 p[0] * z_p_prime_base_95[4]
   - 0.5 p[1] * z_p_prime_base_95[4] + 0.5 p[2] * z_p_prime_base_95[4]
   + 0.5 p[3] * z_p_prime_base_95[4] + 0.5 p[4] * z_p_prime_base_95[4] ]
   <= 0
 infeasibility_p_base_96_row_0_leq_0: 5e-06 z_p_base_96[0] + [
   0.5 p[0] * z_p_base_96[0] + 0.5 p[1] * z_p_base_96[0]
   + 0.5 p[2] * z_p_base_96[0] - 0.5 p[3] * z_p_base_96[0]
   - 0.5 p[4] * z_p_base_96[0] ] <= 0
 infeasibility_p_base_96_row_1_leq_0: 5e-06 z_p_base_96[1] + [
   0.5 p[0] * z_p_base_96[1] - 0.5 p[1] * z_p_base_96[1]
   - 0.5 p[2] * z_p_base_96[1] + 0.5 p[3] * z_p_base_96[1]
   + 0.5 p[4] * z_p_base_96[1] ] <= 0
 infeasibility_p_base_96_row_2_leq_0: 5e-06 z_p_base_96[2] + [
   - 0.5 p[0] * z_p_base_96[2] + 0.5 p[1] * z_p_base_96[2]
   + 0.5 p[2] * z_p_base_96[2] + 0.5 p[3] * z_p_base_96[2]
   - 0.5 p[4] * z_p_base_96[2] ] <= 0
 infeasibility_p_base_96_row_3_leq_0: 5e-06 z_p_base_96[3] + [
   0.5 p[0] * z_p_base_96[3] + 0.5 p[1] * z_p_base_96[3]
   - 0.5 p[2] * z_p_base_96[3] - 0.5 p[3] * z_p_base_96[3]
   + 0.5 p[4] * z_p_base_96[3] ] <= 0
 infeasibility_p_base_96_row_4_leq_0: 5e-06 z_p_base_96[4] + [
   - 0.5 p[0] * z_p_base_96[4] - 0.5 p[1] * z_p_base_96[4]
   + 0.5 p[2] * z_p_base_96[4] + 0.5 p[3] * z_p_base_96[4]
   + 0.5 p[4] * z_p_base_96[4] ] <= 0
 infeasibility_p_prime_base_96_row_0_leq_0: 5e-06 z_p_prime_base_96[0] + [
   0.5 alpha[0] * z_p_prime_base_96[0]
   - 0.5 alpha[1] * z_p_prime_base_96[0]
   - 0.5 alpha[2] * z_p_prime_base_96[0]
   + 0.5 alpha[3] * z_p_prime_base_96[0]
   - 0.5 alpha[4] * z_p_prime_base_96[0] + 0.5 p[0] * z_p_prime_base_96[0]
   + 0.5 p[1] * z_p_prime_base_96[0] + 0.5 p[2] * z_p_prime_base_96[0]
   - 0.5 p[3] * z_p_prime_base_96[0] - 0.5 p[4] * z_p_prime_base_96[0] ]
   <= 0
 infeasibility_p_prime_base_96_row_1_leq_0: 5e-06 z_p_prime_base_96[1] + [
   0.5 alpha[0] * z_p_prime_base_96[1]
   + 0.5 alpha[1] * z_p_prime_base_96[1]
   + 0.5 alpha[2] * z_p_prime_base_96[1]
   - 0.5 alpha[3] * z_p_prime_base_96[1]
   + 0.5 alpha[4] * z_p_prime_base_96[1] + 0.5 p[0] * z_p_prime_base_96[1]
   - 0.5 p[1] * z_p_prime_base_96[1] - 0.5 p[2] * z_p_prime_base_96[1]
   + 0.5 p[3] * z_p_prime_base_96[1] + 0.5 p[4] * z_p_prime_base_96[1] ]
   <= 0
 infeasibility_p_prime_base_96_row_2_leq_0: 5e-06 z_p_prime_base_96[2] + [
   - 0.5 alpha[0] * z_p_prime_base_96[2]
   - 0.5 alpha[1] * z_p_prime_base_96[2]
   - 0.5 alpha[2] * z_p_prime_base_96[2]
   - 0.5 alpha[3] * z_p_prime_base_96[2]
   - 0.5 alpha[4] * z_p_prime_base_96[2] - 0.5 p[0] * z_p_prime_base_96[2]
   + 0.5 p[1] * z_p_prime_base_96[2] + 0.5 p[2] * z_p_prime_base_96[2]
   + 0.5 p[3] * z_p_prime_base_96[2] - 0.5 p[4] * z_p_prime_base_96[2] ]
   <= 0
 infeasibility_p_prime_base_96_row_3_leq_0: 5e-06 z_p_prime_base_96[3] + [
   0.5 alpha[0] * z_p_prime_base_96[3]
   - 0.5 alpha[1] * z_p_prime_base_96[3]
   + 0.5 alpha[2] * z_p_prime_base_96[3]
   + 0.5 alpha[3] * z_p_prime_base_96[3]
   + 0.5 alpha[4] * z_p_prime_base_96[3] + 0.5 p[0] * z_p_prime_base_96[3]
   + 0.5 p[1] * z_p_prime_base_96[3] - 0.5 p[2] * z_p_prime_base_96[3]
   - 0.5 p[3] * z_p_prime_base_96[3] + 0.5 p[4] * z_p_prime_base_96[3] ]
   <= 0
 infeasibility_p_prime_base_96_row_4_leq_0: 5e-06 z_p_prime_base_96[4] + [
   - 0.5 alpha[0] * z_p_prime_base_96[4]
   + 0.5 alpha[1] * z_p_prime_base_96[4]
   - 0.5 alpha[2] * z_p_prime_base_96[4]
   - 0.5 alpha[3] * z_p_prime_base_96[4]
   + 0.5 alpha[4] * z_p_prime_base_96[4] - 0.5 p[0] * z_p_prime_base_96[4]
   - 0.5 p[1] * z_p_prime_base_96[4] + 0.5 p[2] * z_p_prime_base_96[4]
   + 0.5 p[3] * z_p_prime_base_96[4] + 0.5 p[4] * z_p_prime_base_96[4] ]
   <= 0
 infeasibility_p_base_97_row_0_leq_0: 5e-06 z_p_base_97[0] + [
   p[2] * z_p_base_97[0] ] <= 0
 infeasibility_p_base_97_row_1_leq_0: 5e-06 z_p_base_97[1] + [
   p[0] * z_p_base_97[1] - p[2] * z_p_base_97[1] ] <= 0
 infeasibility_p_base_97_row_2_leq_0: 5e-06 z_p_base_97[2] + [
   - 0.5 p[0] * z_p_base_97[2] + 0.5 p[1] * z_p_base_97[2]
   + 0.5 p[2] * z_p_base_97[2] + 0.5 p[3] * z_p_base_97[2]
   - 0.5 p[4] * z_p_base_97[2] ] <= 0
 infeasibility_p_base_97_row_3_leq_0: 5e-06 z_p_base_97[3] + [
   0.5 p[0] * z_p_base_97[3] + 0.5 p[1] * z_p_base_97[3]
   - 0.5 p[2] * z_p_base_97[3] - 0.5 p[3] * z_p_base_97[3]
   + 0.5 p[4] * z_p_base_97[3] ] <= 0
 infeasibility_p_base_97_row_4_leq_0: 5e-06 z_p_base_97[4] + [
   - 0.5 p[0] * z_p_base_97[4] - 0.5 p[1] * z_p_base_97[4]
   + 0.5 p[2] * z_p_base_97[4] + 0.5 p[3] * z_p_base_97[4]
   + 0.5 p[4] * z_p_base_97[4] ] <= 0
 infeasibility_p_prime_base_97_row_0_leq_0: 5e-06 z_p_prime_base_97[0] + [
   - alpha[2] * z_p_prime_base_97[0] + p[2] * z_p_prime_base_97[0] ] <= 0
 infeasibility_p_prime_base_97_row_1_leq_0: 5e-06 z_p_prime_base_97[1] + [
   alpha[0] * z_p_prime_base_97[1] + alpha[2] * z_p_prime_base_97[1]
   + p[0] * z_p_prime_base_97[1] - p[2] * z_p_prime_base_97[1] ] <= 0
 infeasibility_p_prime_base_97_row_2_leq_0: 5e-06 z_p_prime_base_97[2] + [
   - 0.5 alpha[0] * z_p_prime_base_97[2]
   - 0.5 alpha[1] * z_p_prime_base_97[2]
   - 0.5 alpha[2] * z_p_prime_base_97[2]
   - 0.5 alpha[3] * z_p_prime_base_97[2]
   - 0.5 alpha[4] * z_p_prime_base_97[2] - 0.5 p[0] * z_p_prime_base_97[2]
   + 0.5 p[1] * z_p_prime_base_97[2] + 0.5 p[2] * z_p_prime_base_97[2]
   + 0.5 p[3] * z_p_prime_base_97[2] - 0.5 p[4] * z_p_prime_base_97[2] ]
   <= 0
 infeasibility_p_prime_base_97_row_3_leq_0: 5e-06 z_p_prime_base_97[3] + [
   0.5 alpha[0] * z_p_prime_base_97[3]
   - 0.5 alpha[1] * z_p_prime_base_97[3]
   + 0.5 alpha[2] * z_p_prime_base_97[3]
   + 0.5 alpha[3] * z_p_prime_base_97[3]
   + 0.5 alpha[4] * z_p_prime_base_97[3] + 0.5 p[0] * z_p_prime_base_97[3]
   + 0.5 p[1] * z_p_prime_base_97[3] - 0.5 p[2] * z_p_prime_base_97[3]
   - 0.5 p[3] * z_p_prime_base_97[3] + 0.5 p[4] * z_p_prime_base_97[3] ]
   <= 0
 infeasibility_p_prime_base_97_row_4_leq_0: 5e-06 z_p_prime_base_97[4] + [
   - 0.5 alpha[0] * z_p_prime_base_97[4]
   + 0.5 alpha[1] * z_p_prime_base_97[4]
   - 0.5 alpha[2] * z_p_prime_base_97[4]
   - 0.5 alpha[3] * z_p_prime_base_97[4]
   + 0.5 alpha[4] * z_p_prime_base_97[4] - 0.5 p[0] * z_p_prime_base_97[4]
   - 0.5 p[1] * z_p_prime_base_97[4] + 0.5 p[2] * z_p_prime_base_97[4]
   + 0.5 p[3] * z_p_prime_base_97[4] + 0.5 p[4] * z_p_prime_base_97[4] ]
   <= 0
 infeasibility_p_base_98_row_0_leq_0: 5e-06 z_p_base_98[0] + [
   0.5 p[0] * z_p_base_98[0] + 0.5 p[1] * z_p_base_98[0]
   + 0.5 p[2] * z_p_base_98[0] - 0.5 p[3] * z_p_base_98[0]
   - 0.5 p[4] * z_p_base_98[0] ] <= 0
 infeasibility_p_base_98_row_1_leq_0: 5e-06 z_p_base_98[1] + [
   0.5 p[0] * z_p_base_98[1] - 0.5 p[1] * z_p_base_98[1]
   - 0.5 p[2] * z_p_base_98[1] + 0.5 p[3] * z_p_base_98[1]
   + 0.5 p[4] * z_p_base_98[1] ] <= 0
 infeasibility_p_base_98_row_2_leq_0: 5e-06 z_p_base_98[2] + [
   p[1] * z_p_base_98[2] ] <= 0
 infeasibility_p_base_98_row_3_leq_0: 5e-06 z_p_base_98[3] + [
   - 0.5 p[0] * z_p_base_98[3] - 0.5 p[1] * z_p_base_98[3]
   + 0.5 p[2] * z_p_base_98[3] + 0.5 p[3] * z_p_base_98[3]
   - 0.5 p[4] * z_p_base_98[3] ] <= 0
 infeasibility_p_base_98_row_4_leq_0: 5e-06 z_p_base_98[4] + [
   p[4] * z_p_base_98[4] ] <= 0
 infeasibility_p_prime_base_98_row_0_leq_0: 5e-06 z_p_prime_base_98[0] + [
   0.5 alpha[0] * z_p_prime_base_98[0]
   - 0.5 alpha[1] * z_p_prime_base_98[0]
   - 0.5 alpha[2] * z_p_prime_base_98[0]
   + 0.5 alpha[3] * z_p_prime_base_98[0]
   - 0.5 alpha[4] * z_p_prime_base_98[0] + 0.5 p[0] * z_p_prime_base_98[0]
   + 0.5 p[1] * z_p_prime_base_98[0] + 0.5 p[2] * z_p_prime_base_98[0]
   - 0.5 p[3] * z_p_prime_base_98[0] - 0.5 p[4] * z_p_prime_base_98[0] ]
   <= 0
 infeasibility_p_prime_base_98_row_1_leq_0: 5e-06 z_p_prime_base_98[1] + [
   0.5 alpha[0] * z_p_prime_base_98[1]
   + 0.5 alpha[1] * z_p_prime_base_98[1]
   + 0.5 alpha[2] * z_p_prime_base_98[1]
   - 0.5 alpha[3] * z_p_prime_base_98[1]
   + 0.5 alpha[4] * z_p_prime_base_98[1] + 0.5 p[0] * z_p_prime_base_98[1]
   - 0.5 p[1] * z_p_prime_base_98[1] - 0.5 p[2] * z_p_prime_base_98[1]
   + 0.5 p[3] * z_p_prime_base_98[1] + 0.5 p[4] * z_p_prime_base_98[1] ]
   <= 0
 infeasibility_p_prime_base_98_row_2_leq_0: 5e-06 z_p_prime_base_98[2] + [
   - alpha[1] * z_p_prime_base_98[2] + p[1] * z_p_prime_base_98[2] ] <= 0
 infeasibility_p_prime_base_98_row_3_leq_0: 5e-06 z_p_prime_base_98[3] + [
   - 0.5 alpha[0] * z_p_prime_base_98[3]
   + 0.5 alpha[1] * z_p_prime_base_98[3]
   - 0.5 alpha[2] * z_p_prime_base_98[3]
   - 0.5 alpha[3] * z_p_prime_base_98[3]
   - 0.5 alpha[4] * z_p_prime_base_98[3] - 0.5 p[0] * z_p_prime_base_98[3]
   - 0.5 p[1] * z_p_prime_base_98[3] + 0.5 p[2] * z_p_prime_base_98[3]
   + 0.5 p[3] * z_p_prime_base_98[3] - 0.5 p[4] * z_p_prime_base_98[3] ]
   <= 0
 infeasibility_p_prime_base_98_row_4_leq_0: 5e-06 z_p_prime_base_98[4] + [
   alpha[4] * z_p_prime_base_98[4] + p[4] * z_p_prime_base_98[4] ] <= 0
 infeasibility_p_base_99_row_0_leq_0: 5e-06 z_p_base_99[0] + [
   0.5 p[0] * z_p_base_99[0] + 0.5 p[1] * z_p_base_99[0]
   + 0.5 p[2] * z_p_base_99[0] - 0.5 p[3] * z_p_base_99[0]
   + 0.5 p[4] * z_p_base_99[0] ] <= 0
 infeasibility_p_base_99_row_1_leq_0: 5e-06 z_p_base_99[1] + [
   0.5 p[0] * z_p_base_99[1] - 0.5 p[1] * z_p_base_99[1]
   - 0.5 p[2] * z_p_base_99[1] + 0.5 p[3] * z_p_base_99[1]
   - 0.5 p[4] * z_p_base_99[1] ] <= 0
 infeasibility_p_base_99_row_2_leq_0: 5e-06 z_p_base_99[2] + [
   p[1] * z_p_base_99[2] ] <= 0
 infeasibility_p_base_99_row_3_leq_0: 5e-06 z_p_base_99[3] + [
   - 0.5 p[0] * z_p_base_99[3] - 0.5 p[1] * z_p_base_99[3]
   + 0.5 p[2] * z_p_base_99[3] + 0.5 p[3] * z_p_base_99[3]
   - 0.5 p[4] * z_p_base_99[3] ] <= 0
 infeasibility_p_base_99_row_4_leq_0: 5e-06 z_p_base_99[4] + [
   p[4] * z_p_base_99[4] ] <= 0
 infeasibility_p_prime_base_99_row_0_leq_0: 5e-06 z_p_prime_base_99[0] + [
   0.5 alpha[0] * z_p_prime_base_99[0]
   - 0.5 alpha[1] * z_p_prime_base_99[0]
   - 0.5 alpha[2] * z_p_prime_base_99[0]
   + 0.5 alpha[3] * z_p_prime_base_99[0]
   + 0.5 alpha[4] * z_p_prime_base_99[0] + 0.5 p[0] * z_p_prime_base_99[0]
   + 0.5 p[1] * z_p_prime_base_99[0] + 0.5 p[2] * z_p_prime_base_99[0]
   - 0.5 p[3] * z_p_prime_base_99[0] + 0.5 p[4] * z_p_prime_base_99[0] ]
   <= 0
 infeasibility_p_prime_base_99_row_1_leq_0: 5e-06 z_p_prime_base_99[1] + [
   0.5 alpha[0] * z_p_prime_base_99[1]
   + 0.5 alpha[1] * z_p_prime_base_99[1]
   + 0.5 alpha[2] * z_p_prime_base_99[1]
   - 0.5 alpha[3] * z_p_prime_base_99[1]
   - 0.5 alpha[4] * z_p_prime_base_99[1] + 0.5 p[0] * z_p_prime_base_99[1]
   - 0.5 p[1] * z_p_prime_base_99[1] - 0.5 p[2] * z_p_prime_base_99[1]
   + 0.5 p[3] * z_p_prime_base_99[1] - 0.5 p[4] * z_p_prime_base_99[1] ]
   <= 0
 infeasibility_p_prime_base_99_row_2_leq_0: 5e-06 z_p_prime_base_99[2] + [
   - alpha[1] * z_p_prime_base_99[2] + p[1] * z_p_prime_base_99[2] ] <= 0
 infeasibility_p_prime_base_99_row_3_leq_0: 5e-06 z_p_prime_base_99[3] + [
   - 0.5 alpha[0] * z_p_prime_base_99[3]
   + 0.5 alpha[1] * z_p_prime_base_99[3]
   - 0.5 alpha[2] * z_p_prime_base_99[3]
   - 0.5 alpha[3] * z_p_prime_base_99[3]
   - 0.5 alpha[4] * z_p_prime_base_99[3] - 0.5 p[0] * z_p_prime_base_99[3]
   - 0.5 p[1] * z_p_prime_base_99[3] + 0.5 p[2] * z_p_prime_base_99[3]
   + 0.5 p[3] * z_p_prime_base_99[3] - 0.5 p[4] * z_p_prime_base_99[3] ]
   <= 0
 infeasibility_p_prime_base_99_row_4_leq_0: 5e-06 z_p_prime_base_99[4] + [
   alpha[4] * z_p_prime_base_99[4] + p[4] * z_p_prime_base_99[4] ] <= 0
 infeasibility_p_base_100_row_0_leq_0: 5e-06 z_p_base_100[0] + [
   0.5 p[0] * z_p_base_100[0] + 0.5 p[1] * z_p_base_100[0]
   + 0.5 p[2] * z_p_base_100[0] - 0.5 p[3] * z_p_base_100[0]
   - 0.5 p[4] * z_p_base_100[0] ] <= 0
 infeasibility_p_base_100_row_1_leq_0: 5e-06 z_p_base_100[1] + [
   0.5 p[0] * z_p_base_100[1] - 0.5 p[1] * z_p_base_100[1]
   - 0.5 p[2] * z_p_base_100[1] + 0.5 p[3] * z_p_base_100[1]
   + 0.5 p[4] * z_p_base_100[1] ] <= 0
 infeasibility_p_base_100_row_2_leq_0: 5e-06 z_p_base_100[2] + [
   p[1] * z_p_base_100[2] ] <= 0
 infeasibility_p_base_100_row_3_leq_0: 5e-06 z_p_base_100[3] + [
   - 0.5 p[0] * z_p_base_100[3] + 0.5 p[1] * z_p_base_100[3]
   + 0.5 p[2] * z_p_base_100[3] + 0.5 p[3] * z_p_base_100[3]
   - 0.5 p[4] * z_p_base_100[3] ] <= 0
 infeasibility_p_base_100_row_4_leq_0: 5e-06 z_p_base_100[4] + [
   - p[1] * z_p_base_100[4] + p[4] * z_p_base_100[4] ] <= 0
 infeasibility_p_prime_base_100_row_0_leq_0: 5e-06 z_p_prime_base_100[0]
   + [ 0.5 alpha[0] * z_p_prime_base_100[0]
   - 0.5 alpha[1] * z_p_prime_base_100[0]
   - 0.5 alpha[2] * z_p_prime_base_100[0]
   + 0.5 alpha[3] * z_p_prime_base_100[0]
   - 0.5 alpha[4] * z_p_prime_base_100[0]
   + 0.5 p[0] * z_p_prime_base_100[0] + 0.5 p[1] * z_p_prime_base_100[0]
   + 0.5 p[2] * z_p_prime_base_100[0] - 0.5 p[3] * z_p_prime_base_100[0]
   - 0.5 p[4] * z_p_prime_base_100[0] ] <= 0
 infeasibility_p_prime_base_100_row_1_leq_0: 5e-06 z_p_prime_base_100[1]
   + [ 0.5 alpha[0] * z_p_prime_base_100[1]
   + 0.5 alpha[1] * z_p_prime_base_100[1]
   + 0.5 alpha[2] * z_p_prime_base_100[1]
   - 0.5 alpha[3] * z_p_prime_base_100[1]
   + 0.5 alpha[4] * z_p_prime_base_100[1]
   + 0.5 p[0] * z_p_prime_base_100[1] - 0.5 p[1] * z_p_prime_base_100[1]
   - 0.5 p[2] * z_p_prime_base_100[1] + 0.5 p[3] * z_p_prime_base_100[1]
   + 0.5 p[4] * z_p_prime_base_100[1] ] <= 0
 infeasibility_p_prime_base_100_row_2_leq_0: 5e-06 z_p_prime_base_100[2]
   + [ - alpha[1] * z_p_prime_base_100[2] + p[1] * z_p_prime_base_100[2] ]
   <= 0
 infeasibility_p_prime_base_100_row_3_leq_0: 5e-06 z_p_prime_base_100[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_100[3]
   - 0.5 alpha[1] * z_p_prime_base_100[3]
   - 0.5 alpha[2] * z_p_prime_base_100[3]
   - 0.5 alpha[3] * z_p_prime_base_100[3]
   - 0.5 alpha[4] * z_p_prime_base_100[3]
   - 0.5 p[0] * z_p_prime_base_100[3] + 0.5 p[1] * z_p_prime_base_100[3]
   + 0.5 p[2] * z_p_prime_base_100[3] + 0.5 p[3] * z_p_prime_base_100[3]
   - 0.5 p[4] * z_p_prime_base_100[3] ] <= 0
 infeasibility_p_prime_base_100_row_4_leq_0: 5e-06 z_p_prime_base_100[4]
   + [ alpha[1] * z_p_prime_base_100[4] + alpha[4] * z_p_prime_base_100[4]
   - p[1] * z_p_prime_base_100[4] + p[4] * z_p_prime_base_100[4] ] <= 0
 infeasibility_p_base_101_row_0_leq_0: 5e-06 z_p_base_101[0] + [
   0.5 p[0] * z_p_base_101[0] - 0.5 p[1] * z_p_base_101[0]
   + 0.5 p[2] * z_p_base_101[0] - 0.5 p[3] * z_p_base_101[0]
   + 0.5 p[4] * z_p_base_101[0] ] <= 0
 infeasibility_p_base_101_row_1_leq_0: 5e-06 z_p_base_101[1] + [
   0.5 p[0] * z_p_base_101[1] + 0.5 p[1] * z_p_base_101[1]
   - 0.5 p[2] * z_p_base_101[1] + 0.5 p[3] * z_p_base_101[1]
   - 0.5 p[4] * z_p_base_101[1] ] <= 0
 infeasibility_p_base_101_row_2_leq_0: 5e-06 z_p_base_101[2] + [
   p[1] * z_p_base_101[2] ] <= 0
 infeasibility_p_base_101_row_3_leq_0: 5e-06 z_p_base_101[3] + [
   - 0.5 p[0] * z_p_base_101[3] + 0.5 p[1] * z_p_base_101[3]
   + 0.5 p[2] * z_p_base_101[3] + 0.5 p[3] * z_p_base_101[3]
   - 0.5 p[4] * z_p_base_101[3] ] <= 0
 infeasibility_p_base_101_row_4_leq_0: 5e-06 z_p_base_101[4] + [
   - p[1] * z_p_base_101[4] + p[4] * z_p_base_101[4] ] <= 0
 infeasibility_p_prime_base_101_row_0_leq_0: 5e-06 z_p_prime_base_101[0]
   + [ 0.5 alpha[0] * z_p_prime_base_101[0]
   + 0.5 alpha[1] * z_p_prime_base_101[0]
   - 0.5 alpha[2] * z_p_prime_base_101[0]
   + 0.5 alpha[3] * z_p_prime_base_101[0]
   + 0.5 alpha[4] * z_p_prime_base_101[0]
   + 0.5 p[0] * z_p_prime_base_101[0] - 0.5 p[1] * z_p_prime_base_101[0]
   + 0.5 p[2] * z_p_prime_base_101[0] - 0.5 p[3] * z_p_prime_base_101[0]
   + 0.5 p[4] * z_p_prime_base_101[0] ] <= 0
 infeasibility_p_prime_base_101_row_1_leq_0: 5e-06 z_p_prime_base_101[1]
   + [ 0.5 alpha[0] * z_p_prime_base_101[1]
   - 0.5 alpha[1] * z_p_prime_base_101[1]
   + 0.5 alpha[2] * z_p_prime_base_101[1]
   - 0.5 alpha[3] * z_p_prime_base_101[1]
   - 0.5 alpha[4] * z_p_prime_base_101[1]
   + 0.5 p[0] * z_p_prime_base_101[1] + 0.5 p[1] * z_p_prime_base_101[1]
   - 0.5 p[2] * z_p_prime_base_101[1] + 0.5 p[3] * z_p_prime_base_101[1]
   - 0.5 p[4] * z_p_prime_base_101[1] ] <= 0
 infeasibility_p_prime_base_101_row_2_leq_0: 5e-06 z_p_prime_base_101[2]
   + [ - alpha[1] * z_p_prime_base_101[2] + p[1] * z_p_prime_base_101[2] ]
   <= 0
 infeasibility_p_prime_base_101_row_3_leq_0: 5e-06 z_p_prime_base_101[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_101[3]
   - 0.5 alpha[1] * z_p_prime_base_101[3]
   - 0.5 alpha[2] * z_p_prime_base_101[3]
   - 0.5 alpha[3] * z_p_prime_base_101[3]
   - 0.5 alpha[4] * z_p_prime_base_101[3]
   - 0.5 p[0] * z_p_prime_base_101[3] + 0.5 p[1] * z_p_prime_base_101[3]
   + 0.5 p[2] * z_p_prime_base_101[3] + 0.5 p[3] * z_p_prime_base_101[3]
   - 0.5 p[4] * z_p_prime_base_101[3] ] <= 0
 infeasibility_p_prime_base_101_row_4_leq_0: 5e-06 z_p_prime_base_101[4]
   + [ alpha[1] * z_p_prime_base_101[4] + alpha[4] * z_p_prime_base_101[4]
   - p[1] * z_p_prime_base_101[4] + p[4] * z_p_prime_base_101[4] ] <= 0
 infeasibility_p_base_102_row_0_leq_0: 5e-06 z_p_base_102[0] + [
   p[0] * z_p_base_102[0] - p[4] * z_p_base_102[0] ] <= 0
 infeasibility_p_base_102_row_1_leq_0: 5e-06 z_p_base_102[1] + [
   p[4] * z_p_base_102[1] ] <= 0
 infeasibility_p_base_102_row_2_leq_0: 5e-06 z_p_base_102[2] + [
   - 0.5 p[0] * z_p_base_102[2] + 0.5 p[1] * z_p_base_102[2]
   + 0.5 p[2] * z_p_base_102[2] - 0.5 p[3] * z_p_base_102[2]
   + 0.5 p[4] * z_p_base_102[2] ] <= 0
 infeasibility_p_base_102_row_3_leq_0: 5e-06 z_p_base_102[3] + [
   0.5 p[0] * z_p_base_102[3] + 0.5 p[1] * z_p_base_102[3]
   - 0.5 p[2] * z_p_base_102[3] + 0.5 p[3] * z_p_base_102[3]
   - 0.5 p[4] * z_p_base_102[3] ] <= 0
 infeasibility_p_base_102_row_4_leq_0: 5e-06 z_p_base_102[4] + [
   - 0.5 p[0] * z_p_base_102[4] - 0.5 p[1] * z_p_base_102[4]
   + 0.5 p[2] * z_p_base_102[4] + 0.5 p[3] * z_p_base_102[4]
   + 0.5 p[4] * z_p_base_102[4] ] <= 0
 infeasibility_p_prime_base_102_row_0_leq_0: 5e-06 z_p_prime_base_102[0]
   + [ alpha[0] * z_p_prime_base_102[0] - alpha[4] * z_p_prime_base_102[0]
   + p[0] * z_p_prime_base_102[0] - p[4] * z_p_prime_base_102[0] ] <= 0
 infeasibility_p_prime_base_102_row_1_leq_0: 5e-06 z_p_prime_base_102[1]
   + [ alpha[4] * z_p_prime_base_102[1] + p[4] * z_p_prime_base_102[1] ]
   <= 0
 infeasibility_p_prime_base_102_row_2_leq_0: 5e-06 z_p_prime_base_102[2]
   + [ - 0.5 alpha[0] * z_p_prime_base_102[2]
   - 0.5 alpha[1] * z_p_prime_base_102[2]
   - 0.5 alpha[2] * z_p_prime_base_102[2]
   + 0.5 alpha[3] * z_p_prime_base_102[2]
   + 0.5 alpha[4] * z_p_prime_base_102[2]
   - 0.5 p[0] * z_p_prime_base_102[2] + 0.5 p[1] * z_p_prime_base_102[2]
   + 0.5 p[2] * z_p_prime_base_102[2] - 0.5 p[3] * z_p_prime_base_102[2]
   + 0.5 p[4] * z_p_prime_base_102[2] ] <= 0
 infeasibility_p_prime_base_102_row_3_leq_0: 5e-06 z_p_prime_base_102[3]
   + [ 0.5 alpha[0] * z_p_prime_base_102[3]
   - 0.5 alpha[1] * z_p_prime_base_102[3]
   + 0.5 alpha[2] * z_p_prime_base_102[3]
   - 0.5 alpha[3] * z_p_prime_base_102[3]
   - 0.5 alpha[4] * z_p_prime_base_102[3]
   + 0.5 p[0] * z_p_prime_base_102[3] + 0.5 p[1] * z_p_prime_base_102[3]
   - 0.5 p[2] * z_p_prime_base_102[3] + 0.5 p[3] * z_p_prime_base_102[3]
   - 0.5 p[4] * z_p_prime_base_102[3] ] <= 0
 infeasibility_p_prime_base_102_row_4_leq_0: 5e-06 z_p_prime_base_102[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_102[4]
   + 0.5 alpha[1] * z_p_prime_base_102[4]
   - 0.5 alpha[2] * z_p_prime_base_102[4]
   - 0.5 alpha[3] * z_p_prime_base_102[4]
   + 0.5 alpha[4] * z_p_prime_base_102[4]
   - 0.5 p[0] * z_p_prime_base_102[4] - 0.5 p[1] * z_p_prime_base_102[4]
   + 0.5 p[2] * z_p_prime_base_102[4] + 0.5 p[3] * z_p_prime_base_102[4]
   + 0.5 p[4] * z_p_prime_base_102[4] ] <= 0
 infeasibility_p_base_103_row_0_leq_0: 5e-06 z_p_base_103[0] + [
   0.5 p[0] * z_p_base_103[0] - 0.5 p[1] * z_p_base_103[0]
   + 0.5 p[2] * z_p_base_103[0] + 0.5 p[3] * z_p_base_103[0]
   - 0.5 p[4] * z_p_base_103[0] ] <= 0
 infeasibility_p_base_103_row_1_leq_0: 5e-06 z_p_base_103[1] + [
   0.5 p[0] * z_p_base_103[1] + 0.5 p[1] * z_p_base_103[1]
   - 0.5 p[2] * z_p_base_103[1] - 0.5 p[3] * z_p_base_103[1]
   + 0.5 p[4] * z_p_base_103[1] ] <= 0
 infeasibility_p_base_103_row_2_leq_0: 5e-06 z_p_base_103[2] + [
   p[1] * z_p_base_103[2] - p[3] * z_p_base_103[2] ] <= 0
 infeasibility_p_base_103_row_3_leq_0: 5e-06 z_p_base_103[3] + [
   p[3] * z_p_base_103[3] ] <= 0
 infeasibility_p_base_103_row_4_leq_0: 5e-06 z_p_base_103[4] + [
   - 0.5 p[0] * z_p_base_103[4] - 0.5 p[1] * z_p_base_103[4]
   + 0.5 p[2] * z_p_base_103[4] + 0.5 p[3] * z_p_base_103[4]
   + 0.5 p[4] * z_p_base_103[4] ] <= 0
 infeasibility_p_prime_base_103_row_0_leq_0: 5e-06 z_p_prime_base_103[0]
   + [ 0.5 alpha[0] * z_p_prime_base_103[0]
   + 0.5 alpha[1] * z_p_prime_base_103[0]
   - 0.5 alpha[2] * z_p_prime_base_103[0]
   - 0.5 alpha[3] * z_p_prime_base_103[0]
   - 0.5 alpha[4] * z_p_prime_base_103[0]
   + 0.5 p[0] * z_p_prime_base_103[0] - 0.5 p[1] * z_p_prime_base_103[0]
   + 0.5 p[2] * z_p_prime_base_103[0] + 0.5 p[3] * z_p_prime_base_103[0]
   - 0.5 p[4] * z_p_prime_base_103[0] ] <= 0
 infeasibility_p_prime_base_103_row_1_leq_0: 5e-06 z_p_prime_base_103[1]
   + [ 0.5 alpha[0] * z_p_prime_base_103[1]
   - 0.5 alpha[1] * z_p_prime_base_103[1]
   + 0.5 alpha[2] * z_p_prime_base_103[1]
   + 0.5 alpha[3] * z_p_prime_base_103[1]
   + 0.5 alpha[4] * z_p_prime_base_103[1]
   + 0.5 p[0] * z_p_prime_base_103[1] + 0.5 p[1] * z_p_prime_base_103[1]
   - 0.5 p[2] * z_p_prime_base_103[1] - 0.5 p[3] * z_p_prime_base_103[1]
   + 0.5 p[4] * z_p_prime_base_103[1] ] <= 0
 infeasibility_p_prime_base_103_row_2_leq_0: 5e-06 z_p_prime_base_103[2]
   + [ - alpha[1] * z_p_prime_base_103[2]
   + alpha[3] * z_p_prime_base_103[2] + p[1] * z_p_prime_base_103[2]
   - p[3] * z_p_prime_base_103[2] ] <= 0
 infeasibility_p_prime_base_103_row_3_leq_0: 5e-06 z_p_prime_base_103[3]
   + [ - alpha[3] * z_p_prime_base_103[3] + p[3] * z_p_prime_base_103[3] ]
   <= 0
 infeasibility_p_prime_base_103_row_4_leq_0: 5e-06 z_p_prime_base_103[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_103[4]
   + 0.5 alpha[1] * z_p_prime_base_103[4]
   - 0.5 alpha[2] * z_p_prime_base_103[4]
   - 0.5 alpha[3] * z_p_prime_base_103[4]
   + 0.5 alpha[4] * z_p_prime_base_103[4]
   - 0.5 p[0] * z_p_prime_base_103[4] - 0.5 p[1] * z_p_prime_base_103[4]
   + 0.5 p[2] * z_p_prime_base_103[4] + 0.5 p[3] * z_p_prime_base_103[4]
   + 0.5 p[4] * z_p_prime_base_103[4] ] <= 0
 infeasibility_p_base_104_row_0_leq_0: 5e-06 z_p_base_104[0] + [
   0.5 p[0] * z_p_base_104[0] - 0.5 p[1] * z_p_base_104[0]
   + 0.5 p[2] * z_p_base_104[0] + 0.5 p[3] * z_p_base_104[0]
   - 0.5 p[4] * z_p_base_104[0] ] <= 0
 infeasibility_p_base_104_row_1_leq_0: 5e-06 z_p_base_104[1] + [
   0.5 p[0] * z_p_base_104[1] + 0.5 p[1] * z_p_base_104[1]
   - 0.5 p[2] * z_p_base_104[1] - 0.5 p[3] * z_p_base_104[1]
   + 0.5 p[4] * z_p_base_104[1] ] <= 0
 infeasibility_p_base_104_row_2_leq_0: 5e-06 z_p_base_104[2] + [
   - 0.5 p[0] * z_p_base_104[2] + 0.5 p[1] * z_p_base_104[2]
   + 0.5 p[2] * z_p_base_104[2] - 0.5 p[3] * z_p_base_104[2]
   + 0.5 p[4] * z_p_base_104[2] ] <= 0
 infeasibility_p_base_104_row_3_leq_0: 5e-06 z_p_base_104[3] + [
   0.5 p[0] * z_p_base_104[3] + 0.5 p[1] * z_p_base_104[3]
   - 0.5 p[2] * z_p_base_104[3] + 0.5 p[3] * z_p_base_104[3]
   - 0.5 p[4] * z_p_base_104[3] ] <= 0
 infeasibility_p_base_104_row_4_leq_0: 5e-06 z_p_base_104[4] + [
   - 0.5 p[0] * z_p_base_104[4] - 0.5 p[1] * z_p_base_104[4]
   + 0.5 p[2] * z_p_base_104[4] + 0.5 p[3] * z_p_base_104[4]
   + 0.5 p[4] * z_p_base_104[4] ] <= 0
 infeasibility_p_prime_base_104_row_0_leq_0: 5e-06 z_p_prime_base_104[0]
   + [ 0.5 alpha[0] * z_p_prime_base_104[0]
   + 0.5 alpha[1] * z_p_prime_base_104[0]
   - 0.5 alpha[2] * z_p_prime_base_104[0]
   - 0.5 alpha[3] * z_p_prime_base_104[0]
   - 0.5 alpha[4] * z_p_prime_base_104[0]
   + 0.5 p[0] * z_p_prime_base_104[0] - 0.5 p[1] * z_p_prime_base_104[0]
   + 0.5 p[2] * z_p_prime_base_104[0] + 0.5 p[3] * z_p_prime_base_104[0]
   - 0.5 p[4] * z_p_prime_base_104[0] ] <= 0
 infeasibility_p_prime_base_104_row_1_leq_0: 5e-06 z_p_prime_base_104[1]
   + [ 0.5 alpha[0] * z_p_prime_base_104[1]
   - 0.5 alpha[1] * z_p_prime_base_104[1]
   + 0.5 alpha[2] * z_p_prime_base_104[1]
   + 0.5 alpha[3] * z_p_prime_base_104[1]
   + 0.5 alpha[4] * z_p_prime_base_104[1]
   + 0.5 p[0] * z_p_prime_base_104[1] + 0.5 p[1] * z_p_prime_base_104[1]
   - 0.5 p[2] * z_p_prime_base_104[1] - 0.5 p[3] * z_p_prime_base_104[1]
   + 0.5 p[4] * z_p_prime_base_104[1] ] <= 0
 infeasibility_p_prime_base_104_row_2_leq_0: 5e-06 z_p_prime_base_104[2]
   + [ - 0.5 alpha[0] * z_p_prime_base_104[2]
   - 0.5 alpha[1] * z_p_prime_base_104[2]
   - 0.5 alpha[2] * z_p_prime_base_104[2]
   + 0.5 alpha[3] * z_p_prime_base_104[2]
   + 0.5 alpha[4] * z_p_prime_base_104[2]
   - 0.5 p[0] * z_p_prime_base_104[2] + 0.5 p[1] * z_p_prime_base_104[2]
   + 0.5 p[2] * z_p_prime_base_104[2] - 0.5 p[3] * z_p_prime_base_104[2]
   + 0.5 p[4] * z_p_prime_base_104[2] ] <= 0
 infeasibility_p_prime_base_104_row_3_leq_0: 5e-06 z_p_prime_base_104[3]
   + [ 0.5 alpha[0] * z_p_prime_base_104[3]
   - 0.5 alpha[1] * z_p_prime_base_104[3]
   + 0.5 alpha[2] * z_p_prime_base_104[3]
   - 0.5 alpha[3] * z_p_prime_base_104[3]
   - 0.5 alpha[4] * z_p_prime_base_104[3]
   + 0.5 p[0] * z_p_prime_base_104[3] + 0.5 p[1] * z_p_prime_base_104[3]
   - 0.5 p[2] * z_p_prime_base_104[3] + 0.5 p[3] * z_p_prime_base_104[3]
   - 0.5 p[4] * z_p_prime_base_104[3] ] <= 0
 infeasibility_p_prime_base_104_row_4_leq_0: 5e-06 z_p_prime_base_104[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_104[4]
   + 0.5 alpha[1] * z_p_prime_base_104[4]
   - 0.5 alpha[2] * z_p_prime_base_104[4]
   - 0.5 alpha[3] * z_p_prime_base_104[4]
   + 0.5 alpha[4] * z_p_prime_base_104[4]
   - 0.5 p[0] * z_p_prime_base_104[4] - 0.5 p[1] * z_p_prime_base_104[4]
   + 0.5 p[2] * z_p_prime_base_104[4] + 0.5 p[3] * z_p_prime_base_104[4]
   + 0.5 p[4] * z_p_prime_base_104[4] ] <= 0
 infeasibility_p_base_105_row_0_leq_0: 5e-06 z_p_base_105[0] + [
   0.5 p[0] * z_p_base_105[0] - 0.5 p[1] * z_p_base_105[0]
   + 0.5 p[2] * z_p_base_105[0] - 0.5 p[3] * z_p_base_105[0]
   - 0.5 p[4] * z_p_base_105[0] ] <= 0
 infeasibility_p_base_105_row_1_leq_0: 5e-06 z_p_base_105[1] + [
   0.5 p[0] * z_p_base_105[1] + 0.5 p[1] * z_p_base_105[1]
   - 0.5 p[2] * z_p_base_105[1] + 0.5 p[3] * z_p_base_105[1]
   + 0.5 p[4] * z_p_base_105[1] ] <= 0
 infeasibility_p_base_105_row_2_leq_0: 5e-06 z_p_base_105[2] + [
   p[1] * z_p_base_105[2] ] <= 0
 infeasibility_p_base_105_row_3_leq_0: 5e-06 z_p_base_105[3] + [
   p[3] * z_p_base_105[3] ] <= 0
 infeasibility_p_base_105_row_4_leq_0: 5e-06 z_p_base_105[4] + [
   - 0.5 p[0] * z_p_base_105[4] - 0.5 p[1] * z_p_base_105[4]
   + 0.5 p[2] * z_p_base_105[4] - 0.5 p[3] * z_p_base_105[4]
   + 0.5 p[4] * z_p_base_105[4] ] <= 0
 infeasibility_p_prime_base_105_row_0_leq_0: 5e-06 z_p_prime_base_105[0]
   + [ 0.5 alpha[0] * z_p_prime_base_105[0]
   + 0.5 alpha[1] * z_p_prime_base_105[0]
   - 0.5 alpha[2] * z_p_prime_base_105[0]
   + 0.5 alpha[3] * z_p_prime_base_105[0]
   - 0.5 alpha[4] * z_p_prime_base_105[0]
   + 0.5 p[0] * z_p_prime_base_105[0] - 0.5 p[1] * z_p_prime_base_105[0]
   + 0.5 p[2] * z_p_prime_base_105[0] - 0.5 p[3] * z_p_prime_base_105[0]
   - 0.5 p[4] * z_p_prime_base_105[0] ] <= 0
 infeasibility_p_prime_base_105_row_1_leq_0: 5e-06 z_p_prime_base_105[1]
   + [ 0.5 alpha[0] * z_p_prime_base_105[1]
   - 0.5 alpha[1] * z_p_prime_base_105[1]
   + 0.5 alpha[2] * z_p_prime_base_105[1]
   - 0.5 alpha[3] * z_p_prime_base_105[1]
   + 0.5 alpha[4] * z_p_prime_base_105[1]
   + 0.5 p[0] * z_p_prime_base_105[1] + 0.5 p[1] * z_p_prime_base_105[1]
   - 0.5 p[2] * z_p_prime_base_105[1] + 0.5 p[3] * z_p_prime_base_105[1]
   + 0.5 p[4] * z_p_prime_base_105[1] ] <= 0
 infeasibility_p_prime_base_105_row_2_leq_0: 5e-06 z_p_prime_base_105[2]
   + [ - alpha[1] * z_p_prime_base_105[2] + p[1] * z_p_prime_base_105[2] ]
   <= 0
 infeasibility_p_prime_base_105_row_3_leq_0: 5e-06 z_p_prime_base_105[3]
   + [ - alpha[3] * z_p_prime_base_105[3] + p[3] * z_p_prime_base_105[3] ]
   <= 0
 infeasibility_p_prime_base_105_row_4_leq_0: 5e-06 z_p_prime_base_105[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_105[4]
   + 0.5 alpha[1] * z_p_prime_base_105[4]
   - 0.5 alpha[2] * z_p_prime_base_105[4]
   + 0.5 alpha[3] * z_p_prime_base_105[4]
   + 0.5 alpha[4] * z_p_prime_base_105[4]
   - 0.5 p[0] * z_p_prime_base_105[4] - 0.5 p[1] * z_p_prime_base_105[4]
   + 0.5 p[2] * z_p_prime_base_105[4] - 0.5 p[3] * z_p_prime_base_105[4]
   + 0.5 p[4] * z_p_prime_base_105[4] ] <= 0
 infeasibility_p_base_106_row_0_leq_0: 5e-06 z_p_base_106[0] + [
   0.5 p[0] * z_p_base_106[0] - 0.5 p[1] * z_p_base_106[0]
   + 0.5 p[2] * z_p_base_106[0] + 0.5 p[3] * z_p_base_106[0]
   - 0.5 p[4] * z_p_base_106[0] ] <= 0
 infeasibility_p_base_106_row_1_leq_0: 5e-06 z_p_base_106[1] + [
   0.5 p[0] * z_p_base_106[1] + 0.5 p[1] * z_p_base_106[1]
   - 0.5 p[2] * z_p_base_106[1] - 0.5 p[3] * z_p_base_106[1]
   + 0.5 p[4] * z_p_base_106[1] ] <= 0
 infeasibility_p_base_106_row_2_leq_0: 5e-06 z_p_base_106[2] + [
   p[1] * z_p_base_106[2] ] <= 0
 infeasibility_p_base_106_row_3_leq_0: 5e-06 z_p_base_106[3] + [
   - 0.5 p[0] * z_p_base_106[3] - 0.5 p[1] * z_p_base_106[3]
   + 0.5 p[2] * z_p_base_106[3] - 0.5 p[3] * z_p_base_106[3]
   + 0.5 p[4] * z_p_base_106[3] ] <= 0
 infeasibility_p_base_106_row_4_leq_0: 5e-06 z_p_base_106[4] + [
   p[3] * z_p_base_106[4] ] <= 0
 infeasibility_p_prime_base_106_row_0_leq_0: 5e-06 z_p_prime_base_106[0]
   + [ 0.5 alpha[0] * z_p_prime_base_106[0]
   + 0.5 alpha[1] * z_p_prime_base_106[0]
   - 0.5 alpha[2] * z_p_prime_base_106[0]
   - 0.5 alpha[3] * z_p_prime_base_106[0]
   - 0.5 alpha[4] * z_p_prime_base_106[0]
   + 0.5 p[0] * z_p_prime_base_106[0] - 0.5 p[1] * z_p_prime_base_106[0]
   + 0.5 p[2] * z_p_prime_base_106[0] + 0.5 p[3] * z_p_prime_base_106[0]
   - 0.5 p[4] * z_p_prime_base_106[0] ] <= 0
 infeasibility_p_prime_base_106_row_1_leq_0: 5e-06 z_p_prime_base_106[1]
   + [ 0.5 alpha[0] * z_p_prime_base_106[1]
   - 0.5 alpha[1] * z_p_prime_base_106[1]
   + 0.5 alpha[2] * z_p_prime_base_106[1]
   + 0.5 alpha[3] * z_p_prime_base_106[1]
   + 0.5 alpha[4] * z_p_prime_base_106[1]
   + 0.5 p[0] * z_p_prime_base_106[1] + 0.5 p[1] * z_p_prime_base_106[1]
   - 0.5 p[2] * z_p_prime_base_106[1] - 0.5 p[3] * z_p_prime_base_106[1]
   + 0.5 p[4] * z_p_prime_base_106[1] ] <= 0
 infeasibility_p_prime_base_106_row_2_leq_0: 5e-06 z_p_prime_base_106[2]
   + [ - alpha[1] * z_p_prime_base_106[2] + p[1] * z_p_prime_base_106[2] ]
   <= 0
 infeasibility_p_prime_base_106_row_3_leq_0: 5e-06 z_p_prime_base_106[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_106[3]
   + 0.5 alpha[1] * z_p_prime_base_106[3]
   - 0.5 alpha[2] * z_p_prime_base_106[3]
   + 0.5 alpha[3] * z_p_prime_base_106[3]
   + 0.5 alpha[4] * z_p_prime_base_106[3]
   - 0.5 p[0] * z_p_prime_base_106[3] - 0.5 p[1] * z_p_prime_base_106[3]
   + 0.5 p[2] * z_p_prime_base_106[3] - 0.5 p[3] * z_p_prime_base_106[3]
   + 0.5 p[4] * z_p_prime_base_106[3] ] <= 0
 infeasibility_p_prime_base_106_row_4_leq_0: 5e-06 z_p_prime_base_106[4]
   + [ - alpha[3] * z_p_prime_base_106[4] + p[3] * z_p_prime_base_106[4] ]
   <= 0
 infeasibility_p_base_107_row_0_leq_0: 5e-06 z_p_base_107[0] + [
   0.5 p[0] * z_p_base_107[0] + 0.5 p[1] * z_p_base_107[0]
   + 0.5 p[2] * z_p_base_107[0] - 0.5 p[3] * z_p_base_107[0]
   - 0.5 p[4] * z_p_base_107[0] ] <= 0
 infeasibility_p_base_107_row_1_leq_0: 5e-06 z_p_base_107[1] + [
   0.5 p[0] * z_p_base_107[1] - 0.5 p[1] * z_p_base_107[1]
   - 0.5 p[2] * z_p_base_107[1] + 0.5 p[3] * z_p_base_107[1]
   + 0.5 p[4] * z_p_base_107[1] ] <= 0
 infeasibility_p_base_107_row_2_leq_0: 5e-06 z_p_base_107[2] + [
   0.5 p[0] * z_p_base_107[2] + 0.5 p[1] * z_p_base_107[2]
   - 0.5 p[2] * z_p_base_107[2] + 0.5 p[3] * z_p_base_107[2]
   - 0.5 p[4] * z_p_base_107[2] ] <= 0
 infeasibility_p_base_107_row_3_leq_0: 5e-06 z_p_base_107[3] + [
   - 0.5 p[0] * z_p_base_107[3] + 0.5 p[1] * z_p_base_107[3]
   + 0.5 p[2] * z_p_base_107[3] - 0.5 p[3] * z_p_base_107[3]
   + 0.5 p[4] * z_p_base_107[3] ] <= 0
 infeasibility_p_base_107_row_4_leq_0: 5e-06 z_p_base_107[4] + [
   - 0.5 p[0] * z_p_base_107[4] - 0.5 p[1] * z_p_base_107[4]
   + 0.5 p[2] * z_p_base_107[4] + 0.5 p[3] * z_p_base_107[4]
   + 0.5 p[4] * z_p_base_107[4] ] <= 0
 infeasibility_p_prime_base_107_row_0_leq_0: 5e-06 z_p_prime_base_107[0]
   + [ 0.5 alpha[0] * z_p_prime_base_107[0]
   - 0.5 alpha[1] * z_p_prime_base_107[0]
   - 0.5 alpha[2] * z_p_prime_base_107[0]
   + 0.5 alpha[3] * z_p_prime_base_107[0]
   - 0.5 alpha[4] * z_p_prime_base_107[0]
   + 0.5 p[0] * z_p_prime_base_107[0] + 0.5 p[1] * z_p_prime_base_107[0]
   + 0.5 p[2] * z_p_prime_base_107[0] - 0.5 p[3] * z_p_prime_base_107[0]
   - 0.5 p[4] * z_p_prime_base_107[0] ] <= 0
 infeasibility_p_prime_base_107_row_1_leq_0: 5e-06 z_p_prime_base_107[1]
   + [ 0.5 alpha[0] * z_p_prime_base_107[1]
   + 0.5 alpha[1] * z_p_prime_base_107[1]
   + 0.5 alpha[2] * z_p_prime_base_107[1]
   - 0.5 alpha[3] * z_p_prime_base_107[1]
   + 0.5 alpha[4] * z_p_prime_base_107[1]
   + 0.5 p[0] * z_p_prime_base_107[1] - 0.5 p[1] * z_p_prime_base_107[1]
   - 0.5 p[2] * z_p_prime_base_107[1] + 0.5 p[3] * z_p_prime_base_107[1]
   + 0.5 p[4] * z_p_prime_base_107[1] ] <= 0
 infeasibility_p_prime_base_107_row_2_leq_0: 5e-06 z_p_prime_base_107[2]
   + [ 0.5 alpha[0] * z_p_prime_base_107[2]
   - 0.5 alpha[1] * z_p_prime_base_107[2]
   + 0.5 alpha[2] * z_p_prime_base_107[2]
   - 0.5 alpha[3] * z_p_prime_base_107[2]
   - 0.5 alpha[4] * z_p_prime_base_107[2]
   + 0.5 p[0] * z_p_prime_base_107[2] + 0.5 p[1] * z_p_prime_base_107[2]
   - 0.5 p[2] * z_p_prime_base_107[2] + 0.5 p[3] * z_p_prime_base_107[2]
   - 0.5 p[4] * z_p_prime_base_107[2] ] <= 0
 infeasibility_p_prime_base_107_row_3_leq_0: 5e-06 z_p_prime_base_107[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_107[3]
   - 0.5 alpha[1] * z_p_prime_base_107[3]
   - 0.5 alpha[2] * z_p_prime_base_107[3]
   + 0.5 alpha[3] * z_p_prime_base_107[3]
   + 0.5 alpha[4] * z_p_prime_base_107[3]
   - 0.5 p[0] * z_p_prime_base_107[3] + 0.5 p[1] * z_p_prime_base_107[3]
   + 0.5 p[2] * z_p_prime_base_107[3] - 0.5 p[3] * z_p_prime_base_107[3]
   + 0.5 p[4] * z_p_prime_base_107[3] ] <= 0
 infeasibility_p_prime_base_107_row_4_leq_0: 5e-06 z_p_prime_base_107[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_107[4]
   + 0.5 alpha[1] * z_p_prime_base_107[4]
   - 0.5 alpha[2] * z_p_prime_base_107[4]
   - 0.5 alpha[3] * z_p_prime_base_107[4]
   + 0.5 alpha[4] * z_p_prime_base_107[4]
   - 0.5 p[0] * z_p_prime_base_107[4] - 0.5 p[1] * z_p_prime_base_107[4]
   + 0.5 p[2] * z_p_prime_base_107[4] + 0.5 p[3] * z_p_prime_base_107[4]
   + 0.5 p[4] * z_p_prime_base_107[4] ] <= 0
 infeasibility_p_base_108_row_0_leq_0: 5e-06 z_p_base_108[0] + [
   0.5 p[0] * z_p_base_108[0] + 0.5 p[1] * z_p_base_108[0]
   + 0.5 p[2] * z_p_base_108[0] - 0.5 p[3] * z_p_base_108[0]
   - 0.5 p[4] * z_p_base_108[0] ] <= 0
 infeasibility_p_base_108_row_1_leq_0: 5e-06 z_p_base_108[1] + [
   0.5 p[0] * z_p_base_108[1] - 0.5 p[1] * z_p_base_108[1]
   - 0.5 p[2] * z_p_base_108[1] + 0.5 p[3] * z_p_base_108[1]
   + 0.5 p[4] * z_p_base_108[1] ] <= 0
 infeasibility_p_base_108_row_2_leq_0: 5e-06 z_p_base_108[2] + [
   p[3] * z_p_base_108[2] ] <= 0
 infeasibility_p_base_108_row_3_leq_0: 5e-06 z_p_base_108[3] + [
   p[1] * z_p_base_108[3] - p[3] * z_p_base_108[3] ] <= 0
 infeasibility_p_base_108_row_4_leq_0: 5e-06 z_p_base_108[4] + [
   - 0.5 p[0] * z_p_base_108[4] - 0.5 p[1] * z_p_base_108[4]
   + 0.5 p[2] * z_p_base_108[4] + 0.5 p[3] * z_p_base_108[4]
   + 0.5 p[4] * z_p_base_108[4] ] <= 0
 infeasibility_p_prime_base_108_row_0_leq_0: 5e-06 z_p_prime_base_108[0]
   + [ 0.5 alpha[0] * z_p_prime_base_108[0]
   - 0.5 alpha[1] * z_p_prime_base_108[0]
   - 0.5 alpha[2] * z_p_prime_base_108[0]
   + 0.5 alpha[3] * z_p_prime_base_108[0]
   - 0.5 alpha[4] * z_p_prime_base_108[0]
   + 0.5 p[0] * z_p_prime_base_108[0] + 0.5 p[1] * z_p_prime_base_108[0]
   + 0.5 p[2] * z_p_prime_base_108[0] - 0.5 p[3] * z_p_prime_base_108[0]
   - 0.5 p[4] * z_p_prime_base_108[0] ] <= 0
 infeasibility_p_prime_base_108_row_1_leq_0: 5e-06 z_p_prime_base_108[1]
   + [ 0.5 alpha[0] * z_p_prime_base_108[1]
   + 0.5 alpha[1] * z_p_prime_base_108[1]
   + 0.5 alpha[2] * z_p_prime_base_108[1]
   - 0.5 alpha[3] * z_p_prime_base_108[1]
   + 0.5 alpha[4] * z_p_prime_base_108[1]
   + 0.5 p[0] * z_p_prime_base_108[1] - 0.5 p[1] * z_p_prime_base_108[1]
   - 0.5 p[2] * z_p_prime_base_108[1] + 0.5 p[3] * z_p_prime_base_108[1]
   + 0.5 p[4] * z_p_prime_base_108[1] ] <= 0
 infeasibility_p_prime_base_108_row_2_leq_0: 5e-06 z_p_prime_base_108[2]
   + [ - alpha[3] * z_p_prime_base_108[2] + p[3] * z_p_prime_base_108[2] ]
   <= 0
 infeasibility_p_prime_base_108_row_3_leq_0: 5e-06 z_p_prime_base_108[3]
   + [ - alpha[1] * z_p_prime_base_108[3]
   + alpha[3] * z_p_prime_base_108[3] + p[1] * z_p_prime_base_108[3]
   - p[3] * z_p_prime_base_108[3] ] <= 0
 infeasibility_p_prime_base_108_row_4_leq_0: 5e-06 z_p_prime_base_108[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_108[4]
   + 0.5 alpha[1] * z_p_prime_base_108[4]
   - 0.5 alpha[2] * z_p_prime_base_108[4]
   - 0.5 alpha[3] * z_p_prime_base_108[4]
   + 0.5 alpha[4] * z_p_prime_base_108[4]
   - 0.5 p[0] * z_p_prime_base_108[4] - 0.5 p[1] * z_p_prime_base_108[4]
   + 0.5 p[2] * z_p_prime_base_108[4] + 0.5 p[3] * z_p_prime_base_108[4]
   + 0.5 p[4] * z_p_prime_base_108[4] ] <= 0
 infeasibility_p_base_109_row_0_leq_0: 5e-06 z_p_base_109[0] + [
   p[2] * z_p_base_109[0] ] <= 0
 infeasibility_p_base_109_row_1_leq_0: 5e-06 z_p_base_109[1] + [
   p[0] * z_p_base_109[1] - p[2] * z_p_base_109[1] ] <= 0
 infeasibility_p_base_109_row_2_leq_0: 5e-06 z_p_base_109[2] + [
   0.5 p[0] * z_p_base_109[2] + 0.5 p[1] * z_p_base_109[2]
   - 0.5 p[2] * z_p_base_109[2] + 0.5 p[3] * z_p_base_109[2]
   - 0.5 p[4] * z_p_base_109[2] ] <= 0
 infeasibility_p_base_109_row_3_leq_0: 5e-06 z_p_base_109[3] + [
   - 0.5 p[0] * z_p_base_109[3] + 0.5 p[1] * z_p_base_109[3]
   + 0.5 p[2] * z_p_base_109[3] - 0.5 p[3] * z_p_base_109[3]
   + 0.5 p[4] * z_p_base_109[3] ] <= 0
 infeasibility_p_base_109_row_4_leq_0: 5e-06 z_p_base_109[4] + [
   - 0.5 p[0] * z_p_base_109[4] - 0.5 p[1] * z_p_base_109[4]
   + 0.5 p[2] * z_p_base_109[4] + 0.5 p[3] * z_p_base_109[4]
   + 0.5 p[4] * z_p_base_109[4] ] <= 0
 infeasibility_p_prime_base_109_row_0_leq_0: 5e-06 z_p_prime_base_109[0]
   + [ - alpha[2] * z_p_prime_base_109[0] + p[2] * z_p_prime_base_109[0] ]
   <= 0
 infeasibility_p_prime_base_109_row_1_leq_0: 5e-06 z_p_prime_base_109[1]
   + [ alpha[0] * z_p_prime_base_109[1] + alpha[2] * z_p_prime_base_109[1]
   + p[0] * z_p_prime_base_109[1] - p[2] * z_p_prime_base_109[1] ] <= 0
 infeasibility_p_prime_base_109_row_2_leq_0: 5e-06 z_p_prime_base_109[2]
   + [ 0.5 alpha[0] * z_p_prime_base_109[2]
   - 0.5 alpha[1] * z_p_prime_base_109[2]
   + 0.5 alpha[2] * z_p_prime_base_109[2]
   - 0.5 alpha[3] * z_p_prime_base_109[2]
   - 0.5 alpha[4] * z_p_prime_base_109[2]
   + 0.5 p[0] * z_p_prime_base_109[2] + 0.5 p[1] * z_p_prime_base_109[2]
   - 0.5 p[2] * z_p_prime_base_109[2] + 0.5 p[3] * z_p_prime_base_109[2]
   - 0.5 p[4] * z_p_prime_base_109[2] ] <= 0
 infeasibility_p_prime_base_109_row_3_leq_0: 5e-06 z_p_prime_base_109[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_109[3]
   - 0.5 alpha[1] * z_p_prime_base_109[3]
   - 0.5 alpha[2] * z_p_prime_base_109[3]
   + 0.5 alpha[3] * z_p_prime_base_109[3]
   + 0.5 alpha[4] * z_p_prime_base_109[3]
   - 0.5 p[0] * z_p_prime_base_109[3] + 0.5 p[1] * z_p_prime_base_109[3]
   + 0.5 p[2] * z_p_prime_base_109[3] - 0.5 p[3] * z_p_prime_base_109[3]
   + 0.5 p[4] * z_p_prime_base_109[3] ] <= 0
 infeasibility_p_prime_base_109_row_4_leq_0: 5e-06 z_p_prime_base_109[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_109[4]
   + 0.5 alpha[1] * z_p_prime_base_109[4]
   - 0.5 alpha[2] * z_p_prime_base_109[4]
   - 0.5 alpha[3] * z_p_prime_base_109[4]
   + 0.5 alpha[4] * z_p_prime_base_109[4]
   - 0.5 p[0] * z_p_prime_base_109[4] - 0.5 p[1] * z_p_prime_base_109[4]
   + 0.5 p[2] * z_p_prime_base_109[4] + 0.5 p[3] * z_p_prime_base_109[4]
   + 0.5 p[4] * z_p_prime_base_109[4] ] <= 0
 infeasibility_p_base_110_row_0_leq_0: 5e-06 z_p_base_110[0] + [
   0.5 p[0] * z_p_base_110[0] + 0.5 p[1] * z_p_base_110[0]
   + 0.5 p[2] * z_p_base_110[0] - 0.5 p[3] * z_p_base_110[0]
   - 0.5 p[4] * z_p_base_110[0] ] <= 0
 infeasibility_p_base_110_row_1_leq_0: 5e-06 z_p_base_110[1] + [
   0.5 p[0] * z_p_base_110[1] - 0.5 p[1] * z_p_base_110[1]
   - 0.5 p[2] * z_p_base_110[1] + 0.5 p[3] * z_p_base_110[1]
   + 0.5 p[4] * z_p_base_110[1] ] <= 0
 infeasibility_p_base_110_row_2_leq_0: 5e-06 z_p_base_110[2] + [
   p[1] * z_p_base_110[2] ] <= 0
 infeasibility_p_base_110_row_3_leq_0: 5e-06 z_p_base_110[3] + [
   - p[1] * z_p_base_110[3] + p[3] * z_p_base_110[3] ] <= 0
 infeasibility_p_base_110_row_4_leq_0: 5e-06 z_p_base_110[4] + [
   - 0.5 p[0] * z_p_base_110[4] + 0.5 p[1] * z_p_base_110[4]
   + 0.5 p[2] * z_p_base_110[4] - 0.5 p[3] * z_p_base_110[4]
   + 0.5 p[4] * z_p_base_110[4] ] <= 0
 infeasibility_p_prime_base_110_row_0_leq_0: 5e-06 z_p_prime_base_110[0]
   + [ 0.5 alpha[0] * z_p_prime_base_110[0]
   - 0.5 alpha[1] * z_p_prime_base_110[0]
   - 0.5 alpha[2] * z_p_prime_base_110[0]
   + 0.5 alpha[3] * z_p_prime_base_110[0]
   - 0.5 alpha[4] * z_p_prime_base_110[0]
   + 0.5 p[0] * z_p_prime_base_110[0] + 0.5 p[1] * z_p_prime_base_110[0]
   + 0.5 p[2] * z_p_prime_base_110[0] - 0.5 p[3] * z_p_prime_base_110[0]
   - 0.5 p[4] * z_p_prime_base_110[0] ] <= 0
 infeasibility_p_prime_base_110_row_1_leq_0: 5e-06 z_p_prime_base_110[1]
   + [ 0.5 alpha[0] * z_p_prime_base_110[1]
   + 0.5 alpha[1] * z_p_prime_base_110[1]
   + 0.5 alpha[2] * z_p_prime_base_110[1]
   - 0.5 alpha[3] * z_p_prime_base_110[1]
   + 0.5 alpha[4] * z_p_prime_base_110[1]
   + 0.5 p[0] * z_p_prime_base_110[1] - 0.5 p[1] * z_p_prime_base_110[1]
   - 0.5 p[2] * z_p_prime_base_110[1] + 0.5 p[3] * z_p_prime_base_110[1]
   + 0.5 p[4] * z_p_prime_base_110[1] ] <= 0
 infeasibility_p_prime_base_110_row_2_leq_0: 5e-06 z_p_prime_base_110[2]
   + [ - alpha[1] * z_p_prime_base_110[2] + p[1] * z_p_prime_base_110[2] ]
   <= 0
 infeasibility_p_prime_base_110_row_3_leq_0: 5e-06 z_p_prime_base_110[3]
   + [ alpha[1] * z_p_prime_base_110[3] - alpha[3] * z_p_prime_base_110[3]
   - p[1] * z_p_prime_base_110[3] + p[3] * z_p_prime_base_110[3] ] <= 0
 infeasibility_p_prime_base_110_row_4_leq_0: 5e-06 z_p_prime_base_110[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_110[4]
   - 0.5 alpha[1] * z_p_prime_base_110[4]
   - 0.5 alpha[2] * z_p_prime_base_110[4]
   + 0.5 alpha[3] * z_p_prime_base_110[4]
   + 0.5 alpha[4] * z_p_prime_base_110[4]
   - 0.5 p[0] * z_p_prime_base_110[4] + 0.5 p[1] * z_p_prime_base_110[4]
   + 0.5 p[2] * z_p_prime_base_110[4] - 0.5 p[3] * z_p_prime_base_110[4]
   + 0.5 p[4] * z_p_prime_base_110[4] ] <= 0
 infeasibility_p_base_111_row_0_leq_0: 5e-06 z_p_base_111[0] + [
   0.5 p[0] * z_p_base_111[0] - 0.5 p[1] * z_p_base_111[0]
   + 0.5 p[2] * z_p_base_111[0] + 0.5 p[3] * z_p_base_111[0]
   - 0.5 p[4] * z_p_base_111[0] ] <= 0
 infeasibility_p_base_111_row_1_leq_0: 5e-06 z_p_base_111[1] + [
   0.5 p[0] * z_p_base_111[1] + 0.5 p[1] * z_p_base_111[1]
   - 0.5 p[2] * z_p_base_111[1] - 0.5 p[3] * z_p_base_111[1]
   + 0.5 p[4] * z_p_base_111[1] ] <= 0
 infeasibility_p_base_111_row_2_leq_0: 5e-06 z_p_base_111[2] + [
   p[1] * z_p_base_111[2] ] <= 0
 infeasibility_p_base_111_row_3_leq_0: 5e-06 z_p_base_111[3] + [
   - 0.5 p[0] * z_p_base_111[3] + 0.5 p[1] * z_p_base_111[3]
   + 0.5 p[2] * z_p_base_111[3] - 0.5 p[3] * z_p_base_111[3]
   + 0.5 p[4] * z_p_base_111[3] ] <= 0
 infeasibility_p_base_111_row_4_leq_0: 5e-06 z_p_base_111[4] + [
   - p[1] * z_p_base_111[4] + p[3] * z_p_base_111[4] ] <= 0
 infeasibility_p_prime_base_111_row_0_leq_0: 5e-06 z_p_prime_base_111[0]
   + [ 0.5 alpha[0] * z_p_prime_base_111[0]
   + 0.5 alpha[1] * z_p_prime_base_111[0]
   - 0.5 alpha[2] * z_p_prime_base_111[0]
   - 0.5 alpha[3] * z_p_prime_base_111[0]
   - 0.5 alpha[4] * z_p_prime_base_111[0]
   + 0.5 p[0] * z_p_prime_base_111[0] - 0.5 p[1] * z_p_prime_base_111[0]
   + 0.5 p[2] * z_p_prime_base_111[0] + 0.5 p[3] * z_p_prime_base_111[0]
   - 0.5 p[4] * z_p_prime_base_111[0] ] <= 0
 infeasibility_p_prime_base_111_row_1_leq_0: 5e-06 z_p_prime_base_111[1]
   + [ 0.5 alpha[0] * z_p_prime_base_111[1]
   - 0.5 alpha[1] * z_p_prime_base_111[1]
   + 0.5 alpha[2] * z_p_prime_base_111[1]
   + 0.5 alpha[3] * z_p_prime_base_111[1]
   + 0.5 alpha[4] * z_p_prime_base_111[1]
   + 0.5 p[0] * z_p_prime_base_111[1] + 0.5 p[1] * z_p_prime_base_111[1]
   - 0.5 p[2] * z_p_prime_base_111[1] - 0.5 p[3] * z_p_prime_base_111[1]
   + 0.5 p[4] * z_p_prime_base_111[1] ] <= 0
 infeasibility_p_prime_base_111_row_2_leq_0: 5e-06 z_p_prime_base_111[2]
   + [ - alpha[1] * z_p_prime_base_111[2] + p[1] * z_p_prime_base_111[2] ]
   <= 0
 infeasibility_p_prime_base_111_row_3_leq_0: 5e-06 z_p_prime_base_111[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_111[3]
   - 0.5 alpha[1] * z_p_prime_base_111[3]
   - 0.5 alpha[2] * z_p_prime_base_111[3]
   + 0.5 alpha[3] * z_p_prime_base_111[3]
   + 0.5 alpha[4] * z_p_prime_base_111[3]
   - 0.5 p[0] * z_p_prime_base_111[3] + 0.5 p[1] * z_p_prime_base_111[3]
   + 0.5 p[2] * z_p_prime_base_111[3] - 0.5 p[3] * z_p_prime_base_111[3]
   + 0.5 p[4] * z_p_prime_base_111[3] ] <= 0
 infeasibility_p_prime_base_111_row_4_leq_0: 5e-06 z_p_prime_base_111[4]
   + [ alpha[1] * z_p_prime_base_111[4] - alpha[3] * z_p_prime_base_111[4]
   - p[1] * z_p_prime_base_111[4] + p[3] * z_p_prime_base_111[4] ] <= 0
 infeasibility_p_base_112_row_0_leq_0: 5e-06 z_p_base_112[0] + [
   0.5 p[0] * z_p_base_112[0] + 0.5 p[1] * z_p_base_112[0]
   + 0.5 p[2] * z_p_base_112[0] - 0.5 p[3] * z_p_base_112[0]
   - 0.5 p[4] * z_p_base_112[0] ] <= 0
 infeasibility_p_base_112_row_1_leq_0: 5e-06 z_p_base_112[1] + [
   0.5 p[0] * z_p_base_112[1] - 0.5 p[1] * z_p_base_112[1]
   - 0.5 p[2] * z_p_base_112[1] + 0.5 p[3] * z_p_base_112[1]
   + 0.5 p[4] * z_p_base_112[1] ] <= 0
 infeasibility_p_base_112_row_2_leq_0: 5e-06 z_p_base_112[2] + [
   p[1] * z_p_base_112[2] ] <= 0
 infeasibility_p_base_112_row_3_leq_0: 5e-06 z_p_base_112[3] + [
   p[3] * z_p_base_112[3] ] <= 0
 infeasibility_p_base_112_row_4_leq_0: 5e-06 z_p_base_112[4] + [
   - 0.5 p[0] * z_p_base_112[4] - 0.5 p[1] * z_p_base_112[4]
   + 0.5 p[2] * z_p_base_112[4] - 0.5 p[3] * z_p_base_112[4]
   + 0.5 p[4] * z_p_base_112[4] ] <= 0
 infeasibility_p_prime_base_112_row_0_leq_0: 5e-06 z_p_prime_base_112[0]
   + [ 0.5 alpha[0] * z_p_prime_base_112[0]
   - 0.5 alpha[1] * z_p_prime_base_112[0]
   - 0.5 alpha[2] * z_p_prime_base_112[0]
   + 0.5 alpha[3] * z_p_prime_base_112[0]
   - 0.5 alpha[4] * z_p_prime_base_112[0]
   + 0.5 p[0] * z_p_prime_base_112[0] + 0.5 p[1] * z_p_prime_base_112[0]
   + 0.5 p[2] * z_p_prime_base_112[0] - 0.5 p[3] * z_p_prime_base_112[0]
   - 0.5 p[4] * z_p_prime_base_112[0] ] <= 0
 infeasibility_p_prime_base_112_row_1_leq_0: 5e-06 z_p_prime_base_112[1]
   + [ 0.5 alpha[0] * z_p_prime_base_112[1]
   + 0.5 alpha[1] * z_p_prime_base_112[1]
   + 0.5 alpha[2] * z_p_prime_base_112[1]
   - 0.5 alpha[3] * z_p_prime_base_112[1]
   + 0.5 alpha[4] * z_p_prime_base_112[1]
   + 0.5 p[0] * z_p_prime_base_112[1] - 0.5 p[1] * z_p_prime_base_112[1]
   - 0.5 p[2] * z_p_prime_base_112[1] + 0.5 p[3] * z_p_prime_base_112[1]
   + 0.5 p[4] * z_p_prime_base_112[1] ] <= 0
 infeasibility_p_prime_base_112_row_2_leq_0: 5e-06 z_p_prime_base_112[2]
   + [ - alpha[1] * z_p_prime_base_112[2] + p[1] * z_p_prime_base_112[2] ]
   <= 0
 infeasibility_p_prime_base_112_row_3_leq_0: 5e-06 z_p_prime_base_112[3]
   + [ - alpha[3] * z_p_prime_base_112[3] + p[3] * z_p_prime_base_112[3] ]
   <= 0
 infeasibility_p_prime_base_112_row_4_leq_0: 5e-06 z_p_prime_base_112[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_112[4]
   + 0.5 alpha[1] * z_p_prime_base_112[4]
   - 0.5 alpha[2] * z_p_prime_base_112[4]
   + 0.5 alpha[3] * z_p_prime_base_112[4]
   + 0.5 alpha[4] * z_p_prime_base_112[4]
   - 0.5 p[0] * z_p_prime_base_112[4] - 0.5 p[1] * z_p_prime_base_112[4]
   + 0.5 p[2] * z_p_prime_base_112[4] - 0.5 p[3] * z_p_prime_base_112[4]
   + 0.5 p[4] * z_p_prime_base_112[4] ] <= 0
 infeasibility_p_base_113_row_0_leq_0: 5e-06 z_p_base_113[0] + [
   0.5 p[0] * z_p_base_113[0] + 0.5 p[1] * z_p_base_113[0]
   + 0.5 p[2] * z_p_base_113[0] + 0.5 p[3] * z_p_base_113[0]
   - 0.5 p[4] * z_p_base_113[0] ] <= 0
 infeasibility_p_base_113_row_1_leq_0: 5e-06 z_p_base_113[1] + [
   0.5 p[0] * z_p_base_113[1] - 0.5 p[1] * z_p_base_113[1]
   - 0.5 p[2] * z_p_base_113[1] - 0.5 p[3] * z_p_base_113[1]
   + 0.5 p[4] * z_p_base_113[1] ] <= 0
 infeasibility_p_base_113_row_2_leq_0: 5e-06 z_p_base_113[2] + [
   p[1] * z_p_base_113[2] ] <= 0
 infeasibility_p_base_113_row_3_leq_0: 5e-06 z_p_base_113[3] + [
   - 0.5 p[0] * z_p_base_113[3] - 0.5 p[1] * z_p_base_113[3]
   + 0.5 p[2] * z_p_base_113[3] - 0.5 p[3] * z_p_base_113[3]
   + 0.5 p[4] * z_p_base_113[3] ] <= 0
 infeasibility_p_base_113_row_4_leq_0: 5e-06 z_p_base_113[4] + [
   p[3] * z_p_base_113[4] ] <= 0
 infeasibility_p_prime_base_113_row_0_leq_0: 5e-06 z_p_prime_base_113[0]
   + [ 0.5 alpha[0] * z_p_prime_base_113[0]
   - 0.5 alpha[1] * z_p_prime_base_113[0]
   - 0.5 alpha[2] * z_p_prime_base_113[0]
   - 0.5 alpha[3] * z_p_prime_base_113[0]
   - 0.5 alpha[4] * z_p_prime_base_113[0]
   + 0.5 p[0] * z_p_prime_base_113[0] + 0.5 p[1] * z_p_prime_base_113[0]
   + 0.5 p[2] * z_p_prime_base_113[0] + 0.5 p[3] * z_p_prime_base_113[0]
   - 0.5 p[4] * z_p_prime_base_113[0] ] <= 0
 infeasibility_p_prime_base_113_row_1_leq_0: 5e-06 z_p_prime_base_113[1]
   + [ 0.5 alpha[0] * z_p_prime_base_113[1]
   + 0.5 alpha[1] * z_p_prime_base_113[1]
   + 0.5 alpha[2] * z_p_prime_base_113[1]
   + 0.5 alpha[3] * z_p_prime_base_113[1]
   + 0.5 alpha[4] * z_p_prime_base_113[1]
   + 0.5 p[0] * z_p_prime_base_113[1] - 0.5 p[1] * z_p_prime_base_113[1]
   - 0.5 p[2] * z_p_prime_base_113[1] - 0.5 p[3] * z_p_prime_base_113[1]
   + 0.5 p[4] * z_p_prime_base_113[1] ] <= 0
 infeasibility_p_prime_base_113_row_2_leq_0: 5e-06 z_p_prime_base_113[2]
   + [ - alpha[1] * z_p_prime_base_113[2] + p[1] * z_p_prime_base_113[2] ]
   <= 0
 infeasibility_p_prime_base_113_row_3_leq_0: 5e-06 z_p_prime_base_113[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_113[3]
   + 0.5 alpha[1] * z_p_prime_base_113[3]
   - 0.5 alpha[2] * z_p_prime_base_113[3]
   + 0.5 alpha[3] * z_p_prime_base_113[3]
   + 0.5 alpha[4] * z_p_prime_base_113[3]
   - 0.5 p[0] * z_p_prime_base_113[3] - 0.5 p[1] * z_p_prime_base_113[3]
   + 0.5 p[2] * z_p_prime_base_113[3] - 0.5 p[3] * z_p_prime_base_113[3]
   + 0.5 p[4] * z_p_prime_base_113[3] ] <= 0
 infeasibility_p_prime_base_113_row_4_leq_0: 5e-06 z_p_prime_base_113[4]
   + [ - alpha[3] * z_p_prime_base_113[4] + p[3] * z_p_prime_base_113[4] ]
   <= 0
 infeasibility_p_base_114_row_0_leq_0: 5e-06 z_p_base_114[0] + [
   p[0] * z_p_base_114[0] ] <= 0
 infeasibility_p_base_114_row_1_leq_0: 5e-06 z_p_base_114[1] + [
   - 0.5 p[0] * z_p_base_114[1] + 0.5 p[1] * z_p_base_114[1]
   + 0.5 p[2] * z_p_base_114[1] - 0.5 p[3] * z_p_base_114[1]
   - 0.5 p[4] * z_p_base_114[1] ] <= 0
 infeasibility_p_base_114_row_2_leq_0: 5e-06 z_p_base_114[2] + [
   0.5 p[0] * z_p_base_114[2] + 0.5 p[1] * z_p_base_114[2]
   - 0.5 p[2] * z_p_base_114[2] + 0.5 p[3] * z_p_base_114[2]
   - 0.5 p[4] * z_p_base_114[2] ] <= 0
 infeasibility_p_base_114_row_3_leq_0: 5e-06 z_p_base_114[3] + [
   p[4] * z_p_base_114[3] ] <= 0
 infeasibility_p_base_114_row_4_leq_0: 5e-06 z_p_base_114[4] + [
   - 0.5 p[0] * z_p_base_114[4] - 0.5 p[1] * z_p_base_114[4]
   + 0.5 p[2] * z_p_base_114[4] + 0.5 p[3] * z_p_base_114[4]
   + 0.5 p[4] * z_p_base_114[4] ] <= 0
 infeasibility_p_prime_base_114_row_0_leq_0: 5e-06 z_p_prime_base_114[0]
   + [ alpha[0] * z_p_prime_base_114[0] + p[0] * z_p_prime_base_114[0] ]
   <= 0
 infeasibility_p_prime_base_114_row_1_leq_0: 5e-06 z_p_prime_base_114[1]
   + [ - 0.5 alpha[0] * z_p_prime_base_114[1]
   - 0.5 alpha[1] * z_p_prime_base_114[1]
   - 0.5 alpha[2] * z_p_prime_base_114[1]
   + 0.5 alpha[3] * z_p_prime_base_114[1]
   - 0.5 alpha[4] * z_p_prime_base_114[1]
   - 0.5 p[0] * z_p_prime_base_114[1] + 0.5 p[1] * z_p_prime_base_114[1]
   + 0.5 p[2] * z_p_prime_base_114[1] - 0.5 p[3] * z_p_prime_base_114[1]
   - 0.5 p[4] * z_p_prime_base_114[1] ] <= 0
 infeasibility_p_prime_base_114_row_2_leq_0: 5e-06 z_p_prime_base_114[2]
   + [ 0.5 alpha[0] * z_p_prime_base_114[2]
   - 0.5 alpha[1] * z_p_prime_base_114[2]
   + 0.5 alpha[2] * z_p_prime_base_114[2]
   - 0.5 alpha[3] * z_p_prime_base_114[2]
   - 0.5 alpha[4] * z_p_prime_base_114[2]
   + 0.5 p[0] * z_p_prime_base_114[2] + 0.5 p[1] * z_p_prime_base_114[2]
   - 0.5 p[2] * z_p_prime_base_114[2] + 0.5 p[3] * z_p_prime_base_114[2]
   - 0.5 p[4] * z_p_prime_base_114[2] ] <= 0
 infeasibility_p_prime_base_114_row_3_leq_0: 5e-06 z_p_prime_base_114[3]
   + [ alpha[4] * z_p_prime_base_114[3] + p[4] * z_p_prime_base_114[3] ]
   <= 0
 infeasibility_p_prime_base_114_row_4_leq_0: 5e-06 z_p_prime_base_114[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_114[4]
   + 0.5 alpha[1] * z_p_prime_base_114[4]
   - 0.5 alpha[2] * z_p_prime_base_114[4]
   - 0.5 alpha[3] * z_p_prime_base_114[4]
   + 0.5 alpha[4] * z_p_prime_base_114[4]
   - 0.5 p[0] * z_p_prime_base_114[4] - 0.5 p[1] * z_p_prime_base_114[4]
   + 0.5 p[2] * z_p_prime_base_114[4] + 0.5 p[3] * z_p_prime_base_114[4]
   + 0.5 p[4] * z_p_prime_base_114[4] ] <= 0
 infeasibility_p_base_115_row_0_leq_0: 5e-06 z_p_base_115[0] + [
   p[0] * z_p_base_115[0] ] <= 0
 infeasibility_p_base_115_row_1_leq_0: 5e-06 z_p_base_115[1] + [
   - 0.5 p[0] * z_p_base_115[1] + 0.5 p[1] * z_p_base_115[1]
   + 0.5 p[2] * z_p_base_115[1] - 0.5 p[3] * z_p_base_115[1]
   - 0.5 p[4] * z_p_base_115[1] ] <= 0
 infeasibility_p_base_115_row_2_leq_0: 5e-06 z_p_base_115[2] + [
   p[3] * z_p_base_115[2] ] <= 0
 infeasibility_p_base_115_row_3_leq_0: 5e-06 z_p_base_115[3] + [
   0.5 p[0] * z_p_base_115[3] + 0.5 p[1] * z_p_base_115[3]
   - 0.5 p[2] * z_p_base_115[3] - 0.5 p[3] * z_p_base_115[3]
   + 0.5 p[4] * z_p_base_115[3] ] <= 0
 infeasibility_p_base_115_row_4_leq_0: 5e-06 z_p_base_115[4] + [
   - 0.5 p[0] * z_p_base_115[4] - 0.5 p[1] * z_p_base_115[4]
   + 0.5 p[2] * z_p_base_115[4] + 0.5 p[3] * z_p_base_115[4]
   + 0.5 p[4] * z_p_base_115[4] ] <= 0
 infeasibility_p_prime_base_115_row_0_leq_0: 5e-06 z_p_prime_base_115[0]
   + [ alpha[0] * z_p_prime_base_115[0] + p[0] * z_p_prime_base_115[0] ]
   <= 0
 infeasibility_p_prime_base_115_row_1_leq_0: 5e-06 z_p_prime_base_115[1]
   + [ - 0.5 alpha[0] * z_p_prime_base_115[1]
   - 0.5 alpha[1] * z_p_prime_base_115[1]
   - 0.5 alpha[2] * z_p_prime_base_115[1]
   + 0.5 alpha[3] * z_p_prime_base_115[1]
   - 0.5 alpha[4] * z_p_prime_base_115[1]
   - 0.5 p[0] * z_p_prime_base_115[1] + 0.5 p[1] * z_p_prime_base_115[1]
   + 0.5 p[2] * z_p_prime_base_115[1] - 0.5 p[3] * z_p_prime_base_115[1]
   - 0.5 p[4] * z_p_prime_base_115[1] ] <= 0
 infeasibility_p_prime_base_115_row_2_leq_0: 5e-06 z_p_prime_base_115[2]
   + [ - alpha[3] * z_p_prime_base_115[2] + p[3] * z_p_prime_base_115[2] ]
   <= 0
 infeasibility_p_prime_base_115_row_3_leq_0: 5e-06 z_p_prime_base_115[3]
   + [ 0.5 alpha[0] * z_p_prime_base_115[3]
   - 0.5 alpha[1] * z_p_prime_base_115[3]
   + 0.5 alpha[2] * z_p_prime_base_115[3]
   + 0.5 alpha[3] * z_p_prime_base_115[3]
   + 0.5 alpha[4] * z_p_prime_base_115[3]
   + 0.5 p[0] * z_p_prime_base_115[3] + 0.5 p[1] * z_p_prime_base_115[3]
   - 0.5 p[2] * z_p_prime_base_115[3] - 0.5 p[3] * z_p_prime_base_115[3]
   + 0.5 p[4] * z_p_prime_base_115[3] ] <= 0
 infeasibility_p_prime_base_115_row_4_leq_0: 5e-06 z_p_prime_base_115[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_115[4]
   + 0.5 alpha[1] * z_p_prime_base_115[4]
   - 0.5 alpha[2] * z_p_prime_base_115[4]
   - 0.5 alpha[3] * z_p_prime_base_115[4]
   + 0.5 alpha[4] * z_p_prime_base_115[4]
   - 0.5 p[0] * z_p_prime_base_115[4] - 0.5 p[1] * z_p_prime_base_115[4]
   + 0.5 p[2] * z_p_prime_base_115[4] + 0.5 p[3] * z_p_prime_base_115[4]
   + 0.5 p[4] * z_p_prime_base_115[4] ] <= 0
 infeasibility_p_base_116_row_0_leq_0: 5e-06 z_p_base_116[0] + [
   p[0] * z_p_base_116[0] ] <= 0
 infeasibility_p_base_116_row_1_leq_0: 5e-06 z_p_base_116[1] + [
   - p[0] * z_p_base_116[1] + p[2] * z_p_base_116[1] ] <= 0
 infeasibility_p_base_116_row_2_leq_0: 5e-06 z_p_base_116[2] + [
   0.5 p[0] * z_p_base_116[2] + 0.5 p[1] * z_p_base_116[2]
   - 0.5 p[2] * z_p_base_116[2] + 0.5 p[3] * z_p_base_116[2]
   - 0.5 p[4] * z_p_base_116[2] ] <= 0
 infeasibility_p_base_116_row_3_leq_0: 5e-06 z_p_base_116[3] + [
   0.5 p[0] * z_p_base_116[3] + 0.5 p[1] * z_p_base_116[3]
   - 0.5 p[2] * z_p_base_116[3] - 0.5 p[3] * z_p_base_116[3]
   + 0.5 p[4] * z_p_base_116[3] ] <= 0
 infeasibility_p_base_116_row_4_leq_0: 5e-06 z_p_base_116[4] + [
   - 0.5 p[0] * z_p_base_116[4] - 0.5 p[1] * z_p_base_116[4]
   + 0.5 p[2] * z_p_base_116[4] + 0.5 p[3] * z_p_base_116[4]
   + 0.5 p[4] * z_p_base_116[4] ] <= 0
 infeasibility_p_prime_base_116_row_0_leq_0: 5e-06 z_p_prime_base_116[0]
   + [ alpha[0] * z_p_prime_base_116[0] + p[0] * z_p_prime_base_116[0] ]
   <= 0
 infeasibility_p_prime_base_116_row_1_leq_0: 5e-06 z_p_prime_base_116[1]
   + [ - alpha[0] * z_p_prime_base_116[1]
   - alpha[2] * z_p_prime_base_116[1] - p[0] * z_p_prime_base_116[1]
   + p[2] * z_p_prime_base_116[1] ] <= 0
 infeasibility_p_prime_base_116_row_2_leq_0: 5e-06 z_p_prime_base_116[2]
   + [ 0.5 alpha[0] * z_p_prime_base_116[2]
   - 0.5 alpha[1] * z_p_prime_base_116[2]
   + 0.5 alpha[2] * z_p_prime_base_116[2]
   - 0.5 alpha[3] * z_p_prime_base_116[2]
   - 0.5 alpha[4] * z_p_prime_base_116[2]
   + 0.5 p[0] * z_p_prime_base_116[2] + 0.5 p[1] * z_p_prime_base_116[2]
   - 0.5 p[2] * z_p_prime_base_116[2] + 0.5 p[3] * z_p_prime_base_116[2]
   - 0.5 p[4] * z_p_prime_base_116[2] ] <= 0
 infeasibility_p_prime_base_116_row_3_leq_0: 5e-06 z_p_prime_base_116[3]
   + [ 0.5 alpha[0] * z_p_prime_base_116[3]
   - 0.5 alpha[1] * z_p_prime_base_116[3]
   + 0.5 alpha[2] * z_p_prime_base_116[3]
   + 0.5 alpha[3] * z_p_prime_base_116[3]
   + 0.5 alpha[4] * z_p_prime_base_116[3]
   + 0.5 p[0] * z_p_prime_base_116[3] + 0.5 p[1] * z_p_prime_base_116[3]
   - 0.5 p[2] * z_p_prime_base_116[3] - 0.5 p[3] * z_p_prime_base_116[3]
   + 0.5 p[4] * z_p_prime_base_116[3] ] <= 0
 infeasibility_p_prime_base_116_row_4_leq_0: 5e-06 z_p_prime_base_116[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_116[4]
   + 0.5 alpha[1] * z_p_prime_base_116[4]
   - 0.5 alpha[2] * z_p_prime_base_116[4]
   - 0.5 alpha[3] * z_p_prime_base_116[4]
   + 0.5 alpha[4] * z_p_prime_base_116[4]
   - 0.5 p[0] * z_p_prime_base_116[4] - 0.5 p[1] * z_p_prime_base_116[4]
   + 0.5 p[2] * z_p_prime_base_116[4] + 0.5 p[3] * z_p_prime_base_116[4]
   + 0.5 p[4] * z_p_prime_base_116[4] ] <= 0
 infeasibility_p_base_117_row_0_leq_0: 5e-06 z_p_base_117[0] + [
   p[0] * z_p_base_117[0] ] <= 0
 infeasibility_p_base_117_row_1_leq_0: 5e-06 z_p_base_117[1] + [
   - 0.5 p[0] * z_p_base_117[1] + 0.5 p[1] * z_p_base_117[1]
   + 0.5 p[2] * z_p_base_117[1] - 0.5 p[3] * z_p_base_117[1]
   - 0.5 p[4] * z_p_base_117[1] ] <= 0
 infeasibility_p_base_117_row_2_leq_0: 5e-06 z_p_base_117[2] + [
   0.5 p[0] * z_p_base_117[2] + 0.5 p[1] * z_p_base_117[2]
   - 0.5 p[2] * z_p_base_117[2] + 0.5 p[3] * z_p_base_117[2]
   + 0.5 p[4] * z_p_base_117[2] ] <= 0
 infeasibility_p_base_117_row_3_leq_0: 5e-06 z_p_base_117[3] + [
   - 0.5 p[0] * z_p_base_117[3] - 0.5 p[1] * z_p_base_117[3]
   + 0.5 p[2] * z_p_base_117[3] + 0.5 p[3] * z_p_base_117[3]
   - 0.5 p[4] * z_p_base_117[3] ] <= 0
 infeasibility_p_base_117_row_4_leq_0: 5e-06 z_p_base_117[4] + [
   p[4] * z_p_base_117[4] ] <= 0
 infeasibility_p_prime_base_117_row_0_leq_0: 5e-06 z_p_prime_base_117[0]
   + [ alpha[0] * z_p_prime_base_117[0] + p[0] * z_p_prime_base_117[0] ]
   <= 0
 infeasibility_p_prime_base_117_row_1_leq_0: 5e-06 z_p_prime_base_117[1]
   + [ - 0.5 alpha[0] * z_p_prime_base_117[1]
   - 0.5 alpha[1] * z_p_prime_base_117[1]
   - 0.5 alpha[2] * z_p_prime_base_117[1]
   + 0.5 alpha[3] * z_p_prime_base_117[1]
   - 0.5 alpha[4] * z_p_prime_base_117[1]
   - 0.5 p[0] * z_p_prime_base_117[1] + 0.5 p[1] * z_p_prime_base_117[1]
   + 0.5 p[2] * z_p_prime_base_117[1] - 0.5 p[3] * z_p_prime_base_117[1]
   - 0.5 p[4] * z_p_prime_base_117[1] ] <= 0
 infeasibility_p_prime_base_117_row_2_leq_0: 5e-06 z_p_prime_base_117[2]
   + [ 0.5 alpha[0] * z_p_prime_base_117[2]
   - 0.5 alpha[1] * z_p_prime_base_117[2]
   + 0.5 alpha[2] * z_p_prime_base_117[2]
   - 0.5 alpha[3] * z_p_prime_base_117[2]
   + 0.5 alpha[4] * z_p_prime_base_117[2]
   + 0.5 p[0] * z_p_prime_base_117[2] + 0.5 p[1] * z_p_prime_base_117[2]
   - 0.5 p[2] * z_p_prime_base_117[2] + 0.5 p[3] * z_p_prime_base_117[2]
   + 0.5 p[4] * z_p_prime_base_117[2] ] <= 0
 infeasibility_p_prime_base_117_row_3_leq_0: 5e-06 z_p_prime_base_117[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_117[3]
   + 0.5 alpha[1] * z_p_prime_base_117[3]
   - 0.5 alpha[2] * z_p_prime_base_117[3]
   - 0.5 alpha[3] * z_p_prime_base_117[3]
   - 0.5 alpha[4] * z_p_prime_base_117[3]
   - 0.5 p[0] * z_p_prime_base_117[3] - 0.5 p[1] * z_p_prime_base_117[3]
   + 0.5 p[2] * z_p_prime_base_117[3] + 0.5 p[3] * z_p_prime_base_117[3]
   - 0.5 p[4] * z_p_prime_base_117[3] ] <= 0
 infeasibility_p_prime_base_117_row_4_leq_0: 5e-06 z_p_prime_base_117[4]
   + [ alpha[4] * z_p_prime_base_117[4] + p[4] * z_p_prime_base_117[4] ]
   <= 0
 infeasibility_p_base_118_row_0_leq_0: 5e-06 z_p_base_118[0] + [
   p[0] * z_p_base_118[0] ] <= 0
 infeasibility_p_base_118_row_1_leq_0: 5e-06 z_p_base_118[1] + [
   - 0.5 p[0] * z_p_base_118[1] + 0.5 p[1] * z_p_base_118[1]
   + 0.5 p[2] * z_p_base_118[1] - 0.5 p[3] * z_p_base_118[1]
   + 0.5 p[4] * z_p_base_118[1] ] <= 0
 infeasibility_p_base_118_row_2_leq_0: 5e-06 z_p_base_118[2] + [
   0.5 p[0] * z_p_base_118[2] + 0.5 p[1] * z_p_base_118[2]
   - 0.5 p[2] * z_p_base_118[2] + 0.5 p[3] * z_p_base_118[2]
   - 0.5 p[4] * z_p_base_118[2] ] <= 0
 infeasibility_p_base_118_row_3_leq_0: 5e-06 z_p_base_118[3] + [
   - 0.5 p[0] * z_p_base_118[3] - 0.5 p[1] * z_p_base_118[3]
   + 0.5 p[2] * z_p_base_118[3] + 0.5 p[3] * z_p_base_118[3]
   - 0.5 p[4] * z_p_base_118[3] ] <= 0
 infeasibility_p_base_118_row_4_leq_0: 5e-06 z_p_base_118[4] + [
   p[4] * z_p_base_118[4] ] <= 0
 infeasibility_p_prime_base_118_row_0_leq_0: 5e-06 z_p_prime_base_118[0]
   + [ alpha[0] * z_p_prime_base_118[0] + p[0] * z_p_prime_base_118[0] ]
   <= 0
 infeasibility_p_prime_base_118_row_1_leq_0: 5e-06 z_p_prime_base_118[1]
   + [ - 0.5 alpha[0] * z_p_prime_base_118[1]
   - 0.5 alpha[1] * z_p_prime_base_118[1]
   - 0.5 alpha[2] * z_p_prime_base_118[1]
   + 0.5 alpha[3] * z_p_prime_base_118[1]
   + 0.5 alpha[4] * z_p_prime_base_118[1]
   - 0.5 p[0] * z_p_prime_base_118[1] + 0.5 p[1] * z_p_prime_base_118[1]
   + 0.5 p[2] * z_p_prime_base_118[1] - 0.5 p[3] * z_p_prime_base_118[1]
   + 0.5 p[4] * z_p_prime_base_118[1] ] <= 0
 infeasibility_p_prime_base_118_row_2_leq_0: 5e-06 z_p_prime_base_118[2]
   + [ 0.5 alpha[0] * z_p_prime_base_118[2]
   - 0.5 alpha[1] * z_p_prime_base_118[2]
   + 0.5 alpha[2] * z_p_prime_base_118[2]
   - 0.5 alpha[3] * z_p_prime_base_118[2]
   - 0.5 alpha[4] * z_p_prime_base_118[2]
   + 0.5 p[0] * z_p_prime_base_118[2] + 0.5 p[1] * z_p_prime_base_118[2]
   - 0.5 p[2] * z_p_prime_base_118[2] + 0.5 p[3] * z_p_prime_base_118[2]
   - 0.5 p[4] * z_p_prime_base_118[2] ] <= 0
 infeasibility_p_prime_base_118_row_3_leq_0: 5e-06 z_p_prime_base_118[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_118[3]
   + 0.5 alpha[1] * z_p_prime_base_118[3]
   - 0.5 alpha[2] * z_p_prime_base_118[3]
   - 0.5 alpha[3] * z_p_prime_base_118[3]
   - 0.5 alpha[4] * z_p_prime_base_118[3]
   - 0.5 p[0] * z_p_prime_base_118[3] - 0.5 p[1] * z_p_prime_base_118[3]
   + 0.5 p[2] * z_p_prime_base_118[3] + 0.5 p[3] * z_p_prime_base_118[3]
   - 0.5 p[4] * z_p_prime_base_118[3] ] <= 0
 infeasibility_p_prime_base_118_row_4_leq_0: 5e-06 z_p_prime_base_118[4]
   + [ alpha[4] * z_p_prime_base_118[4] + p[4] * z_p_prime_base_118[4] ]
   <= 0
 infeasibility_p_base_119_row_0_leq_0: 5e-06 z_p_base_119[0] + [
   p[0] * z_p_base_119[0] ] <= 0
 infeasibility_p_base_119_row_1_leq_0: 5e-06 z_p_base_119[1] + [
   - 0.5 p[0] * z_p_base_119[1] + 0.5 p[1] * z_p_base_119[1]
   + 0.5 p[2] * z_p_base_119[1] - 0.5 p[3] * z_p_base_119[1]
   - 0.5 p[4] * z_p_base_119[1] ] <= 0
 infeasibility_p_base_119_row_2_leq_0: 5e-06 z_p_base_119[2] + [
   0.5 p[0] * z_p_base_119[2] + 0.5 p[1] * z_p_base_119[2]
   - 0.5 p[2] * z_p_base_119[2] + 0.5 p[3] * z_p_base_119[2]
   + 0.5 p[4] * z_p_base_119[2] ] <= 0
 infeasibility_p_base_119_row_3_leq_0: 5e-06 z_p_base_119[3] + [
   p[3] * z_p_base_119[3] ] <= 0
 infeasibility_p_base_119_row_4_leq_0: 5e-06 z_p_base_119[4] + [
   - 0.5 p[0] * z_p_base_119[4] - 0.5 p[1] * z_p_base_119[4]
   + 0.5 p[2] * z_p_base_119[4] - 0.5 p[3] * z_p_base_119[4]
   + 0.5 p[4] * z_p_base_119[4] ] <= 0
 infeasibility_p_prime_base_119_row_0_leq_0: 5e-06 z_p_prime_base_119[0]
   + [ alpha[0] * z_p_prime_base_119[0] + p[0] * z_p_prime_base_119[0] ]
   <= 0
 infeasibility_p_prime_base_119_row_1_leq_0: 5e-06 z_p_prime_base_119[1]
   + [ - 0.5 alpha[0] * z_p_prime_base_119[1]
   - 0.5 alpha[1] * z_p_prime_base_119[1]
   - 0.5 alpha[2] * z_p_prime_base_119[1]
   + 0.5 alpha[3] * z_p_prime_base_119[1]
   - 0.5 alpha[4] * z_p_prime_base_119[1]
   - 0.5 p[0] * z_p_prime_base_119[1] + 0.5 p[1] * z_p_prime_base_119[1]
   + 0.5 p[2] * z_p_prime_base_119[1] - 0.5 p[3] * z_p_prime_base_119[1]
   - 0.5 p[4] * z_p_prime_base_119[1] ] <= 0
 infeasibility_p_prime_base_119_row_2_leq_0: 5e-06 z_p_prime_base_119[2]
   + [ 0.5 alpha[0] * z_p_prime_base_119[2]
   - 0.5 alpha[1] * z_p_prime_base_119[2]
   + 0.5 alpha[2] * z_p_prime_base_119[2]
   - 0.5 alpha[3] * z_p_prime_base_119[2]
   + 0.5 alpha[4] * z_p_prime_base_119[2]
   + 0.5 p[0] * z_p_prime_base_119[2] + 0.5 p[1] * z_p_prime_base_119[2]
   - 0.5 p[2] * z_p_prime_base_119[2] + 0.5 p[3] * z_p_prime_base_119[2]
   + 0.5 p[4] * z_p_prime_base_119[2] ] <= 0
 infeasibility_p_prime_base_119_row_3_leq_0: 5e-06 z_p_prime_base_119[3]
   + [ - alpha[3] * z_p_prime_base_119[3] + p[3] * z_p_prime_base_119[3] ]
   <= 0
 infeasibility_p_prime_base_119_row_4_leq_0: 5e-06 z_p_prime_base_119[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_119[4]
   + 0.5 alpha[1] * z_p_prime_base_119[4]
   - 0.5 alpha[2] * z_p_prime_base_119[4]
   + 0.5 alpha[3] * z_p_prime_base_119[4]
   + 0.5 alpha[4] * z_p_prime_base_119[4]
   - 0.5 p[0] * z_p_prime_base_119[4] - 0.5 p[1] * z_p_prime_base_119[4]
   + 0.5 p[2] * z_p_prime_base_119[4] - 0.5 p[3] * z_p_prime_base_119[4]
   + 0.5 p[4] * z_p_prime_base_119[4] ] <= 0
 infeasibility_p_base_120_row_0_leq_0: 5e-06 z_p_base_120[0] + [
   p[0] * z_p_base_120[0] ] <= 0
 infeasibility_p_base_120_row_1_leq_0: 5e-06 z_p_base_120[1] + [
   - 0.5 p[0] * z_p_base_120[1] + 0.5 p[1] * z_p_base_120[1]
   + 0.5 p[2] * z_p_base_120[1] + 0.5 p[3] * z_p_base_120[1]
   - 0.5 p[4] * z_p_base_120[1] ] <= 0
 infeasibility_p_base_120_row_2_leq_0: 5e-06 z_p_base_120[2] + [
   0.5 p[0] * z_p_base_120[2] + 0.5 p[1] * z_p_base_120[2]
   - 0.5 p[2] * z_p_base_120[2] - 0.5 p[3] * z_p_base_120[2]
   + 0.5 p[4] * z_p_base_120[2] ] <= 0
 infeasibility_p_base_120_row_3_leq_0: 5e-06 z_p_base_120[3] + [
   - 0.5 p[0] * z_p_base_120[3] - 0.5 p[1] * z_p_base_120[3]
   + 0.5 p[2] * z_p_base_120[3] - 0.5 p[3] * z_p_base_120[3]
   + 0.5 p[4] * z_p_base_120[3] ] <= 0
 infeasibility_p_base_120_row_4_leq_0: 5e-06 z_p_base_120[4] + [
   p[3] * z_p_base_120[4] ] <= 0
 infeasibility_p_prime_base_120_row_0_leq_0: 5e-06 z_p_prime_base_120[0]
   + [ alpha[0] * z_p_prime_base_120[0] + p[0] * z_p_prime_base_120[0] ]
   <= 0
 infeasibility_p_prime_base_120_row_1_leq_0: 5e-06 z_p_prime_base_120[1]
   + [ - 0.5 alpha[0] * z_p_prime_base_120[1]
   - 0.5 alpha[1] * z_p_prime_base_120[1]
   - 0.5 alpha[2] * z_p_prime_base_120[1]
   - 0.5 alpha[3] * z_p_prime_base_120[1]
   - 0.5 alpha[4] * z_p_prime_base_120[1]
   - 0.5 p[0] * z_p_prime_base_120[1] + 0.5 p[1] * z_p_prime_base_120[1]
   + 0.5 p[2] * z_p_prime_base_120[1] + 0.5 p[3] * z_p_prime_base_120[1]
   - 0.5 p[4] * z_p_prime_base_120[1] ] <= 0
 infeasibility_p_prime_base_120_row_2_leq_0: 5e-06 z_p_prime_base_120[2]
   + [ 0.5 alpha[0] * z_p_prime_base_120[2]
   - 0.5 alpha[1] * z_p_prime_base_120[2]
   + 0.5 alpha[2] * z_p_prime_base_120[2]
   + 0.5 alpha[3] * z_p_prime_base_120[2]
   + 0.5 alpha[4] * z_p_prime_base_120[2]
   + 0.5 p[0] * z_p_prime_base_120[2] + 0.5 p[1] * z_p_prime_base_120[2]
   - 0.5 p[2] * z_p_prime_base_120[2] - 0.5 p[3] * z_p_prime_base_120[2]
   + 0.5 p[4] * z_p_prime_base_120[2] ] <= 0
 infeasibility_p_prime_base_120_row_3_leq_0: 5e-06 z_p_prime_base_120[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_120[3]
   + 0.5 alpha[1] * z_p_prime_base_120[3]
   - 0.5 alpha[2] * z_p_prime_base_120[3]
   + 0.5 alpha[3] * z_p_prime_base_120[3]
   + 0.5 alpha[4] * z_p_prime_base_120[3]
   - 0.5 p[0] * z_p_prime_base_120[3] - 0.5 p[1] * z_p_prime_base_120[3]
   + 0.5 p[2] * z_p_prime_base_120[3] - 0.5 p[3] * z_p_prime_base_120[3]
   + 0.5 p[4] * z_p_prime_base_120[3] ] <= 0
 infeasibility_p_prime_base_120_row_4_leq_0: 5e-06 z_p_prime_base_120[4]
   + [ - alpha[3] * z_p_prime_base_120[4] + p[3] * z_p_prime_base_120[4] ]
   <= 0
 infeasibility_p_base_121_row_0_leq_0: 5e-06 z_p_base_121[0] + [
   p[0] * z_p_base_121[0] ] <= 0
 infeasibility_p_base_121_row_1_leq_0: 5e-06 z_p_base_121[1] + [
   p[1] * z_p_base_121[1] ] <= 0
 infeasibility_p_base_121_row_2_leq_0: 5e-06 z_p_base_121[2] + [
   - 0.5 p[0] * z_p_base_121[2] - 0.5 p[1] * z_p_base_121[2]
   + 0.5 p[2] * z_p_base_121[2] + 0.5 p[3] * z_p_base_121[2]
   - 0.5 p[4] * z_p_base_121[2] ] <= 0
 infeasibility_p_base_121_row_3_leq_0: 5e-06 z_p_base_121[3] + [
   - 0.5 p[0] * z_p_base_121[3] - 0.5 p[1] * z_p_base_121[3]
   + 0.5 p[2] * z_p_base_121[3] - 0.5 p[3] * z_p_base_121[3]
   + 0.5 p[4] * z_p_base_121[3] ] <= 0
 infeasibility_p_base_121_row_4_leq_0: 5e-06 z_p_base_121[4] + [
   0.5 p[0] * z_p_base_121[4] + 0.5 p[1] * z_p_base_121[4]
   - 0.5 p[2] * z_p_base_121[4] + 0.5 p[3] * z_p_base_121[4]
   + 0.5 p[4] * z_p_base_121[4] ] <= 0
 infeasibility_p_prime_base_121_row_0_leq_0: 5e-06 z_p_prime_base_121[0]
   + [ alpha[0] * z_p_prime_base_121[0] + p[0] * z_p_prime_base_121[0] ]
   <= 0
 infeasibility_p_prime_base_121_row_1_leq_0: 5e-06 z_p_prime_base_121[1]
   + [ - alpha[1] * z_p_prime_base_121[1] + p[1] * z_p_prime_base_121[1] ]
   <= 0
 infeasibility_p_prime_base_121_row_2_leq_0: 5e-06 z_p_prime_base_121[2]
   + [ - 0.5 alpha[0] * z_p_prime_base_121[2]
   + 0.5 alpha[1] * z_p_prime_base_121[2]
   - 0.5 alpha[2] * z_p_prime_base_121[2]
   - 0.5 alpha[3] * z_p_prime_base_121[2]
   - 0.5 alpha[4] * z_p_prime_base_121[2]
   - 0.5 p[0] * z_p_prime_base_121[2] - 0.5 p[1] * z_p_prime_base_121[2]
   + 0.5 p[2] * z_p_prime_base_121[2] + 0.5 p[3] * z_p_prime_base_121[2]
   - 0.5 p[4] * z_p_prime_base_121[2] ] <= 0
 infeasibility_p_prime_base_121_row_3_leq_0: 5e-06 z_p_prime_base_121[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_121[3]
   + 0.5 alpha[1] * z_p_prime_base_121[3]
   - 0.5 alpha[2] * z_p_prime_base_121[3]
   + 0.5 alpha[3] * z_p_prime_base_121[3]
   + 0.5 alpha[4] * z_p_prime_base_121[3]
   - 0.5 p[0] * z_p_prime_base_121[3] - 0.5 p[1] * z_p_prime_base_121[3]
   + 0.5 p[2] * z_p_prime_base_121[3] - 0.5 p[3] * z_p_prime_base_121[3]
   + 0.5 p[4] * z_p_prime_base_121[3] ] <= 0
 infeasibility_p_prime_base_121_row_4_leq_0: 5e-06 z_p_prime_base_121[4]
   + [ 0.5 alpha[0] * z_p_prime_base_121[4]
   - 0.5 alpha[1] * z_p_prime_base_121[4]
   + 0.5 alpha[2] * z_p_prime_base_121[4]
   - 0.5 alpha[3] * z_p_prime_base_121[4]
   + 0.5 alpha[4] * z_p_prime_base_121[4]
   + 0.5 p[0] * z_p_prime_base_121[4] + 0.5 p[1] * z_p_prime_base_121[4]
   - 0.5 p[2] * z_p_prime_base_121[4] + 0.5 p[3] * z_p_prime_base_121[4]
   + 0.5 p[4] * z_p_prime_base_121[4] ] <= 0
 infeasibility_p_base_122_row_0_leq_0: 5e-06 z_p_base_122[0] + [
   p[0] * z_p_base_122[0] ] <= 0
 infeasibility_p_base_122_row_1_leq_0: 5e-06 z_p_base_122[1] + [
   0.5 p[0] * z_p_base_122[1] + 0.5 p[1] * z_p_base_122[1]
   - 0.5 p[2] * z_p_base_122[1] + 0.5 p[3] * z_p_base_122[1]
   - 0.5 p[4] * z_p_base_122[1] ] <= 0
 infeasibility_p_base_122_row_2_leq_0: 5e-06 z_p_base_122[2] + [
   - 0.5 p[0] * z_p_base_122[2] + 0.5 p[1] * z_p_base_122[2]
   + 0.5 p[2] * z_p_base_122[2] - 0.5 p[3] * z_p_base_122[2]
   + 0.5 p[4] * z_p_base_122[2] ] <= 0
 infeasibility_p_base_122_row_3_leq_0: 5e-06 z_p_base_122[3] + [
   - p[0] * z_p_base_122[3] + p[2] * z_p_base_122[3] ] <= 0
 infeasibility_p_base_122_row_4_leq_0: 5e-06 z_p_base_122[4] + [
   0.5 p[0] * z_p_base_122[4] - 0.5 p[1] * z_p_base_122[4]
   - 0.5 p[2] * z_p_base_122[4] + 0.5 p[3] * z_p_base_122[4]
   + 0.5 p[4] * z_p_base_122[4] ] <= 0
 infeasibility_p_prime_base_122_row_0_leq_0: 5e-06 z_p_prime_base_122[0]
   + [ alpha[0] * z_p_prime_base_122[0] + p[0] * z_p_prime_base_122[0] ]
   <= 0
 infeasibility_p_prime_base_122_row_1_leq_0: 5e-06 z_p_prime_base_122[1]
   + [ 0.5 alpha[0] * z_p_prime_base_122[1]
   - 0.5 alpha[1] * z_p_prime_base_122[1]
   + 0.5 alpha[2] * z_p_prime_base_122[1]
   - 0.5 alpha[3] * z_p_prime_base_122[1]
   - 0.5 alpha[4] * z_p_prime_base_122[1]
   + 0.5 p[0] * z_p_prime_base_122[1] + 0.5 p[1] * z_p_prime_base_122[1]
   - 0.5 p[2] * z_p_prime_base_122[1] + 0.5 p[3] * z_p_prime_base_122[1]
   - 0.5 p[4] * z_p_prime_base_122[1] ] <= 0
 infeasibility_p_prime_base_122_row_2_leq_0: 5e-06 z_p_prime_base_122[2]
   + [ - 0.5 alpha[0] * z_p_prime_base_122[2]
   - 0.5 alpha[1] * z_p_prime_base_122[2]
   - 0.5 alpha[2] * z_p_prime_base_122[2]
   + 0.5 alpha[3] * z_p_prime_base_122[2]
   + 0.5 alpha[4] * z_p_prime_base_122[2]
   - 0.5 p[0] * z_p_prime_base_122[2] + 0.5 p[1] * z_p_prime_base_122[2]
   + 0.5 p[2] * z_p_prime_base_122[2] - 0.5 p[3] * z_p_prime_base_122[2]
   + 0.5 p[4] * z_p_prime_base_122[2] ] <= 0
 infeasibility_p_prime_base_122_row_3_leq_0: 5e-06 z_p_prime_base_122[3]
   + [ - alpha[0] * z_p_prime_base_122[3]
   - alpha[2] * z_p_prime_base_122[3] - p[0] * z_p_prime_base_122[3]
   + p[2] * z_p_prime_base_122[3] ] <= 0
 infeasibility_p_prime_base_122_row_4_leq_0: 5e-06 z_p_prime_base_122[4]
   + [ 0.5 alpha[0] * z_p_prime_base_122[4]
   + 0.5 alpha[1] * z_p_prime_base_122[4]
   + 0.5 alpha[2] * z_p_prime_base_122[4]
   - 0.5 alpha[3] * z_p_prime_base_122[4]
   + 0.5 alpha[4] * z_p_prime_base_122[4]
   + 0.5 p[0] * z_p_prime_base_122[4] - 0.5 p[1] * z_p_prime_base_122[4]
   - 0.5 p[2] * z_p_prime_base_122[4] + 0.5 p[3] * z_p_prime_base_122[4]
   + 0.5 p[4] * z_p_prime_base_122[4] ] <= 0
 infeasibility_p_base_123_row_0_leq_0: 5e-06 z_p_base_123[0] + [
   p[0] * z_p_base_123[0] ] <= 0
 infeasibility_p_base_123_row_1_leq_0: 5e-06 z_p_base_123[1] + [
   - 0.5 p[0] * z_p_base_123[1] + 0.5 p[1] * z_p_base_123[1]
   + 0.5 p[2] * z_p_base_123[1] + 0.5 p[3] * z_p_base_123[1]
   - 0.5 p[4] * z_p_base_123[1] ] <= 0
 infeasibility_p_base_123_row_2_leq_0: 5e-06 z_p_base_123[2] + [
   0.5 p[0] * z_p_base_123[2] + 0.5 p[1] * z_p_base_123[2]
   - 0.5 p[2] * z_p_base_123[2] - 0.5 p[3] * z_p_base_123[2]
   + 0.5 p[4] * z_p_base_123[2] ] <= 0
 infeasibility_p_base_123_row_3_leq_0: 5e-06 z_p_base_123[3] + [
   - p[0] * z_p_base_123[3] + p[2] * z_p_base_123[3] ] <= 0
 infeasibility_p_base_123_row_4_leq_0: 5e-06 z_p_base_123[4] + [
   0.5 p[0] * z_p_base_123[4] - 0.5 p[1] * z_p_base_123[4]
   - 0.5 p[2] * z_p_base_123[4] + 0.5 p[3] * z_p_base_123[4]
   + 0.5 p[4] * z_p_base_123[4] ] <= 0
 infeasibility_p_prime_base_123_row_0_leq_0: 5e-06 z_p_prime_base_123[0]
   + [ alpha[0] * z_p_prime_base_123[0] + p[0] * z_p_prime_base_123[0] ]
   <= 0
 infeasibility_p_prime_base_123_row_1_leq_0: 5e-06 z_p_prime_base_123[1]
   + [ - 0.5 alpha[0] * z_p_prime_base_123[1]
   - 0.5 alpha[1] * z_p_prime_base_123[1]
   - 0.5 alpha[2] * z_p_prime_base_123[1]
   - 0.5 alpha[3] * z_p_prime_base_123[1]
   - 0.5 alpha[4] * z_p_prime_base_123[1]
   - 0.5 p[0] * z_p_prime_base_123[1] + 0.5 p[1] * z_p_prime_base_123[1]
   + 0.5 p[2] * z_p_prime_base_123[1] + 0.5 p[3] * z_p_prime_base_123[1]
   - 0.5 p[4] * z_p_prime_base_123[1] ] <= 0
 infeasibility_p_prime_base_123_row_2_leq_0: 5e-06 z_p_prime_base_123[2]
   + [ 0.5 alpha[0] * z_p_prime_base_123[2]
   - 0.5 alpha[1] * z_p_prime_base_123[2]
   + 0.5 alpha[2] * z_p_prime_base_123[2]
   + 0.5 alpha[3] * z_p_prime_base_123[2]
   + 0.5 alpha[4] * z_p_prime_base_123[2]
   + 0.5 p[0] * z_p_prime_base_123[2] + 0.5 p[1] * z_p_prime_base_123[2]
   - 0.5 p[2] * z_p_prime_base_123[2] - 0.5 p[3] * z_p_prime_base_123[2]
   + 0.5 p[4] * z_p_prime_base_123[2] ] <= 0
 infeasibility_p_prime_base_123_row_3_leq_0: 5e-06 z_p_prime_base_123[3]
   + [ - alpha[0] * z_p_prime_base_123[3]
   - alpha[2] * z_p_prime_base_123[3] - p[0] * z_p_prime_base_123[3]
   + p[2] * z_p_prime_base_123[3] ] <= 0
 infeasibility_p_prime_base_123_row_4_leq_0: 5e-06 z_p_prime_base_123[4]
   + [ 0.5 alpha[0] * z_p_prime_base_123[4]
   + 0.5 alpha[1] * z_p_prime_base_123[4]
   + 0.5 alpha[2] * z_p_prime_base_123[4]
   - 0.5 alpha[3] * z_p_prime_base_123[4]
   + 0.5 alpha[4] * z_p_prime_base_123[4]
   + 0.5 p[0] * z_p_prime_base_123[4] - 0.5 p[1] * z_p_prime_base_123[4]
   - 0.5 p[2] * z_p_prime_base_123[4] + 0.5 p[3] * z_p_prime_base_123[4]
   + 0.5 p[4] * z_p_prime_base_123[4] ] <= 0
 infeasibility_p_base_124_row_0_leq_0: 5e-06 z_p_base_124[0] + [
   p[0] * z_p_base_124[0] ] <= 0
 infeasibility_p_base_124_row_1_leq_0: 5e-06 z_p_base_124[1] + [
   p[1] * z_p_base_124[1] ] <= 0
 infeasibility_p_base_124_row_2_leq_0: 5e-06 z_p_base_124[2] + [
   - 0.5 p[0] * z_p_base_124[2] - 0.5 p[1] * z_p_base_124[2]
   + 0.5 p[2] * z_p_base_124[2] + 0.5 p[3] * z_p_base_124[2]
   - 0.5 p[4] * z_p_base_124[2] ] <= 0
 infeasibility_p_base_124_row_3_leq_0: 5e-06 z_p_base_124[3] + [
   - 0.5 p[0] * z_p_base_124[3] + 0.5 p[1] * z_p_base_124[3]
   + 0.5 p[2] * z_p_base_124[3] - 0.5 p[3] * z_p_base_124[3]
   + 0.5 p[4] * z_p_base_124[3] ] <= 0
 infeasibility_p_base_124_row_4_leq_0: 5e-06 z_p_base_124[4] + [
   0.5 p[0] * z_p_base_124[4] - 0.5 p[1] * z_p_base_124[4]
   - 0.5 p[2] * z_p_base_124[4] + 0.5 p[3] * z_p_base_124[4]
   + 0.5 p[4] * z_p_base_124[4] ] <= 0
 infeasibility_p_prime_base_124_row_0_leq_0: 5e-06 z_p_prime_base_124[0]
   + [ alpha[0] * z_p_prime_base_124[0] + p[0] * z_p_prime_base_124[0] ]
   <= 0
 infeasibility_p_prime_base_124_row_1_leq_0: 5e-06 z_p_prime_base_124[1]
   + [ - alpha[1] * z_p_prime_base_124[1] + p[1] * z_p_prime_base_124[1] ]
   <= 0
 infeasibility_p_prime_base_124_row_2_leq_0: 5e-06 z_p_prime_base_124[2]
   + [ - 0.5 alpha[0] * z_p_prime_base_124[2]
   + 0.5 alpha[1] * z_p_prime_base_124[2]
   - 0.5 alpha[2] * z_p_prime_base_124[2]
   - 0.5 alpha[3] * z_p_prime_base_124[2]
   - 0.5 alpha[4] * z_p_prime_base_124[2]
   - 0.5 p[0] * z_p_prime_base_124[2] - 0.5 p[1] * z_p_prime_base_124[2]
   + 0.5 p[2] * z_p_prime_base_124[2] + 0.5 p[3] * z_p_prime_base_124[2]
   - 0.5 p[4] * z_p_prime_base_124[2] ] <= 0
 infeasibility_p_prime_base_124_row_3_leq_0: 5e-06 z_p_prime_base_124[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_124[3]
   - 0.5 alpha[1] * z_p_prime_base_124[3]
   - 0.5 alpha[2] * z_p_prime_base_124[3]
   + 0.5 alpha[3] * z_p_prime_base_124[3]
   + 0.5 alpha[4] * z_p_prime_base_124[3]
   - 0.5 p[0] * z_p_prime_base_124[3] + 0.5 p[1] * z_p_prime_base_124[3]
   + 0.5 p[2] * z_p_prime_base_124[3] - 0.5 p[3] * z_p_prime_base_124[3]
   + 0.5 p[4] * z_p_prime_base_124[3] ] <= 0
 infeasibility_p_prime_base_124_row_4_leq_0: 5e-06 z_p_prime_base_124[4]
   + [ 0.5 alpha[0] * z_p_prime_base_124[4]
   + 0.5 alpha[1] * z_p_prime_base_124[4]
   + 0.5 alpha[2] * z_p_prime_base_124[4]
   - 0.5 alpha[3] * z_p_prime_base_124[4]
   + 0.5 alpha[4] * z_p_prime_base_124[4]
   + 0.5 p[0] * z_p_prime_base_124[4] - 0.5 p[1] * z_p_prime_base_124[4]
   - 0.5 p[2] * z_p_prime_base_124[4] + 0.5 p[3] * z_p_prime_base_124[4]
   + 0.5 p[4] * z_p_prime_base_124[4] ] <= 0
 infeasibility_p_base_125_row_0_leq_0: 5e-06 z_p_base_125[0] + [
   p[0] * z_p_base_125[0] ] <= 0
 infeasibility_p_base_125_row_1_leq_0: 5e-06 z_p_base_125[1] + [
   p[1] * z_p_base_125[1] ] <= 0
 infeasibility_p_base_125_row_2_leq_0: 5e-06 z_p_base_125[2] + [
   - 0.5 p[0] * z_p_base_125[2] + 0.5 p[1] * z_p_base_125[2]
   + 0.5 p[2] * z_p_base_125[2] + 0.5 p[3] * z_p_base_125[2]
   - 0.5 p[4] * z_p_base_125[2] ] <= 0
 infeasibility_p_base_125_row_3_leq_0: 5e-06 z_p_base_125[3] + [
   - 0.5 p[0] * z_p_base_125[3] - 0.5 p[1] * z_p_base_125[3]
   + 0.5 p[2] * z_p_base_125[3] - 0.5 p[3] * z_p_base_125[3]
   + 0.5 p[4] * z_p_base_125[3] ] <= 0
 infeasibility_p_base_125_row_4_leq_0: 5e-06 z_p_base_125[4] + [
   0.5 p[0] * z_p_base_125[4] - 0.5 p[1] * z_p_base_125[4]
   - 0.5 p[2] * z_p_base_125[4] + 0.5 p[3] * z_p_base_125[4]
   + 0.5 p[4] * z_p_base_125[4] ] <= 0
 infeasibility_p_prime_base_125_row_0_leq_0: 5e-06 z_p_prime_base_125[0]
   + [ alpha[0] * z_p_prime_base_125[0] + p[0] * z_p_prime_base_125[0] ]
   <= 0
 infeasibility_p_prime_base_125_row_1_leq_0: 5e-06 z_p_prime_base_125[1]
   + [ - alpha[1] * z_p_prime_base_125[1] + p[1] * z_p_prime_base_125[1] ]
   <= 0
 infeasibility_p_prime_base_125_row_2_leq_0: 5e-06 z_p_prime_base_125[2]
   + [ - 0.5 alpha[0] * z_p_prime_base_125[2]
   - 0.5 alpha[1] * z_p_prime_base_125[2]
   - 0.5 alpha[2] * z_p_prime_base_125[2]
   - 0.5 alpha[3] * z_p_prime_base_125[2]
   - 0.5 alpha[4] * z_p_prime_base_125[2]
   - 0.5 p[0] * z_p_prime_base_125[2] + 0.5 p[1] * z_p_prime_base_125[2]
   + 0.5 p[2] * z_p_prime_base_125[2] + 0.5 p[3] * z_p_prime_base_125[2]
   - 0.5 p[4] * z_p_prime_base_125[2] ] <= 0
 infeasibility_p_prime_base_125_row_3_leq_0: 5e-06 z_p_prime_base_125[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_125[3]
   + 0.5 alpha[1] * z_p_prime_base_125[3]
   - 0.5 alpha[2] * z_p_prime_base_125[3]
   + 0.5 alpha[3] * z_p_prime_base_125[3]
   + 0.5 alpha[4] * z_p_prime_base_125[3]
   - 0.5 p[0] * z_p_prime_base_125[3] - 0.5 p[1] * z_p_prime_base_125[3]
   + 0.5 p[2] * z_p_prime_base_125[3] - 0.5 p[3] * z_p_prime_base_125[3]
   + 0.5 p[4] * z_p_prime_base_125[3] ] <= 0
 infeasibility_p_prime_base_125_row_4_leq_0: 5e-06 z_p_prime_base_125[4]
   + [ 0.5 alpha[0] * z_p_prime_base_125[4]
   + 0.5 alpha[1] * z_p_prime_base_125[4]
   + 0.5 alpha[2] * z_p_prime_base_125[4]
   - 0.5 alpha[3] * z_p_prime_base_125[4]
   + 0.5 alpha[4] * z_p_prime_base_125[4]
   + 0.5 p[0] * z_p_prime_base_125[4] - 0.5 p[1] * z_p_prime_base_125[4]
   - 0.5 p[2] * z_p_prime_base_125[4] + 0.5 p[3] * z_p_prime_base_125[4]
   + 0.5 p[4] * z_p_prime_base_125[4] ] <= 0
 infeasibility_p_base_126_row_0_leq_0: 5e-06 z_p_base_126[0] + [
   p[0] * z_p_base_126[0] - p[4] * z_p_base_126[0] ] <= 0
 infeasibility_p_base_126_row_1_leq_0: 5e-06 z_p_base_126[1] + [
   p[4] * z_p_base_126[1] ] <= 0
 infeasibility_p_base_126_row_2_leq_0: 5e-06 z_p_base_126[2] + [
   0.5 p[0] * z_p_base_126[2] + 0.5 p[1] * z_p_base_126[2]
   + 0.5 p[2] * z_p_base_126[2] - 0.5 p[3] * z_p_base_126[2]
   - 0.5 p[4] * z_p_base_126[2] ] <= 0
 infeasibility_p_base_126_row_3_leq_0: 5e-06 z_p_base_126[3] + [
   - 0.5 p[0] * z_p_base_126[3] + 0.5 p[1] * z_p_base_126[3]
   - 0.5 p[2] * z_p_base_126[3] + 0.5 p[3] * z_p_base_126[3]
   + 0.5 p[4] * z_p_base_126[3] ] <= 0
 infeasibility_p_base_126_row_4_leq_0: 5e-06 z_p_base_126[4] + [
   - 0.5 p[0] * z_p_base_126[4] - 0.5 p[1] * z_p_base_126[4]
   + 0.5 p[2] * z_p_base_126[4] + 0.5 p[3] * z_p_base_126[4]
   + 0.5 p[4] * z_p_base_126[4] ] <= 0
 infeasibility_p_prime_base_126_row_0_leq_0: 5e-06 z_p_prime_base_126[0]
   + [ alpha[0] * z_p_prime_base_126[0] - alpha[4] * z_p_prime_base_126[0]
   + p[0] * z_p_prime_base_126[0] - p[4] * z_p_prime_base_126[0] ] <= 0
 infeasibility_p_prime_base_126_row_1_leq_0: 5e-06 z_p_prime_base_126[1]
   + [ alpha[4] * z_p_prime_base_126[1] + p[4] * z_p_prime_base_126[1] ]
   <= 0
 infeasibility_p_prime_base_126_row_2_leq_0: 5e-06 z_p_prime_base_126[2]
   + [ 0.5 alpha[0] * z_p_prime_base_126[2]
   - 0.5 alpha[1] * z_p_prime_base_126[2]
   - 0.5 alpha[2] * z_p_prime_base_126[2]
   + 0.5 alpha[3] * z_p_prime_base_126[2]
   - 0.5 alpha[4] * z_p_prime_base_126[2]
   + 0.5 p[0] * z_p_prime_base_126[2] + 0.5 p[1] * z_p_prime_base_126[2]
   + 0.5 p[2] * z_p_prime_base_126[2] - 0.5 p[3] * z_p_prime_base_126[2]
   - 0.5 p[4] * z_p_prime_base_126[2] ] <= 0
 infeasibility_p_prime_base_126_row_3_leq_0: 5e-06 z_p_prime_base_126[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_126[3]
   - 0.5 alpha[1] * z_p_prime_base_126[3]
   + 0.5 alpha[2] * z_p_prime_base_126[3]
   - 0.5 alpha[3] * z_p_prime_base_126[3]
   + 0.5 alpha[4] * z_p_prime_base_126[3]
   - 0.5 p[0] * z_p_prime_base_126[3] + 0.5 p[1] * z_p_prime_base_126[3]
   - 0.5 p[2] * z_p_prime_base_126[3] + 0.5 p[3] * z_p_prime_base_126[3]
   + 0.5 p[4] * z_p_prime_base_126[3] ] <= 0
 infeasibility_p_prime_base_126_row_4_leq_0: 5e-06 z_p_prime_base_126[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_126[4]
   + 0.5 alpha[1] * z_p_prime_base_126[4]
   - 0.5 alpha[2] * z_p_prime_base_126[4]
   - 0.5 alpha[3] * z_p_prime_base_126[4]
   + 0.5 alpha[4] * z_p_prime_base_126[4]
   - 0.5 p[0] * z_p_prime_base_126[4] - 0.5 p[1] * z_p_prime_base_126[4]
   + 0.5 p[2] * z_p_prime_base_126[4] + 0.5 p[3] * z_p_prime_base_126[4]
   + 0.5 p[4] * z_p_prime_base_126[4] ] <= 0
 infeasibility_p_base_127_row_0_leq_0: 5e-06 z_p_base_127[0] + [
   0.5 p[0] * z_p_base_127[0] - 0.5 p[1] * z_p_base_127[0]
   + 0.5 p[2] * z_p_base_127[0] + 0.5 p[3] * z_p_base_127[0]
   - 0.5 p[4] * z_p_base_127[0] ] <= 0
 infeasibility_p_base_127_row_1_leq_0: 5e-06 z_p_base_127[1] + [
   0.5 p[0] * z_p_base_127[1] + 0.5 p[1] * z_p_base_127[1]
   - 0.5 p[2] * z_p_base_127[1] - 0.5 p[3] * z_p_base_127[1]
   + 0.5 p[4] * z_p_base_127[1] ] <= 0
 infeasibility_p_base_127_row_2_leq_0: 5e-06 z_p_base_127[2] + [
   0.5 p[0] * z_p_base_127[2] + 0.5 p[1] * z_p_base_127[2]
   + 0.5 p[2] * z_p_base_127[2] - 0.5 p[3] * z_p_base_127[2]
   - 0.5 p[4] * z_p_base_127[2] ] <= 0
 infeasibility_p_base_127_row_3_leq_0: 5e-06 z_p_base_127[3] + [
   - 0.5 p[0] * z_p_base_127[3] + 0.5 p[1] * z_p_base_127[3]
   - 0.5 p[2] * z_p_base_127[3] + 0.5 p[3] * z_p_base_127[3]
   + 0.5 p[4] * z_p_base_127[3] ] <= 0
 infeasibility_p_base_127_row_4_leq_0: 5e-06 z_p_base_127[4] + [
   - 0.5 p[0] * z_p_base_127[4] - 0.5 p[1] * z_p_base_127[4]
   + 0.5 p[2] * z_p_base_127[4] + 0.5 p[3] * z_p_base_127[4]
   + 0.5 p[4] * z_p_base_127[4] ] <= 0
 infeasibility_p_prime_base_127_row_0_leq_0: 5e-06 z_p_prime_base_127[0]
   + [ 0.5 alpha[0] * z_p_prime_base_127[0]
   + 0.5 alpha[1] * z_p_prime_base_127[0]
   - 0.5 alpha[2] * z_p_prime_base_127[0]
   - 0.5 alpha[3] * z_p_prime_base_127[0]
   - 0.5 alpha[4] * z_p_prime_base_127[0]
   + 0.5 p[0] * z_p_prime_base_127[0] - 0.5 p[1] * z_p_prime_base_127[0]
   + 0.5 p[2] * z_p_prime_base_127[0] + 0.5 p[3] * z_p_prime_base_127[0]
   - 0.5 p[4] * z_p_prime_base_127[0] ] <= 0
 infeasibility_p_prime_base_127_row_1_leq_0: 5e-06 z_p_prime_base_127[1]
   + [ 0.5 alpha[0] * z_p_prime_base_127[1]
   - 0.5 alpha[1] * z_p_prime_base_127[1]
   + 0.5 alpha[2] * z_p_prime_base_127[1]
   + 0.5 alpha[3] * z_p_prime_base_127[1]
   + 0.5 alpha[4] * z_p_prime_base_127[1]
   + 0.5 p[0] * z_p_prime_base_127[1] + 0.5 p[1] * z_p_prime_base_127[1]
   - 0.5 p[2] * z_p_prime_base_127[1] - 0.5 p[3] * z_p_prime_base_127[1]
   + 0.5 p[4] * z_p_prime_base_127[1] ] <= 0
 infeasibility_p_prime_base_127_row_2_leq_0: 5e-06 z_p_prime_base_127[2]
   + [ 0.5 alpha[0] * z_p_prime_base_127[2]
   - 0.5 alpha[1] * z_p_prime_base_127[2]
   - 0.5 alpha[2] * z_p_prime_base_127[2]
   + 0.5 alpha[3] * z_p_prime_base_127[2]
   - 0.5 alpha[4] * z_p_prime_base_127[2]
   + 0.5 p[0] * z_p_prime_base_127[2] + 0.5 p[1] * z_p_prime_base_127[2]
   + 0.5 p[2] * z_p_prime_base_127[2] - 0.5 p[3] * z_p_prime_base_127[2]
   - 0.5 p[4] * z_p_prime_base_127[2] ] <= 0
 infeasibility_p_prime_base_127_row_3_leq_0: 5e-06 z_p_prime_base_127[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_127[3]
   - 0.5 alpha[1] * z_p_prime_base_127[3]
   + 0.5 alpha[2] * z_p_prime_base_127[3]
   - 0.5 alpha[3] * z_p_prime_base_127[3]
   + 0.5 alpha[4] * z_p_prime_base_127[3]
   - 0.5 p[0] * z_p_prime_base_127[3] + 0.5 p[1] * z_p_prime_base_127[3]
   - 0.5 p[2] * z_p_prime_base_127[3] + 0.5 p[3] * z_p_prime_base_127[3]
   + 0.5 p[4] * z_p_prime_base_127[3] ] <= 0
 infeasibility_p_prime_base_127_row_4_leq_0: 5e-06 z_p_prime_base_127[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_127[4]
   + 0.5 alpha[1] * z_p_prime_base_127[4]
   - 0.5 alpha[2] * z_p_prime_base_127[4]
   - 0.5 alpha[3] * z_p_prime_base_127[4]
   + 0.5 alpha[4] * z_p_prime_base_127[4]
   - 0.5 p[0] * z_p_prime_base_127[4] - 0.5 p[1] * z_p_prime_base_127[4]
   + 0.5 p[2] * z_p_prime_base_127[4] + 0.5 p[3] * z_p_prime_base_127[4]
   + 0.5 p[4] * z_p_prime_base_127[4] ] <= 0
 infeasibility_p_base_128_row_0_leq_0: 5e-06 z_p_base_128[0] + [
   0.5 p[0] * z_p_base_128[0] - 0.5 p[1] * z_p_base_128[0]
   + 0.5 p[2] * z_p_base_128[0] + 0.5 p[3] * z_p_base_128[0]
   - 0.5 p[4] * z_p_base_128[0] ] <= 0
 infeasibility_p_base_128_row_1_leq_0: 5e-06 z_p_base_128[1] + [
   0.5 p[0] * z_p_base_128[1] + 0.5 p[1] * z_p_base_128[1]
   - 0.5 p[2] * z_p_base_128[1] - 0.5 p[3] * z_p_base_128[1]
   + 0.5 p[4] * z_p_base_128[1] ] <= 0
 infeasibility_p_base_128_row_2_leq_0: 5e-06 z_p_base_128[2] + [
   p[2] * z_p_base_128[2] ] <= 0
 infeasibility_p_base_128_row_3_leq_0: 5e-06 z_p_base_128[3] + [
   p[1] * z_p_base_128[3] - p[2] * z_p_base_128[3] ] <= 0
 infeasibility_p_base_128_row_4_leq_0: 5e-06 z_p_base_128[4] + [
   - 0.5 p[0] * z_p_base_128[4] - 0.5 p[1] * z_p_base_128[4]
   + 0.5 p[2] * z_p_base_128[4] + 0.5 p[3] * z_p_base_128[4]
   + 0.5 p[4] * z_p_base_128[4] ] <= 0
 infeasibility_p_prime_base_128_row_0_leq_0: 5e-06 z_p_prime_base_128[0]
   + [ 0.5 alpha[0] * z_p_prime_base_128[0]
   + 0.5 alpha[1] * z_p_prime_base_128[0]
   - 0.5 alpha[2] * z_p_prime_base_128[0]
   - 0.5 alpha[3] * z_p_prime_base_128[0]
   - 0.5 alpha[4] * z_p_prime_base_128[0]
   + 0.5 p[0] * z_p_prime_base_128[0] - 0.5 p[1] * z_p_prime_base_128[0]
   + 0.5 p[2] * z_p_prime_base_128[0] + 0.5 p[3] * z_p_prime_base_128[0]
   - 0.5 p[4] * z_p_prime_base_128[0] ] <= 0
 infeasibility_p_prime_base_128_row_1_leq_0: 5e-06 z_p_prime_base_128[1]
   + [ 0.5 alpha[0] * z_p_prime_base_128[1]
   - 0.5 alpha[1] * z_p_prime_base_128[1]
   + 0.5 alpha[2] * z_p_prime_base_128[1]
   + 0.5 alpha[3] * z_p_prime_base_128[1]
   + 0.5 alpha[4] * z_p_prime_base_128[1]
   + 0.5 p[0] * z_p_prime_base_128[1] + 0.5 p[1] * z_p_prime_base_128[1]
   - 0.5 p[2] * z_p_prime_base_128[1] - 0.5 p[3] * z_p_prime_base_128[1]
   + 0.5 p[4] * z_p_prime_base_128[1] ] <= 0
 infeasibility_p_prime_base_128_row_2_leq_0: 5e-06 z_p_prime_base_128[2]
   + [ - alpha[2] * z_p_prime_base_128[2] + p[2] * z_p_prime_base_128[2] ]
   <= 0
 infeasibility_p_prime_base_128_row_3_leq_0: 5e-06 z_p_prime_base_128[3]
   + [ - alpha[1] * z_p_prime_base_128[3]
   + alpha[2] * z_p_prime_base_128[3] + p[1] * z_p_prime_base_128[3]
   - p[2] * z_p_prime_base_128[3] ] <= 0
 infeasibility_p_prime_base_128_row_4_leq_0: 5e-06 z_p_prime_base_128[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_128[4]
   + 0.5 alpha[1] * z_p_prime_base_128[4]
   - 0.5 alpha[2] * z_p_prime_base_128[4]
   - 0.5 alpha[3] * z_p_prime_base_128[4]
   + 0.5 alpha[4] * z_p_prime_base_128[4]
   - 0.5 p[0] * z_p_prime_base_128[4] - 0.5 p[1] * z_p_prime_base_128[4]
   + 0.5 p[2] * z_p_prime_base_128[4] + 0.5 p[3] * z_p_prime_base_128[4]
   + 0.5 p[4] * z_p_prime_base_128[4] ] <= 0
 infeasibility_p_base_129_row_0_leq_0: 5e-06 z_p_base_129[0] + [
   0.5 p[0] * z_p_base_129[0] + 0.5 p[1] * z_p_base_129[0]
   - 0.5 p[2] * z_p_base_129[0] + 0.5 p[3] * z_p_base_129[0]
   - 0.5 p[4] * z_p_base_129[0] ] <= 0
 infeasibility_p_base_129_row_1_leq_0: 5e-06 z_p_base_129[1] + [
   0.5 p[0] * z_p_base_129[1] - 0.5 p[1] * z_p_base_129[1]
   + 0.5 p[2] * z_p_base_129[1] - 0.5 p[3] * z_p_base_129[1]
   + 0.5 p[4] * z_p_base_129[1] ] <= 0
 infeasibility_p_base_129_row_2_leq_0: 5e-06 z_p_base_129[2] + [
   0.5 p[0] * z_p_base_129[2] + 0.5 p[1] * z_p_base_129[2]
   + 0.5 p[2] * z_p_base_129[2] - 0.5 p[3] * z_p_base_129[2]
   - 0.5 p[4] * z_p_base_129[2] ] <= 0
 infeasibility_p_base_129_row_3_leq_0: 5e-06 z_p_base_129[3] + [
   - 0.5 p[0] * z_p_base_129[3] + 0.5 p[1] * z_p_base_129[3]
   - 0.5 p[2] * z_p_base_129[3] + 0.5 p[3] * z_p_base_129[3]
   + 0.5 p[4] * z_p_base_129[3] ] <= 0
 infeasibility_p_base_129_row_4_leq_0: 5e-06 z_p_base_129[4] + [
   - 0.5 p[0] * z_p_base_129[4] - 0.5 p[1] * z_p_base_129[4]
   + 0.5 p[2] * z_p_base_129[4] + 0.5 p[3] * z_p_base_129[4]
   + 0.5 p[4] * z_p_base_129[4] ] <= 0
 infeasibility_p_prime_base_129_row_0_leq_0: 5e-06 z_p_prime_base_129[0]
   + [ 0.5 alpha[0] * z_p_prime_base_129[0]
   - 0.5 alpha[1] * z_p_prime_base_129[0]
   + 0.5 alpha[2] * z_p_prime_base_129[0]
   - 0.5 alpha[3] * z_p_prime_base_129[0]
   - 0.5 alpha[4] * z_p_prime_base_129[0]
   + 0.5 p[0] * z_p_prime_base_129[0] + 0.5 p[1] * z_p_prime_base_129[0]
   - 0.5 p[2] * z_p_prime_base_129[0] + 0.5 p[3] * z_p_prime_base_129[0]
   - 0.5 p[4] * z_p_prime_base_129[0] ] <= 0
 infeasibility_p_prime_base_129_row_1_leq_0: 5e-06 z_p_prime_base_129[1]
   + [ 0.5 alpha[0] * z_p_prime_base_129[1]
   + 0.5 alpha[1] * z_p_prime_base_129[1]
   - 0.5 alpha[2] * z_p_prime_base_129[1]
   + 0.5 alpha[3] * z_p_prime_base_129[1]
   + 0.5 alpha[4] * z_p_prime_base_129[1]
   + 0.5 p[0] * z_p_prime_base_129[1] - 0.5 p[1] * z_p_prime_base_129[1]
   + 0.5 p[2] * z_p_prime_base_129[1] - 0.5 p[3] * z_p_prime_base_129[1]
   + 0.5 p[4] * z_p_prime_base_129[1] ] <= 0
 infeasibility_p_prime_base_129_row_2_leq_0: 5e-06 z_p_prime_base_129[2]
   + [ 0.5 alpha[0] * z_p_prime_base_129[2]
   - 0.5 alpha[1] * z_p_prime_base_129[2]
   - 0.5 alpha[2] * z_p_prime_base_129[2]
   + 0.5 alpha[3] * z_p_prime_base_129[2]
   - 0.5 alpha[4] * z_p_prime_base_129[2]
   + 0.5 p[0] * z_p_prime_base_129[2] + 0.5 p[1] * z_p_prime_base_129[2]
   + 0.5 p[2] * z_p_prime_base_129[2] - 0.5 p[3] * z_p_prime_base_129[2]
   - 0.5 p[4] * z_p_prime_base_129[2] ] <= 0
 infeasibility_p_prime_base_129_row_3_leq_0: 5e-06 z_p_prime_base_129[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_129[3]
   - 0.5 alpha[1] * z_p_prime_base_129[3]
   + 0.5 alpha[2] * z_p_prime_base_129[3]
   - 0.5 alpha[3] * z_p_prime_base_129[3]
   + 0.5 alpha[4] * z_p_prime_base_129[3]
   - 0.5 p[0] * z_p_prime_base_129[3] + 0.5 p[1] * z_p_prime_base_129[3]
   - 0.5 p[2] * z_p_prime_base_129[3] + 0.5 p[3] * z_p_prime_base_129[3]
   + 0.5 p[4] * z_p_prime_base_129[3] ] <= 0
 infeasibility_p_prime_base_129_row_4_leq_0: 5e-06 z_p_prime_base_129[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_129[4]
   + 0.5 alpha[1] * z_p_prime_base_129[4]
   - 0.5 alpha[2] * z_p_prime_base_129[4]
   - 0.5 alpha[3] * z_p_prime_base_129[4]
   + 0.5 alpha[4] * z_p_prime_base_129[4]
   - 0.5 p[0] * z_p_prime_base_129[4] - 0.5 p[1] * z_p_prime_base_129[4]
   + 0.5 p[2] * z_p_prime_base_129[4] + 0.5 p[3] * z_p_prime_base_129[4]
   + 0.5 p[4] * z_p_prime_base_129[4] ] <= 0
 infeasibility_p_base_130_row_0_leq_0: 5e-06 z_p_base_130[0] + [
   p[3] * z_p_base_130[0] ] <= 0
 infeasibility_p_base_130_row_1_leq_0: 5e-06 z_p_base_130[1] + [
   p[0] * z_p_base_130[1] - p[3] * z_p_base_130[1] ] <= 0
 infeasibility_p_base_130_row_2_leq_0: 5e-06 z_p_base_130[2] + [
   0.5 p[0] * z_p_base_130[2] + 0.5 p[1] * z_p_base_130[2]
   + 0.5 p[2] * z_p_base_130[2] - 0.5 p[3] * z_p_base_130[2]
   - 0.5 p[4] * z_p_base_130[2] ] <= 0
 infeasibility_p_base_130_row_3_leq_0: 5e-06 z_p_base_130[3] + [
   - 0.5 p[0] * z_p_base_130[3] + 0.5 p[1] * z_p_base_130[3]
   - 0.5 p[2] * z_p_base_130[3] + 0.5 p[3] * z_p_base_130[3]
   + 0.5 p[4] * z_p_base_130[3] ] <= 0
 infeasibility_p_base_130_row_4_leq_0: 5e-06 z_p_base_130[4] + [
   - 0.5 p[0] * z_p_base_130[4] - 0.5 p[1] * z_p_base_130[4]
   + 0.5 p[2] * z_p_base_130[4] + 0.5 p[3] * z_p_base_130[4]
   + 0.5 p[4] * z_p_base_130[4] ] <= 0
 infeasibility_p_prime_base_130_row_0_leq_0: 5e-06 z_p_prime_base_130[0]
   + [ - alpha[3] * z_p_prime_base_130[0] + p[3] * z_p_prime_base_130[0] ]
   <= 0
 infeasibility_p_prime_base_130_row_1_leq_0: 5e-06 z_p_prime_base_130[1]
   + [ alpha[0] * z_p_prime_base_130[1] + alpha[3] * z_p_prime_base_130[1]
   + p[0] * z_p_prime_base_130[1] - p[3] * z_p_prime_base_130[1] ] <= 0
 infeasibility_p_prime_base_130_row_2_leq_0: 5e-06 z_p_prime_base_130[2]
   + [ 0.5 alpha[0] * z_p_prime_base_130[2]
   - 0.5 alpha[1] * z_p_prime_base_130[2]
   - 0.5 alpha[2] * z_p_prime_base_130[2]
   + 0.5 alpha[3] * z_p_prime_base_130[2]
   - 0.5 alpha[4] * z_p_prime_base_130[2]
   + 0.5 p[0] * z_p_prime_base_130[2] + 0.5 p[1] * z_p_prime_base_130[2]
   + 0.5 p[2] * z_p_prime_base_130[2] - 0.5 p[3] * z_p_prime_base_130[2]
   - 0.5 p[4] * z_p_prime_base_130[2] ] <= 0
 infeasibility_p_prime_base_130_row_3_leq_0: 5e-06 z_p_prime_base_130[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_130[3]
   - 0.5 alpha[1] * z_p_prime_base_130[3]
   + 0.5 alpha[2] * z_p_prime_base_130[3]
   - 0.5 alpha[3] * z_p_prime_base_130[3]
   + 0.5 alpha[4] * z_p_prime_base_130[3]
   - 0.5 p[0] * z_p_prime_base_130[3] + 0.5 p[1] * z_p_prime_base_130[3]
   - 0.5 p[2] * z_p_prime_base_130[3] + 0.5 p[3] * z_p_prime_base_130[3]
   + 0.5 p[4] * z_p_prime_base_130[3] ] <= 0
 infeasibility_p_prime_base_130_row_4_leq_0: 5e-06 z_p_prime_base_130[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_130[4]
   + 0.5 alpha[1] * z_p_prime_base_130[4]
   - 0.5 alpha[2] * z_p_prime_base_130[4]
   - 0.5 alpha[3] * z_p_prime_base_130[4]
   + 0.5 alpha[4] * z_p_prime_base_130[4]
   - 0.5 p[0] * z_p_prime_base_130[4] - 0.5 p[1] * z_p_prime_base_130[4]
   + 0.5 p[2] * z_p_prime_base_130[4] + 0.5 p[3] * z_p_prime_base_130[4]
   + 0.5 p[4] * z_p_prime_base_130[4] ] <= 0
 infeasibility_p_base_131_row_0_leq_0: 5e-06 z_p_base_131[0] + [
   0.5 p[0] * z_p_base_131[0] + 0.5 p[1] * z_p_base_131[0]
   - 0.5 p[2] * z_p_base_131[0] + 0.5 p[3] * z_p_base_131[0]
   - 0.5 p[4] * z_p_base_131[0] ] <= 0
 infeasibility_p_base_131_row_1_leq_0: 5e-06 z_p_base_131[1] + [
   0.5 p[0] * z_p_base_131[1] - 0.5 p[1] * z_p_base_131[1]
   + 0.5 p[2] * z_p_base_131[1] - 0.5 p[3] * z_p_base_131[1]
   + 0.5 p[4] * z_p_base_131[1] ] <= 0
 infeasibility_p_base_131_row_2_leq_0: 5e-06 z_p_base_131[2] + [
   p[2] * z_p_base_131[2] ] <= 0
 infeasibility_p_base_131_row_3_leq_0: 5e-06 z_p_base_131[3] + [
   p[1] * z_p_base_131[3] - p[2] * z_p_base_131[3] ] <= 0
 infeasibility_p_base_131_row_4_leq_0: 5e-06 z_p_base_131[4] + [
   - 0.5 p[0] * z_p_base_131[4] - 0.5 p[1] * z_p_base_131[4]
   + 0.5 p[2] * z_p_base_131[4] + 0.5 p[3] * z_p_base_131[4]
   + 0.5 p[4] * z_p_base_131[4] ] <= 0
 infeasibility_p_prime_base_131_row_0_leq_0: 5e-06 z_p_prime_base_131[0]
   + [ 0.5 alpha[0] * z_p_prime_base_131[0]
   - 0.5 alpha[1] * z_p_prime_base_131[0]
   + 0.5 alpha[2] * z_p_prime_base_131[0]
   - 0.5 alpha[3] * z_p_prime_base_131[0]
   - 0.5 alpha[4] * z_p_prime_base_131[0]
   + 0.5 p[0] * z_p_prime_base_131[0] + 0.5 p[1] * z_p_prime_base_131[0]
   - 0.5 p[2] * z_p_prime_base_131[0] + 0.5 p[3] * z_p_prime_base_131[0]
   - 0.5 p[4] * z_p_prime_base_131[0] ] <= 0
 infeasibility_p_prime_base_131_row_1_leq_0: 5e-06 z_p_prime_base_131[1]
   + [ 0.5 alpha[0] * z_p_prime_base_131[1]
   + 0.5 alpha[1] * z_p_prime_base_131[1]
   - 0.5 alpha[2] * z_p_prime_base_131[1]
   + 0.5 alpha[3] * z_p_prime_base_131[1]
   + 0.5 alpha[4] * z_p_prime_base_131[1]
   + 0.5 p[0] * z_p_prime_base_131[1] - 0.5 p[1] * z_p_prime_base_131[1]
   + 0.5 p[2] * z_p_prime_base_131[1] - 0.5 p[3] * z_p_prime_base_131[1]
   + 0.5 p[4] * z_p_prime_base_131[1] ] <= 0
 infeasibility_p_prime_base_131_row_2_leq_0: 5e-06 z_p_prime_base_131[2]
   + [ - alpha[2] * z_p_prime_base_131[2] + p[2] * z_p_prime_base_131[2] ]
   <= 0
 infeasibility_p_prime_base_131_row_3_leq_0: 5e-06 z_p_prime_base_131[3]
   + [ - alpha[1] * z_p_prime_base_131[3]
   + alpha[2] * z_p_prime_base_131[3] + p[1] * z_p_prime_base_131[3]
   - p[2] * z_p_prime_base_131[3] ] <= 0
 infeasibility_p_prime_base_131_row_4_leq_0: 5e-06 z_p_prime_base_131[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_131[4]
   + 0.5 alpha[1] * z_p_prime_base_131[4]
   - 0.5 alpha[2] * z_p_prime_base_131[4]
   - 0.5 alpha[3] * z_p_prime_base_131[4]
   + 0.5 alpha[4] * z_p_prime_base_131[4]
   - 0.5 p[0] * z_p_prime_base_131[4] - 0.5 p[1] * z_p_prime_base_131[4]
   + 0.5 p[2] * z_p_prime_base_131[4] + 0.5 p[3] * z_p_prime_base_131[4]
   + 0.5 p[4] * z_p_prime_base_131[4] ] <= 0
 infeasibility_p_base_132_row_0_leq_0: 5e-06 z_p_base_132[0] + [
   0.5 p[0] * z_p_base_132[0] + 0.5 p[1] * z_p_base_132[0]
   - 0.5 p[2] * z_p_base_132[0] + 0.5 p[3] * z_p_base_132[0]
   - 0.5 p[4] * z_p_base_132[0] ] <= 0
 infeasibility_p_base_132_row_1_leq_0: 5e-06 z_p_base_132[1] + [
   0.5 p[0] * z_p_base_132[1] - 0.5 p[1] * z_p_base_132[1]
   + 0.5 p[2] * z_p_base_132[1] - 0.5 p[3] * z_p_base_132[1]
   + 0.5 p[4] * z_p_base_132[1] ] <= 0
 infeasibility_p_base_132_row_2_leq_0: 5e-06 z_p_base_132[2] + [
   p[1] * z_p_base_132[2] ] <= 0
 infeasibility_p_base_132_row_3_leq_0: 5e-06 z_p_base_132[3] + [
   - p[1] * z_p_base_132[3] + p[2] * z_p_base_132[3] ] <= 0
 infeasibility_p_base_132_row_4_leq_0: 5e-06 z_p_base_132[4] + [
   - 0.5 p[0] * z_p_base_132[4] + 0.5 p[1] * z_p_base_132[4]
   - 0.5 p[2] * z_p_base_132[4] + 0.5 p[3] * z_p_base_132[4]
   + 0.5 p[4] * z_p_base_132[4] ] <= 0
 infeasibility_p_prime_base_132_row_0_leq_0: 5e-06 z_p_prime_base_132[0]
   + [ 0.5 alpha[0] * z_p_prime_base_132[0]
   - 0.5 alpha[1] * z_p_prime_base_132[0]
   + 0.5 alpha[2] * z_p_prime_base_132[0]
   - 0.5 alpha[3] * z_p_prime_base_132[0]
   - 0.5 alpha[4] * z_p_prime_base_132[0]
   + 0.5 p[0] * z_p_prime_base_132[0] + 0.5 p[1] * z_p_prime_base_132[0]
   - 0.5 p[2] * z_p_prime_base_132[0] + 0.5 p[3] * z_p_prime_base_132[0]
   - 0.5 p[4] * z_p_prime_base_132[0] ] <= 0
 infeasibility_p_prime_base_132_row_1_leq_0: 5e-06 z_p_prime_base_132[1]
   + [ 0.5 alpha[0] * z_p_prime_base_132[1]
   + 0.5 alpha[1] * z_p_prime_base_132[1]
   - 0.5 alpha[2] * z_p_prime_base_132[1]
   + 0.5 alpha[3] * z_p_prime_base_132[1]
   + 0.5 alpha[4] * z_p_prime_base_132[1]
   + 0.5 p[0] * z_p_prime_base_132[1] - 0.5 p[1] * z_p_prime_base_132[1]
   + 0.5 p[2] * z_p_prime_base_132[1] - 0.5 p[3] * z_p_prime_base_132[1]
   + 0.5 p[4] * z_p_prime_base_132[1] ] <= 0
 infeasibility_p_prime_base_132_row_2_leq_0: 5e-06 z_p_prime_base_132[2]
   + [ - alpha[1] * z_p_prime_base_132[2] + p[1] * z_p_prime_base_132[2] ]
   <= 0
 infeasibility_p_prime_base_132_row_3_leq_0: 5e-06 z_p_prime_base_132[3]
   + [ alpha[1] * z_p_prime_base_132[3] - alpha[2] * z_p_prime_base_132[3]
   - p[1] * z_p_prime_base_132[3] + p[2] * z_p_prime_base_132[3] ] <= 0
 infeasibility_p_prime_base_132_row_4_leq_0: 5e-06 z_p_prime_base_132[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_132[4]
   - 0.5 alpha[1] * z_p_prime_base_132[4]
   + 0.5 alpha[2] * z_p_prime_base_132[4]
   - 0.5 alpha[3] * z_p_prime_base_132[4]
   + 0.5 alpha[4] * z_p_prime_base_132[4]
   - 0.5 p[0] * z_p_prime_base_132[4] + 0.5 p[1] * z_p_prime_base_132[4]
   - 0.5 p[2] * z_p_prime_base_132[4] + 0.5 p[3] * z_p_prime_base_132[4]
   + 0.5 p[4] * z_p_prime_base_132[4] ] <= 0
 infeasibility_p_base_133_row_0_leq_0: 5e-06 z_p_base_133[0] + [
   0.5 p[0] * z_p_base_133[0] - 0.5 p[1] * z_p_base_133[0]
   + 0.5 p[2] * z_p_base_133[0] + 0.5 p[3] * z_p_base_133[0]
   - 0.5 p[4] * z_p_base_133[0] ] <= 0
 infeasibility_p_base_133_row_1_leq_0: 5e-06 z_p_base_133[1] + [
   0.5 p[0] * z_p_base_133[1] + 0.5 p[1] * z_p_base_133[1]
   - 0.5 p[2] * z_p_base_133[1] - 0.5 p[3] * z_p_base_133[1]
   + 0.5 p[4] * z_p_base_133[1] ] <= 0
 infeasibility_p_base_133_row_2_leq_0: 5e-06 z_p_base_133[2] + [
   p[1] * z_p_base_133[2] ] <= 0
 infeasibility_p_base_133_row_3_leq_0: 5e-06 z_p_base_133[3] + [
   - p[1] * z_p_base_133[3] + p[2] * z_p_base_133[3] ] <= 0
 infeasibility_p_base_133_row_4_leq_0: 5e-06 z_p_base_133[4] + [
   - 0.5 p[0] * z_p_base_133[4] + 0.5 p[1] * z_p_base_133[4]
   - 0.5 p[2] * z_p_base_133[4] + 0.5 p[3] * z_p_base_133[4]
   + 0.5 p[4] * z_p_base_133[4] ] <= 0
 infeasibility_p_prime_base_133_row_0_leq_0: 5e-06 z_p_prime_base_133[0]
   + [ 0.5 alpha[0] * z_p_prime_base_133[0]
   + 0.5 alpha[1] * z_p_prime_base_133[0]
   - 0.5 alpha[2] * z_p_prime_base_133[0]
   - 0.5 alpha[3] * z_p_prime_base_133[0]
   - 0.5 alpha[4] * z_p_prime_base_133[0]
   + 0.5 p[0] * z_p_prime_base_133[0] - 0.5 p[1] * z_p_prime_base_133[0]
   + 0.5 p[2] * z_p_prime_base_133[0] + 0.5 p[3] * z_p_prime_base_133[0]
   - 0.5 p[4] * z_p_prime_base_133[0] ] <= 0
 infeasibility_p_prime_base_133_row_1_leq_0: 5e-06 z_p_prime_base_133[1]
   + [ 0.5 alpha[0] * z_p_prime_base_133[1]
   - 0.5 alpha[1] * z_p_prime_base_133[1]
   + 0.5 alpha[2] * z_p_prime_base_133[1]
   + 0.5 alpha[3] * z_p_prime_base_133[1]
   + 0.5 alpha[4] * z_p_prime_base_133[1]
   + 0.5 p[0] * z_p_prime_base_133[1] + 0.5 p[1] * z_p_prime_base_133[1]
   - 0.5 p[2] * z_p_prime_base_133[1] - 0.5 p[3] * z_p_prime_base_133[1]
   + 0.5 p[4] * z_p_prime_base_133[1] ] <= 0
 infeasibility_p_prime_base_133_row_2_leq_0: 5e-06 z_p_prime_base_133[2]
   + [ - alpha[1] * z_p_prime_base_133[2] + p[1] * z_p_prime_base_133[2] ]
   <= 0
 infeasibility_p_prime_base_133_row_3_leq_0: 5e-06 z_p_prime_base_133[3]
   + [ alpha[1] * z_p_prime_base_133[3] - alpha[2] * z_p_prime_base_133[3]
   - p[1] * z_p_prime_base_133[3] + p[2] * z_p_prime_base_133[3] ] <= 0
 infeasibility_p_prime_base_133_row_4_leq_0: 5e-06 z_p_prime_base_133[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_133[4]
   - 0.5 alpha[1] * z_p_prime_base_133[4]
   + 0.5 alpha[2] * z_p_prime_base_133[4]
   - 0.5 alpha[3] * z_p_prime_base_133[4]
   + 0.5 alpha[4] * z_p_prime_base_133[4]
   - 0.5 p[0] * z_p_prime_base_133[4] + 0.5 p[1] * z_p_prime_base_133[4]
   - 0.5 p[2] * z_p_prime_base_133[4] + 0.5 p[3] * z_p_prime_base_133[4]
   + 0.5 p[4] * z_p_prime_base_133[4] ] <= 0
 infeasibility_p_base_134_row_0_leq_0: 5e-06 z_p_base_134[0] + [
   0.5 p[0] * z_p_base_134[0] - 0.5 p[1] * z_p_base_134[0]
   - 0.5 p[2] * z_p_base_134[0] + 0.5 p[3] * z_p_base_134[0]
   - 0.5 p[4] * z_p_base_134[0] ] <= 0
 infeasibility_p_base_134_row_1_leq_0: 5e-06 z_p_base_134[1] + [
   0.5 p[0] * z_p_base_134[1] + 0.5 p[1] * z_p_base_134[1]
   + 0.5 p[2] * z_p_base_134[1] - 0.5 p[3] * z_p_base_134[1]
   + 0.5 p[4] * z_p_base_134[1] ] <= 0
 infeasibility_p_base_134_row_2_leq_0: 5e-06 z_p_base_134[2] + [
   p[1] * z_p_base_134[2] ] <= 0
 infeasibility_p_base_134_row_3_leq_0: 5e-06 z_p_base_134[3] + [
   p[2] * z_p_base_134[3] ] <= 0
 infeasibility_p_base_134_row_4_leq_0: 5e-06 z_p_base_134[4] + [
   - 0.5 p[0] * z_p_base_134[4] - 0.5 p[1] * z_p_base_134[4]
   - 0.5 p[2] * z_p_base_134[4] + 0.5 p[3] * z_p_base_134[4]
   + 0.5 p[4] * z_p_base_134[4] ] <= 0
 infeasibility_p_prime_base_134_row_0_leq_0: 5e-06 z_p_prime_base_134[0]
   + [ 0.5 alpha[0] * z_p_prime_base_134[0]
   + 0.5 alpha[1] * z_p_prime_base_134[0]
   + 0.5 alpha[2] * z_p_prime_base_134[0]
   - 0.5 alpha[3] * z_p_prime_base_134[0]
   - 0.5 alpha[4] * z_p_prime_base_134[0]
   + 0.5 p[0] * z_p_prime_base_134[0] - 0.5 p[1] * z_p_prime_base_134[0]
   - 0.5 p[2] * z_p_prime_base_134[0] + 0.5 p[3] * z_p_prime_base_134[0]
   - 0.5 p[4] * z_p_prime_base_134[0] ] <= 0
 infeasibility_p_prime_base_134_row_1_leq_0: 5e-06 z_p_prime_base_134[1]
   + [ 0.5 alpha[0] * z_p_prime_base_134[1]
   - 0.5 alpha[1] * z_p_prime_base_134[1]
   - 0.5 alpha[2] * z_p_prime_base_134[1]
   + 0.5 alpha[3] * z_p_prime_base_134[1]
   + 0.5 alpha[4] * z_p_prime_base_134[1]
   + 0.5 p[0] * z_p_prime_base_134[1] + 0.5 p[1] * z_p_prime_base_134[1]
   + 0.5 p[2] * z_p_prime_base_134[1] - 0.5 p[3] * z_p_prime_base_134[1]
   + 0.5 p[4] * z_p_prime_base_134[1] ] <= 0
 infeasibility_p_prime_base_134_row_2_leq_0: 5e-06 z_p_prime_base_134[2]
   + [ - alpha[1] * z_p_prime_base_134[2] + p[1] * z_p_prime_base_134[2] ]
   <= 0
 infeasibility_p_prime_base_134_row_3_leq_0: 5e-06 z_p_prime_base_134[3]
   + [ - alpha[2] * z_p_prime_base_134[3] + p[2] * z_p_prime_base_134[3] ]
   <= 0
 infeasibility_p_prime_base_134_row_4_leq_0: 5e-06 z_p_prime_base_134[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_134[4]
   + 0.5 alpha[1] * z_p_prime_base_134[4]
   + 0.5 alpha[2] * z_p_prime_base_134[4]
   - 0.5 alpha[3] * z_p_prime_base_134[4]
   + 0.5 alpha[4] * z_p_prime_base_134[4]
   - 0.5 p[0] * z_p_prime_base_134[4] - 0.5 p[1] * z_p_prime_base_134[4]
   - 0.5 p[2] * z_p_prime_base_134[4] + 0.5 p[3] * z_p_prime_base_134[4]
   + 0.5 p[4] * z_p_prime_base_134[4] ] <= 0
 infeasibility_p_base_135_row_0_leq_0: 5e-06 z_p_base_135[0] + [
   0.5 p[0] * z_p_base_135[0] - 0.5 p[1] * z_p_base_135[0]
   + 0.5 p[2] * z_p_base_135[0] + 0.5 p[3] * z_p_base_135[0]
   - 0.5 p[4] * z_p_base_135[0] ] <= 0
 infeasibility_p_base_135_row_1_leq_0: 5e-06 z_p_base_135[1] + [
   0.5 p[0] * z_p_base_135[1] + 0.5 p[1] * z_p_base_135[1]
   - 0.5 p[2] * z_p_base_135[1] - 0.5 p[3] * z_p_base_135[1]
   + 0.5 p[4] * z_p_base_135[1] ] <= 0
 infeasibility_p_base_135_row_2_leq_0: 5e-06 z_p_base_135[2] + [
   p[1] * z_p_base_135[2] ] <= 0
 infeasibility_p_base_135_row_3_leq_0: 5e-06 z_p_base_135[3] + [
   p[2] * z_p_base_135[3] ] <= 0
 infeasibility_p_base_135_row_4_leq_0: 5e-06 z_p_base_135[4] + [
   - 0.5 p[0] * z_p_base_135[4] - 0.5 p[1] * z_p_base_135[4]
   - 0.5 p[2] * z_p_base_135[4] + 0.5 p[3] * z_p_base_135[4]
   + 0.5 p[4] * z_p_base_135[4] ] <= 0
 infeasibility_p_prime_base_135_row_0_leq_0: 5e-06 z_p_prime_base_135[0]
   + [ 0.5 alpha[0] * z_p_prime_base_135[0]
   + 0.5 alpha[1] * z_p_prime_base_135[0]
   - 0.5 alpha[2] * z_p_prime_base_135[0]
   - 0.5 alpha[3] * z_p_prime_base_135[0]
   - 0.5 alpha[4] * z_p_prime_base_135[0]
   + 0.5 p[0] * z_p_prime_base_135[0] - 0.5 p[1] * z_p_prime_base_135[0]
   + 0.5 p[2] * z_p_prime_base_135[0] + 0.5 p[3] * z_p_prime_base_135[0]
   - 0.5 p[4] * z_p_prime_base_135[0] ] <= 0
 infeasibility_p_prime_base_135_row_1_leq_0: 5e-06 z_p_prime_base_135[1]
   + [ 0.5 alpha[0] * z_p_prime_base_135[1]
   - 0.5 alpha[1] * z_p_prime_base_135[1]
   + 0.5 alpha[2] * z_p_prime_base_135[1]
   + 0.5 alpha[3] * z_p_prime_base_135[1]
   + 0.5 alpha[4] * z_p_prime_base_135[1]
   + 0.5 p[0] * z_p_prime_base_135[1] + 0.5 p[1] * z_p_prime_base_135[1]
   - 0.5 p[2] * z_p_prime_base_135[1] - 0.5 p[3] * z_p_prime_base_135[1]
   + 0.5 p[4] * z_p_prime_base_135[1] ] <= 0
 infeasibility_p_prime_base_135_row_2_leq_0: 5e-06 z_p_prime_base_135[2]
   + [ - alpha[1] * z_p_prime_base_135[2] + p[1] * z_p_prime_base_135[2] ]
   <= 0
 infeasibility_p_prime_base_135_row_3_leq_0: 5e-06 z_p_prime_base_135[3]
   + [ - alpha[2] * z_p_prime_base_135[3] + p[2] * z_p_prime_base_135[3] ]
   <= 0
 infeasibility_p_prime_base_135_row_4_leq_0: 5e-06 z_p_prime_base_135[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_135[4]
   + 0.5 alpha[1] * z_p_prime_base_135[4]
   + 0.5 alpha[2] * z_p_prime_base_135[4]
   - 0.5 alpha[3] * z_p_prime_base_135[4]
   + 0.5 alpha[4] * z_p_prime_base_135[4]
   - 0.5 p[0] * z_p_prime_base_135[4] - 0.5 p[1] * z_p_prime_base_135[4]
   - 0.5 p[2] * z_p_prime_base_135[4] + 0.5 p[3] * z_p_prime_base_135[4]
   + 0.5 p[4] * z_p_prime_base_135[4] ] <= 0
 infeasibility_p_base_136_row_0_leq_0: 5e-06 z_p_base_136[0] + [
   0.5 p[0] * z_p_base_136[0] + 0.5 p[1] * z_p_base_136[0]
   - 0.5 p[2] * z_p_base_136[0] + 0.5 p[3] * z_p_base_136[0]
   - 0.5 p[4] * z_p_base_136[0] ] <= 0
 infeasibility_p_base_136_row_1_leq_0: 5e-06 z_p_base_136[1] + [
   0.5 p[0] * z_p_base_136[1] - 0.5 p[1] * z_p_base_136[1]
   + 0.5 p[2] * z_p_base_136[1] - 0.5 p[3] * z_p_base_136[1]
   + 0.5 p[4] * z_p_base_136[1] ] <= 0
 infeasibility_p_base_136_row_2_leq_0: 5e-06 z_p_base_136[2] + [
   p[1] * z_p_base_136[2] ] <= 0
 infeasibility_p_base_136_row_3_leq_0: 5e-06 z_p_base_136[3] + [
   p[2] * z_p_base_136[3] ] <= 0
 infeasibility_p_base_136_row_4_leq_0: 5e-06 z_p_base_136[4] + [
   - 0.5 p[0] * z_p_base_136[4] - 0.5 p[1] * z_p_base_136[4]
   - 0.5 p[2] * z_p_base_136[4] + 0.5 p[3] * z_p_base_136[4]
   + 0.5 p[4] * z_p_base_136[4] ] <= 0
 infeasibility_p_prime_base_136_row_0_leq_0: 5e-06 z_p_prime_base_136[0]
   + [ 0.5 alpha[0] * z_p_prime_base_136[0]
   - 0.5 alpha[1] * z_p_prime_base_136[0]
   + 0.5 alpha[2] * z_p_prime_base_136[0]
   - 0.5 alpha[3] * z_p_prime_base_136[0]
   - 0.5 alpha[4] * z_p_prime_base_136[0]
   + 0.5 p[0] * z_p_prime_base_136[0] + 0.5 p[1] * z_p_prime_base_136[0]
   - 0.5 p[2] * z_p_prime_base_136[0] + 0.5 p[3] * z_p_prime_base_136[0]
   - 0.5 p[4] * z_p_prime_base_136[0] ] <= 0
 infeasibility_p_prime_base_136_row_1_leq_0: 5e-06 z_p_prime_base_136[1]
   + [ 0.5 alpha[0] * z_p_prime_base_136[1]
   + 0.5 alpha[1] * z_p_prime_base_136[1]
   - 0.5 alpha[2] * z_p_prime_base_136[1]
   + 0.5 alpha[3] * z_p_prime_base_136[1]
   + 0.5 alpha[4] * z_p_prime_base_136[1]
   + 0.5 p[0] * z_p_prime_base_136[1] - 0.5 p[1] * z_p_prime_base_136[1]
   + 0.5 p[2] * z_p_prime_base_136[1] - 0.5 p[3] * z_p_prime_base_136[1]
   + 0.5 p[4] * z_p_prime_base_136[1] ] <= 0
 infeasibility_p_prime_base_136_row_2_leq_0: 5e-06 z_p_prime_base_136[2]
   + [ - alpha[1] * z_p_prime_base_136[2] + p[1] * z_p_prime_base_136[2] ]
   <= 0
 infeasibility_p_prime_base_136_row_3_leq_0: 5e-06 z_p_prime_base_136[3]
   + [ - alpha[2] * z_p_prime_base_136[3] + p[2] * z_p_prime_base_136[3] ]
   <= 0
 infeasibility_p_prime_base_136_row_4_leq_0: 5e-06 z_p_prime_base_136[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_136[4]
   + 0.5 alpha[1] * z_p_prime_base_136[4]
   + 0.5 alpha[2] * z_p_prime_base_136[4]
   - 0.5 alpha[3] * z_p_prime_base_136[4]
   + 0.5 alpha[4] * z_p_prime_base_136[4]
   - 0.5 p[0] * z_p_prime_base_136[4] - 0.5 p[1] * z_p_prime_base_136[4]
   - 0.5 p[2] * z_p_prime_base_136[4] + 0.5 p[3] * z_p_prime_base_136[4]
   + 0.5 p[4] * z_p_prime_base_136[4] ] <= 0
 infeasibility_p_base_137_row_0_leq_0: 5e-06 z_p_base_137[0] + [
   0.5 p[0] * z_p_base_137[0] + 0.5 p[1] * z_p_base_137[0]
   + 0.5 p[2] * z_p_base_137[0] + 0.5 p[3] * z_p_base_137[0]
   - 0.5 p[4] * z_p_base_137[0] ] <= 0
 infeasibility_p_base_137_row_1_leq_0: 5e-06 z_p_base_137[1] + [
   0.5 p[0] * z_p_base_137[1] - 0.5 p[1] * z_p_base_137[1]
   - 0.5 p[2] * z_p_base_137[1] - 0.5 p[3] * z_p_base_137[1]
   + 0.5 p[4] * z_p_base_137[1] ] <= 0
 infeasibility_p_base_137_row_2_leq_0: 5e-06 z_p_base_137[2] + [
   p[1] * z_p_base_137[2] ] <= 0
 infeasibility_p_base_137_row_3_leq_0: 5e-06 z_p_base_137[3] + [
   p[2] * z_p_base_137[3] ] <= 0
 infeasibility_p_base_137_row_4_leq_0: 5e-06 z_p_base_137[4] + [
   - 0.5 p[0] * z_p_base_137[4] - 0.5 p[1] * z_p_base_137[4]
   - 0.5 p[2] * z_p_base_137[4] + 0.5 p[3] * z_p_base_137[4]
   + 0.5 p[4] * z_p_base_137[4] ] <= 0
 infeasibility_p_prime_base_137_row_0_leq_0: 5e-06 z_p_prime_base_137[0]
   + [ 0.5 alpha[0] * z_p_prime_base_137[0]
   - 0.5 alpha[1] * z_p_prime_base_137[0]
   - 0.5 alpha[2] * z_p_prime_base_137[0]
   - 0.5 alpha[3] * z_p_prime_base_137[0]
   - 0.5 alpha[4] * z_p_prime_base_137[0]
   + 0.5 p[0] * z_p_prime_base_137[0] + 0.5 p[1] * z_p_prime_base_137[0]
   + 0.5 p[2] * z_p_prime_base_137[0] + 0.5 p[3] * z_p_prime_base_137[0]
   - 0.5 p[4] * z_p_prime_base_137[0] ] <= 0
 infeasibility_p_prime_base_137_row_1_leq_0: 5e-06 z_p_prime_base_137[1]
   + [ 0.5 alpha[0] * z_p_prime_base_137[1]
   + 0.5 alpha[1] * z_p_prime_base_137[1]
   + 0.5 alpha[2] * z_p_prime_base_137[1]
   + 0.5 alpha[3] * z_p_prime_base_137[1]
   + 0.5 alpha[4] * z_p_prime_base_137[1]
   + 0.5 p[0] * z_p_prime_base_137[1] - 0.5 p[1] * z_p_prime_base_137[1]
   - 0.5 p[2] * z_p_prime_base_137[1] - 0.5 p[3] * z_p_prime_base_137[1]
   + 0.5 p[4] * z_p_prime_base_137[1] ] <= 0
 infeasibility_p_prime_base_137_row_2_leq_0: 5e-06 z_p_prime_base_137[2]
   + [ - alpha[1] * z_p_prime_base_137[2] + p[1] * z_p_prime_base_137[2] ]
   <= 0
 infeasibility_p_prime_base_137_row_3_leq_0: 5e-06 z_p_prime_base_137[3]
   + [ - alpha[2] * z_p_prime_base_137[3] + p[2] * z_p_prime_base_137[3] ]
   <= 0
 infeasibility_p_prime_base_137_row_4_leq_0: 5e-06 z_p_prime_base_137[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_137[4]
   + 0.5 alpha[1] * z_p_prime_base_137[4]
   + 0.5 alpha[2] * z_p_prime_base_137[4]
   - 0.5 alpha[3] * z_p_prime_base_137[4]
   + 0.5 alpha[4] * z_p_prime_base_137[4]
   - 0.5 p[0] * z_p_prime_base_137[4] - 0.5 p[1] * z_p_prime_base_137[4]
   - 0.5 p[2] * z_p_prime_base_137[4] + 0.5 p[3] * z_p_prime_base_137[4]
   + 0.5 p[4] * z_p_prime_base_137[4] ] <= 0
 infeasibility_p_base_138_row_0_leq_0: 5e-06 z_p_base_138[0] + [
   p[0] * z_p_base_138[0] ] <= 0
 infeasibility_p_base_138_row_1_leq_0: 5e-06 z_p_base_138[1] + [
   0.5 p[0] * z_p_base_138[1] + 0.5 p[1] * z_p_base_138[1]
   + 0.5 p[2] * z_p_base_138[1] - 0.5 p[3] * z_p_base_138[1]
   - 0.5 p[4] * z_p_base_138[1] ] <= 0
 infeasibility_p_base_138_row_2_leq_0: 5e-06 z_p_base_138[2] + [
   - 0.5 p[0] * z_p_base_138[2] + 0.5 p[1] * z_p_base_138[2]
   - 0.5 p[2] * z_p_base_138[2] + 0.5 p[3] * z_p_base_138[2]
   - 0.5 p[4] * z_p_base_138[2] ] <= 0
 infeasibility_p_base_138_row_3_leq_0: 5e-06 z_p_base_138[3] + [
   p[4] * z_p_base_138[3] ] <= 0
 infeasibility_p_base_138_row_4_leq_0: 5e-06 z_p_base_138[4] + [
   - 0.5 p[0] * z_p_base_138[4] - 0.5 p[1] * z_p_base_138[4]
   + 0.5 p[2] * z_p_base_138[4] + 0.5 p[3] * z_p_base_138[4]
   + 0.5 p[4] * z_p_base_138[4] ] <= 0
 infeasibility_p_prime_base_138_row_0_leq_0: 5e-06 z_p_prime_base_138[0]
   + [ alpha[0] * z_p_prime_base_138[0] + p[0] * z_p_prime_base_138[0] ]
   <= 0
 infeasibility_p_prime_base_138_row_1_leq_0: 5e-06 z_p_prime_base_138[1]
   + [ 0.5 alpha[0] * z_p_prime_base_138[1]
   - 0.5 alpha[1] * z_p_prime_base_138[1]
   - 0.5 alpha[2] * z_p_prime_base_138[1]
   + 0.5 alpha[3] * z_p_prime_base_138[1]
   - 0.5 alpha[4] * z_p_prime_base_138[1]
   + 0.5 p[0] * z_p_prime_base_138[1] + 0.5 p[1] * z_p_prime_base_138[1]
   + 0.5 p[2] * z_p_prime_base_138[1] - 0.5 p[3] * z_p_prime_base_138[1]
   - 0.5 p[4] * z_p_prime_base_138[1] ] <= 0
 infeasibility_p_prime_base_138_row_2_leq_0: 5e-06 z_p_prime_base_138[2]
   + [ - 0.5 alpha[0] * z_p_prime_base_138[2]
   - 0.5 alpha[1] * z_p_prime_base_138[2]
   + 0.5 alpha[2] * z_p_prime_base_138[2]
   - 0.5 alpha[3] * z_p_prime_base_138[2]
   - 0.5 alpha[4] * z_p_prime_base_138[2]
   - 0.5 p[0] * z_p_prime_base_138[2] + 0.5 p[1] * z_p_prime_base_138[2]
   - 0.5 p[2] * z_p_prime_base_138[2] + 0.5 p[3] * z_p_prime_base_138[2]
   - 0.5 p[4] * z_p_prime_base_138[2] ] <= 0
 infeasibility_p_prime_base_138_row_3_leq_0: 5e-06 z_p_prime_base_138[3]
   + [ alpha[4] * z_p_prime_base_138[3] + p[4] * z_p_prime_base_138[3] ]
   <= 0
 infeasibility_p_prime_base_138_row_4_leq_0: 5e-06 z_p_prime_base_138[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_138[4]
   + 0.5 alpha[1] * z_p_prime_base_138[4]
   - 0.5 alpha[2] * z_p_prime_base_138[4]
   - 0.5 alpha[3] * z_p_prime_base_138[4]
   + 0.5 alpha[4] * z_p_prime_base_138[4]
   - 0.5 p[0] * z_p_prime_base_138[4] - 0.5 p[1] * z_p_prime_base_138[4]
   + 0.5 p[2] * z_p_prime_base_138[4] + 0.5 p[3] * z_p_prime_base_138[4]
   + 0.5 p[4] * z_p_prime_base_138[4] ] <= 0
 infeasibility_p_base_139_row_0_leq_0: 5e-06 z_p_base_139[0] + [
   p[0] * z_p_base_139[0] ] <= 0
 infeasibility_p_base_139_row_1_leq_0: 5e-06 z_p_base_139[1] + [
   0.5 p[0] * z_p_base_139[1] + 0.5 p[1] * z_p_base_139[1]
   + 0.5 p[2] * z_p_base_139[1] - 0.5 p[3] * z_p_base_139[1]
   - 0.5 p[4] * z_p_base_139[1] ] <= 0
 infeasibility_p_base_139_row_2_leq_0: 5e-06 z_p_base_139[2] + [
   - p[0] * z_p_base_139[2] + p[3] * z_p_base_139[2] ] <= 0
 infeasibility_p_base_139_row_3_leq_0: 5e-06 z_p_base_139[3] + [
   0.5 p[0] * z_p_base_139[3] + 0.5 p[1] * z_p_base_139[3]
   - 0.5 p[2] * z_p_base_139[3] - 0.5 p[3] * z_p_base_139[3]
   + 0.5 p[4] * z_p_base_139[3] ] <= 0
 infeasibility_p_base_139_row_4_leq_0: 5e-06 z_p_base_139[4] + [
   - 0.5 p[0] * z_p_base_139[4] - 0.5 p[1] * z_p_base_139[4]
   + 0.5 p[2] * z_p_base_139[4] + 0.5 p[3] * z_p_base_139[4]
   + 0.5 p[4] * z_p_base_139[4] ] <= 0
 infeasibility_p_prime_base_139_row_0_leq_0: 5e-06 z_p_prime_base_139[0]
   + [ alpha[0] * z_p_prime_base_139[0] + p[0] * z_p_prime_base_139[0] ]
   <= 0
 infeasibility_p_prime_base_139_row_1_leq_0: 5e-06 z_p_prime_base_139[1]
   + [ 0.5 alpha[0] * z_p_prime_base_139[1]
   - 0.5 alpha[1] * z_p_prime_base_139[1]
   - 0.5 alpha[2] * z_p_prime_base_139[1]
   + 0.5 alpha[3] * z_p_prime_base_139[1]
   - 0.5 alpha[4] * z_p_prime_base_139[1]
   + 0.5 p[0] * z_p_prime_base_139[1] + 0.5 p[1] * z_p_prime_base_139[1]
   + 0.5 p[2] * z_p_prime_base_139[1] - 0.5 p[3] * z_p_prime_base_139[1]
   - 0.5 p[4] * z_p_prime_base_139[1] ] <= 0
 infeasibility_p_prime_base_139_row_2_leq_0: 5e-06 z_p_prime_base_139[2]
   + [ - alpha[0] * z_p_prime_base_139[2]
   - alpha[3] * z_p_prime_base_139[2] - p[0] * z_p_prime_base_139[2]
   + p[3] * z_p_prime_base_139[2] ] <= 0
 infeasibility_p_prime_base_139_row_3_leq_0: 5e-06 z_p_prime_base_139[3]
   + [ 0.5 alpha[0] * z_p_prime_base_139[3]
   - 0.5 alpha[1] * z_p_prime_base_139[3]
   + 0.5 alpha[2] * z_p_prime_base_139[3]
   + 0.5 alpha[3] * z_p_prime_base_139[3]
   + 0.5 alpha[4] * z_p_prime_base_139[3]
   + 0.5 p[0] * z_p_prime_base_139[3] + 0.5 p[1] * z_p_prime_base_139[3]
   - 0.5 p[2] * z_p_prime_base_139[3] - 0.5 p[3] * z_p_prime_base_139[3]
   + 0.5 p[4] * z_p_prime_base_139[3] ] <= 0
 infeasibility_p_prime_base_139_row_4_leq_0: 5e-06 z_p_prime_base_139[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_139[4]
   + 0.5 alpha[1] * z_p_prime_base_139[4]
   - 0.5 alpha[2] * z_p_prime_base_139[4]
   - 0.5 alpha[3] * z_p_prime_base_139[4]
   + 0.5 alpha[4] * z_p_prime_base_139[4]
   - 0.5 p[0] * z_p_prime_base_139[4] - 0.5 p[1] * z_p_prime_base_139[4]
   + 0.5 p[2] * z_p_prime_base_139[4] + 0.5 p[3] * z_p_prime_base_139[4]
   + 0.5 p[4] * z_p_prime_base_139[4] ] <= 0
 infeasibility_p_base_140_row_0_leq_0: 5e-06 z_p_base_140[0] + [
   p[0] * z_p_base_140[0] ] <= 0
 infeasibility_p_base_140_row_1_leq_0: 5e-06 z_p_base_140[1] + [
   p[2] * z_p_base_140[1] ] <= 0
 infeasibility_p_base_140_row_2_leq_0: 5e-06 z_p_base_140[2] + [
   - 0.5 p[0] * z_p_base_140[2] + 0.5 p[1] * z_p_base_140[2]
   - 0.5 p[2] * z_p_base_140[2] + 0.5 p[3] * z_p_base_140[2]
   - 0.5 p[4] * z_p_base_140[2] ] <= 0
 infeasibility_p_base_140_row_3_leq_0: 5e-06 z_p_base_140[3] + [
   0.5 p[0] * z_p_base_140[3] + 0.5 p[1] * z_p_base_140[3]
   - 0.5 p[2] * z_p_base_140[3] - 0.5 p[3] * z_p_base_140[3]
   + 0.5 p[4] * z_p_base_140[3] ] <= 0
 infeasibility_p_base_140_row_4_leq_0: 5e-06 z_p_base_140[4] + [
   - 0.5 p[0] * z_p_base_140[4] - 0.5 p[1] * z_p_base_140[4]
   + 0.5 p[2] * z_p_base_140[4] + 0.5 p[3] * z_p_base_140[4]
   + 0.5 p[4] * z_p_base_140[4] ] <= 0
 infeasibility_p_prime_base_140_row_0_leq_0: 5e-06 z_p_prime_base_140[0]
   + [ alpha[0] * z_p_prime_base_140[0] + p[0] * z_p_prime_base_140[0] ]
   <= 0
 infeasibility_p_prime_base_140_row_1_leq_0: 5e-06 z_p_prime_base_140[1]
   + [ - alpha[2] * z_p_prime_base_140[1] + p[2] * z_p_prime_base_140[1] ]
   <= 0
 infeasibility_p_prime_base_140_row_2_leq_0: 5e-06 z_p_prime_base_140[2]
   + [ - 0.5 alpha[0] * z_p_prime_base_140[2]
   - 0.5 alpha[1] * z_p_prime_base_140[2]
   + 0.5 alpha[2] * z_p_prime_base_140[2]
   - 0.5 alpha[3] * z_p_prime_base_140[2]
   - 0.5 alpha[4] * z_p_prime_base_140[2]
   - 0.5 p[0] * z_p_prime_base_140[2] + 0.5 p[1] * z_p_prime_base_140[2]
   - 0.5 p[2] * z_p_prime_base_140[2] + 0.5 p[3] * z_p_prime_base_140[2]
   - 0.5 p[4] * z_p_prime_base_140[2] ] <= 0
 infeasibility_p_prime_base_140_row_3_leq_0: 5e-06 z_p_prime_base_140[3]
   + [ 0.5 alpha[0] * z_p_prime_base_140[3]
   - 0.5 alpha[1] * z_p_prime_base_140[3]
   + 0.5 alpha[2] * z_p_prime_base_140[3]
   + 0.5 alpha[3] * z_p_prime_base_140[3]
   + 0.5 alpha[4] * z_p_prime_base_140[3]
   + 0.5 p[0] * z_p_prime_base_140[3] + 0.5 p[1] * z_p_prime_base_140[3]
   - 0.5 p[2] * z_p_prime_base_140[3] - 0.5 p[3] * z_p_prime_base_140[3]
   + 0.5 p[4] * z_p_prime_base_140[3] ] <= 0
 infeasibility_p_prime_base_140_row_4_leq_0: 5e-06 z_p_prime_base_140[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_140[4]
   + 0.5 alpha[1] * z_p_prime_base_140[4]
   - 0.5 alpha[2] * z_p_prime_base_140[4]
   - 0.5 alpha[3] * z_p_prime_base_140[4]
   + 0.5 alpha[4] * z_p_prime_base_140[4]
   - 0.5 p[0] * z_p_prime_base_140[4] - 0.5 p[1] * z_p_prime_base_140[4]
   + 0.5 p[2] * z_p_prime_base_140[4] + 0.5 p[3] * z_p_prime_base_140[4]
   + 0.5 p[4] * z_p_prime_base_140[4] ] <= 0
 infeasibility_p_base_141_row_0_leq_0: 5e-06 z_p_base_141[0] + [
   p[0] * z_p_base_141[0] ] <= 0
 infeasibility_p_base_141_row_1_leq_0: 5e-06 z_p_base_141[1] + [
   0.5 p[0] * z_p_base_141[1] + 0.5 p[1] * z_p_base_141[1]
   + 0.5 p[2] * z_p_base_141[1] - 0.5 p[3] * z_p_base_141[1]
   - 0.5 p[4] * z_p_base_141[1] ] <= 0
 infeasibility_p_base_141_row_2_leq_0: 5e-06 z_p_base_141[2] + [
   - 0.5 p[0] * z_p_base_141[2] + 0.5 p[1] * z_p_base_141[2]
   - 0.5 p[2] * z_p_base_141[2] + 0.5 p[3] * z_p_base_141[2]
   + 0.5 p[4] * z_p_base_141[2] ] <= 0
 infeasibility_p_base_141_row_3_leq_0: 5e-06 z_p_base_141[3] + [
   - 0.5 p[0] * z_p_base_141[3] - 0.5 p[1] * z_p_base_141[3]
   + 0.5 p[2] * z_p_base_141[3] + 0.5 p[3] * z_p_base_141[3]
   - 0.5 p[4] * z_p_base_141[3] ] <= 0
 infeasibility_p_base_141_row_4_leq_0: 5e-06 z_p_base_141[4] + [
   p[4] * z_p_base_141[4] ] <= 0
 infeasibility_p_prime_base_141_row_0_leq_0: 5e-06 z_p_prime_base_141[0]
   + [ alpha[0] * z_p_prime_base_141[0] + p[0] * z_p_prime_base_141[0] ]
   <= 0
 infeasibility_p_prime_base_141_row_1_leq_0: 5e-06 z_p_prime_base_141[1]
   + [ 0.5 alpha[0] * z_p_prime_base_141[1]
   - 0.5 alpha[1] * z_p_prime_base_141[1]
   - 0.5 alpha[2] * z_p_prime_base_141[1]
   + 0.5 alpha[3] * z_p_prime_base_141[1]
   - 0.5 alpha[4] * z_p_prime_base_141[1]
   + 0.5 p[0] * z_p_prime_base_141[1] + 0.5 p[1] * z_p_prime_base_141[1]
   + 0.5 p[2] * z_p_prime_base_141[1] - 0.5 p[3] * z_p_prime_base_141[1]
   - 0.5 p[4] * z_p_prime_base_141[1] ] <= 0
 infeasibility_p_prime_base_141_row_2_leq_0: 5e-06 z_p_prime_base_141[2]
   + [ - 0.5 alpha[0] * z_p_prime_base_141[2]
   - 0.5 alpha[1] * z_p_prime_base_141[2]
   + 0.5 alpha[2] * z_p_prime_base_141[2]
   - 0.5 alpha[3] * z_p_prime_base_141[2]
   + 0.5 alpha[4] * z_p_prime_base_141[2]
   - 0.5 p[0] * z_p_prime_base_141[2] + 0.5 p[1] * z_p_prime_base_141[2]
   - 0.5 p[2] * z_p_prime_base_141[2] + 0.5 p[3] * z_p_prime_base_141[2]
   + 0.5 p[4] * z_p_prime_base_141[2] ] <= 0
 infeasibility_p_prime_base_141_row_3_leq_0: 5e-06 z_p_prime_base_141[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_141[3]
   + 0.5 alpha[1] * z_p_prime_base_141[3]
   - 0.5 alpha[2] * z_p_prime_base_141[3]
   - 0.5 alpha[3] * z_p_prime_base_141[3]
   - 0.5 alpha[4] * z_p_prime_base_141[3]
   - 0.5 p[0] * z_p_prime_base_141[3] - 0.5 p[1] * z_p_prime_base_141[3]
   + 0.5 p[2] * z_p_prime_base_141[3] + 0.5 p[3] * z_p_prime_base_141[3]
   - 0.5 p[4] * z_p_prime_base_141[3] ] <= 0
 infeasibility_p_prime_base_141_row_4_leq_0: 5e-06 z_p_prime_base_141[4]
   + [ alpha[4] * z_p_prime_base_141[4] + p[4] * z_p_prime_base_141[4] ]
   <= 0
 infeasibility_p_base_142_row_0_leq_0: 5e-06 z_p_base_142[0] + [
   p[0] * z_p_base_142[0] ] <= 0
 infeasibility_p_base_142_row_1_leq_0: 5e-06 z_p_base_142[1] + [
   0.5 p[0] * z_p_base_142[1] + 0.5 p[1] * z_p_base_142[1]
   + 0.5 p[2] * z_p_base_142[1] - 0.5 p[3] * z_p_base_142[1]
   + 0.5 p[4] * z_p_base_142[1] ] <= 0
 infeasibility_p_base_142_row_2_leq_0: 5e-06 z_p_base_142[2] + [
   - 0.5 p[0] * z_p_base_142[2] + 0.5 p[1] * z_p_base_142[2]
   - 0.5 p[2] * z_p_base_142[2] + 0.5 p[3] * z_p_base_142[2]
   - 0.5 p[4] * z_p_base_142[2] ] <= 0
 infeasibility_p_base_142_row_3_leq_0: 5e-06 z_p_base_142[3] + [
   - 0.5 p[0] * z_p_base_142[3] - 0.5 p[1] * z_p_base_142[3]
   + 0.5 p[2] * z_p_base_142[3] + 0.5 p[3] * z_p_base_142[3]
   - 0.5 p[4] * z_p_base_142[3] ] <= 0
 infeasibility_p_base_142_row_4_leq_0: 5e-06 z_p_base_142[4] + [
   p[4] * z_p_base_142[4] ] <= 0
 infeasibility_p_prime_base_142_row_0_leq_0: 5e-06 z_p_prime_base_142[0]
   + [ alpha[0] * z_p_prime_base_142[0] + p[0] * z_p_prime_base_142[0] ]
   <= 0
 infeasibility_p_prime_base_142_row_1_leq_0: 5e-06 z_p_prime_base_142[1]
   + [ 0.5 alpha[0] * z_p_prime_base_142[1]
   - 0.5 alpha[1] * z_p_prime_base_142[1]
   - 0.5 alpha[2] * z_p_prime_base_142[1]
   + 0.5 alpha[3] * z_p_prime_base_142[1]
   + 0.5 alpha[4] * z_p_prime_base_142[1]
   + 0.5 p[0] * z_p_prime_base_142[1] + 0.5 p[1] * z_p_prime_base_142[1]
   + 0.5 p[2] * z_p_prime_base_142[1] - 0.5 p[3] * z_p_prime_base_142[1]
   + 0.5 p[4] * z_p_prime_base_142[1] ] <= 0
 infeasibility_p_prime_base_142_row_2_leq_0: 5e-06 z_p_prime_base_142[2]
   + [ - 0.5 alpha[0] * z_p_prime_base_142[2]
   - 0.5 alpha[1] * z_p_prime_base_142[2]
   + 0.5 alpha[2] * z_p_prime_base_142[2]
   - 0.5 alpha[3] * z_p_prime_base_142[2]
   - 0.5 alpha[4] * z_p_prime_base_142[2]
   - 0.5 p[0] * z_p_prime_base_142[2] + 0.5 p[1] * z_p_prime_base_142[2]
   - 0.5 p[2] * z_p_prime_base_142[2] + 0.5 p[3] * z_p_prime_base_142[2]
   - 0.5 p[4] * z_p_prime_base_142[2] ] <= 0
 infeasibility_p_prime_base_142_row_3_leq_0: 5e-06 z_p_prime_base_142[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_142[3]
   + 0.5 alpha[1] * z_p_prime_base_142[3]
   - 0.5 alpha[2] * z_p_prime_base_142[3]
   - 0.5 alpha[3] * z_p_prime_base_142[3]
   - 0.5 alpha[4] * z_p_prime_base_142[3]
   - 0.5 p[0] * z_p_prime_base_142[3] - 0.5 p[1] * z_p_prime_base_142[3]
   + 0.5 p[2] * z_p_prime_base_142[3] + 0.5 p[3] * z_p_prime_base_142[3]
   - 0.5 p[4] * z_p_prime_base_142[3] ] <= 0
 infeasibility_p_prime_base_142_row_4_leq_0: 5e-06 z_p_prime_base_142[4]
   + [ alpha[4] * z_p_prime_base_142[4] + p[4] * z_p_prime_base_142[4] ]
   <= 0
 infeasibility_p_base_143_row_0_leq_0: 5e-06 z_p_base_143[0] + [
   p[0] * z_p_base_143[0] ] <= 0
 infeasibility_p_base_143_row_1_leq_0: 5e-06 z_p_base_143[1] + [
   0.5 p[0] * z_p_base_143[1] + 0.5 p[1] * z_p_base_143[1]
   + 0.5 p[2] * z_p_base_143[1] - 0.5 p[3] * z_p_base_143[1]
   - 0.5 p[4] * z_p_base_143[1] ] <= 0
 infeasibility_p_base_143_row_2_leq_0: 5e-06 z_p_base_143[2] + [
   - 0.5 p[0] * z_p_base_143[2] + 0.5 p[1] * z_p_base_143[2]
   - 0.5 p[2] * z_p_base_143[2] + 0.5 p[3] * z_p_base_143[2]
   + 0.5 p[4] * z_p_base_143[2] ] <= 0
 infeasibility_p_base_143_row_3_leq_0: 5e-06 z_p_base_143[3] + [
   - p[0] * z_p_base_143[3] + p[3] * z_p_base_143[3] ] <= 0
 infeasibility_p_base_143_row_4_leq_0: 5e-06 z_p_base_143[4] + [
   0.5 p[0] * z_p_base_143[4] - 0.5 p[1] * z_p_base_143[4]
   + 0.5 p[2] * z_p_base_143[4] - 0.5 p[3] * z_p_base_143[4]
   + 0.5 p[4] * z_p_base_143[4] ] <= 0
 infeasibility_p_prime_base_143_row_0_leq_0: 5e-06 z_p_prime_base_143[0]
   + [ alpha[0] * z_p_prime_base_143[0] + p[0] * z_p_prime_base_143[0] ]
   <= 0
 infeasibility_p_prime_base_143_row_1_leq_0: 5e-06 z_p_prime_base_143[1]
   + [ 0.5 alpha[0] * z_p_prime_base_143[1]
   - 0.5 alpha[1] * z_p_prime_base_143[1]
   - 0.5 alpha[2] * z_p_prime_base_143[1]
   + 0.5 alpha[3] * z_p_prime_base_143[1]
   - 0.5 alpha[4] * z_p_prime_base_143[1]
   + 0.5 p[0] * z_p_prime_base_143[1] + 0.5 p[1] * z_p_prime_base_143[1]
   + 0.5 p[2] * z_p_prime_base_143[1] - 0.5 p[3] * z_p_prime_base_143[1]
   - 0.5 p[4] * z_p_prime_base_143[1] ] <= 0
 infeasibility_p_prime_base_143_row_2_leq_0: 5e-06 z_p_prime_base_143[2]
   + [ - 0.5 alpha[0] * z_p_prime_base_143[2]
   - 0.5 alpha[1] * z_p_prime_base_143[2]
   + 0.5 alpha[2] * z_p_prime_base_143[2]
   - 0.5 alpha[3] * z_p_prime_base_143[2]
   + 0.5 alpha[4] * z_p_prime_base_143[2]
   - 0.5 p[0] * z_p_prime_base_143[2] + 0.5 p[1] * z_p_prime_base_143[2]
   - 0.5 p[2] * z_p_prime_base_143[2] + 0.5 p[3] * z_p_prime_base_143[2]
   + 0.5 p[4] * z_p_prime_base_143[2] ] <= 0
 infeasibility_p_prime_base_143_row_3_leq_0: 5e-06 z_p_prime_base_143[3]
   + [ - alpha[0] * z_p_prime_base_143[3]
   - alpha[3] * z_p_prime_base_143[3] - p[0] * z_p_prime_base_143[3]
   + p[3] * z_p_prime_base_143[3] ] <= 0
 infeasibility_p_prime_base_143_row_4_leq_0: 5e-06 z_p_prime_base_143[4]
   + [ 0.5 alpha[0] * z_p_prime_base_143[4]
   + 0.5 alpha[1] * z_p_prime_base_143[4]
   - 0.5 alpha[2] * z_p_prime_base_143[4]
   + 0.5 alpha[3] * z_p_prime_base_143[4]
   + 0.5 alpha[4] * z_p_prime_base_143[4]
   + 0.5 p[0] * z_p_prime_base_143[4] - 0.5 p[1] * z_p_prime_base_143[4]
   + 0.5 p[2] * z_p_prime_base_143[4] - 0.5 p[3] * z_p_prime_base_143[4]
   + 0.5 p[4] * z_p_prime_base_143[4] ] <= 0
 infeasibility_p_base_144_row_0_leq_0: 5e-06 z_p_base_144[0] + [
   p[0] * z_p_base_144[0] ] <= 0
 infeasibility_p_base_144_row_1_leq_0: 5e-06 z_p_base_144[1] + [
   - 0.5 p[0] * z_p_base_144[1] + 0.5 p[1] * z_p_base_144[1]
   + 0.5 p[2] * z_p_base_144[1] + 0.5 p[3] * z_p_base_144[1]
   - 0.5 p[4] * z_p_base_144[1] ] <= 0
 infeasibility_p_base_144_row_2_leq_0: 5e-06 z_p_base_144[2] + [
   0.5 p[0] * z_p_base_144[2] + 0.5 p[1] * z_p_base_144[2]
   - 0.5 p[2] * z_p_base_144[2] - 0.5 p[3] * z_p_base_144[2]
   + 0.5 p[4] * z_p_base_144[2] ] <= 0
 infeasibility_p_base_144_row_3_leq_0: 5e-06 z_p_base_144[3] + [
   0.5 p[0] * z_p_base_144[3] - 0.5 p[1] * z_p_base_144[3]
   + 0.5 p[2] * z_p_base_144[3] - 0.5 p[3] * z_p_base_144[3]
   + 0.5 p[4] * z_p_base_144[3] ] <= 0
 infeasibility_p_base_144_row_4_leq_0: 5e-06 z_p_base_144[4] + [
   - p[0] * z_p_base_144[4] + p[3] * z_p_base_144[4] ] <= 0
 infeasibility_p_prime_base_144_row_0_leq_0: 5e-06 z_p_prime_base_144[0]
   + [ alpha[0] * z_p_prime_base_144[0] + p[0] * z_p_prime_base_144[0] ]
   <= 0
 infeasibility_p_prime_base_144_row_1_leq_0: 5e-06 z_p_prime_base_144[1]
   + [ - 0.5 alpha[0] * z_p_prime_base_144[1]
   - 0.5 alpha[1] * z_p_prime_base_144[1]
   - 0.5 alpha[2] * z_p_prime_base_144[1]
   - 0.5 alpha[3] * z_p_prime_base_144[1]
   - 0.5 alpha[4] * z_p_prime_base_144[1]
   - 0.5 p[0] * z_p_prime_base_144[1] + 0.5 p[1] * z_p_prime_base_144[1]
   + 0.5 p[2] * z_p_prime_base_144[1] + 0.5 p[3] * z_p_prime_base_144[1]
   - 0.5 p[4] * z_p_prime_base_144[1] ] <= 0
 infeasibility_p_prime_base_144_row_2_leq_0: 5e-06 z_p_prime_base_144[2]
   + [ 0.5 alpha[0] * z_p_prime_base_144[2]
   - 0.5 alpha[1] * z_p_prime_base_144[2]
   + 0.5 alpha[2] * z_p_prime_base_144[2]
   + 0.5 alpha[3] * z_p_prime_base_144[2]
   + 0.5 alpha[4] * z_p_prime_base_144[2]
   + 0.5 p[0] * z_p_prime_base_144[2] + 0.5 p[1] * z_p_prime_base_144[2]
   - 0.5 p[2] * z_p_prime_base_144[2] - 0.5 p[3] * z_p_prime_base_144[2]
   + 0.5 p[4] * z_p_prime_base_144[2] ] <= 0
 infeasibility_p_prime_base_144_row_3_leq_0: 5e-06 z_p_prime_base_144[3]
   + [ 0.5 alpha[0] * z_p_prime_base_144[3]
   + 0.5 alpha[1] * z_p_prime_base_144[3]
   - 0.5 alpha[2] * z_p_prime_base_144[3]
   + 0.5 alpha[3] * z_p_prime_base_144[3]
   + 0.5 alpha[4] * z_p_prime_base_144[3]
   + 0.5 p[0] * z_p_prime_base_144[3] - 0.5 p[1] * z_p_prime_base_144[3]
   + 0.5 p[2] * z_p_prime_base_144[3] - 0.5 p[3] * z_p_prime_base_144[3]
   + 0.5 p[4] * z_p_prime_base_144[3] ] <= 0
 infeasibility_p_prime_base_144_row_4_leq_0: 5e-06 z_p_prime_base_144[4]
   + [ - alpha[0] * z_p_prime_base_144[4]
   - alpha[3] * z_p_prime_base_144[4] - p[0] * z_p_prime_base_144[4]
   + p[3] * z_p_prime_base_144[4] ] <= 0
 infeasibility_p_base_145_row_0_leq_0: 5e-06 z_p_base_145[0] + [
   p[0] * z_p_base_145[0] ] <= 0
 infeasibility_p_base_145_row_1_leq_0: 5e-06 z_p_base_145[1] + [
   p[1] * z_p_base_145[1] ] <= 0
 infeasibility_p_base_145_row_2_leq_0: 5e-06 z_p_base_145[2] + [
   - 0.5 p[0] * z_p_base_145[2] - 0.5 p[1] * z_p_base_145[2]
   + 0.5 p[2] * z_p_base_145[2] + 0.5 p[3] * z_p_base_145[2]
   - 0.5 p[4] * z_p_base_145[2] ] <= 0
 infeasibility_p_base_145_row_3_leq_0: 5e-06 z_p_base_145[3] + [
   0.5 p[0] * z_p_base_145[3] - 0.5 p[1] * z_p_base_145[3]
   + 0.5 p[2] * z_p_base_145[3] - 0.5 p[3] * z_p_base_145[3]
   + 0.5 p[4] * z_p_base_145[3] ] <= 0
 infeasibility_p_base_145_row_4_leq_0: 5e-06 z_p_base_145[4] + [
   - 0.5 p[0] * z_p_base_145[4] + 0.5 p[1] * z_p_base_145[4]
   - 0.5 p[2] * z_p_base_145[4] + 0.5 p[3] * z_p_base_145[4]
   + 0.5 p[4] * z_p_base_145[4] ] <= 0
 infeasibility_p_prime_base_145_row_0_leq_0: 5e-06 z_p_prime_base_145[0]
   + [ alpha[0] * z_p_prime_base_145[0] + p[0] * z_p_prime_base_145[0] ]
   <= 0
 infeasibility_p_prime_base_145_row_1_leq_0: 5e-06 z_p_prime_base_145[1]
   + [ - alpha[1] * z_p_prime_base_145[1] + p[1] * z_p_prime_base_145[1] ]
   <= 0
 infeasibility_p_prime_base_145_row_2_leq_0: 5e-06 z_p_prime_base_145[2]
   + [ - 0.5 alpha[0] * z_p_prime_base_145[2]
   + 0.5 alpha[1] * z_p_prime_base_145[2]
   - 0.5 alpha[2] * z_p_prime_base_145[2]
   - 0.5 alpha[3] * z_p_prime_base_145[2]
   - 0.5 alpha[4] * z_p_prime_base_145[2]
   - 0.5 p[0] * z_p_prime_base_145[2] - 0.5 p[1] * z_p_prime_base_145[2]
   + 0.5 p[2] * z_p_prime_base_145[2] + 0.5 p[3] * z_p_prime_base_145[2]
   - 0.5 p[4] * z_p_prime_base_145[2] ] <= 0
 infeasibility_p_prime_base_145_row_3_leq_0: 5e-06 z_p_prime_base_145[3]
   + [ 0.5 alpha[0] * z_p_prime_base_145[3]
   + 0.5 alpha[1] * z_p_prime_base_145[3]
   - 0.5 alpha[2] * z_p_prime_base_145[3]
   + 0.5 alpha[3] * z_p_prime_base_145[3]
   + 0.5 alpha[4] * z_p_prime_base_145[3]
   + 0.5 p[0] * z_p_prime_base_145[3] - 0.5 p[1] * z_p_prime_base_145[3]
   + 0.5 p[2] * z_p_prime_base_145[3] - 0.5 p[3] * z_p_prime_base_145[3]
   + 0.5 p[4] * z_p_prime_base_145[3] ] <= 0
 infeasibility_p_prime_base_145_row_4_leq_0: 5e-06 z_p_prime_base_145[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_145[4]
   - 0.5 alpha[1] * z_p_prime_base_145[4]
   + 0.5 alpha[2] * z_p_prime_base_145[4]
   - 0.5 alpha[3] * z_p_prime_base_145[4]
   + 0.5 alpha[4] * z_p_prime_base_145[4]
   - 0.5 p[0] * z_p_prime_base_145[4] + 0.5 p[1] * z_p_prime_base_145[4]
   - 0.5 p[2] * z_p_prime_base_145[4] + 0.5 p[3] * z_p_prime_base_145[4]
   + 0.5 p[4] * z_p_prime_base_145[4] ] <= 0
 infeasibility_p_base_146_row_0_leq_0: 5e-06 z_p_base_146[0] + [
   p[0] * z_p_base_146[0] ] <= 0
 infeasibility_p_base_146_row_1_leq_0: 5e-06 z_p_base_146[1] + [
   - 0.5 p[0] * z_p_base_146[1] + 0.5 p[1] * z_p_base_146[1]
   - 0.5 p[2] * z_p_base_146[1] + 0.5 p[3] * z_p_base_146[1]
   - 0.5 p[4] * z_p_base_146[1] ] <= 0
 infeasibility_p_base_146_row_2_leq_0: 5e-06 z_p_base_146[2] + [
   0.5 p[0] * z_p_base_146[2] + 0.5 p[1] * z_p_base_146[2]
   + 0.5 p[2] * z_p_base_146[2] - 0.5 p[3] * z_p_base_146[2]
   + 0.5 p[4] * z_p_base_146[2] ] <= 0
 infeasibility_p_base_146_row_3_leq_0: 5e-06 z_p_base_146[3] + [
   p[2] * z_p_base_146[3] ] <= 0
 infeasibility_p_base_146_row_4_leq_0: 5e-06 z_p_base_146[4] + [
   - 0.5 p[0] * z_p_base_146[4] - 0.5 p[1] * z_p_base_146[4]
   - 0.5 p[2] * z_p_base_146[4] + 0.5 p[3] * z_p_base_146[4]
   + 0.5 p[4] * z_p_base_146[4] ] <= 0
 infeasibility_p_prime_base_146_row_0_leq_0: 5e-06 z_p_prime_base_146[0]
   + [ alpha[0] * z_p_prime_base_146[0] + p[0] * z_p_prime_base_146[0] ]
   <= 0
 infeasibility_p_prime_base_146_row_1_leq_0: 5e-06 z_p_prime_base_146[1]
   + [ - 0.5 alpha[0] * z_p_prime_base_146[1]
   - 0.5 alpha[1] * z_p_prime_base_146[1]
   + 0.5 alpha[2] * z_p_prime_base_146[1]
   - 0.5 alpha[3] * z_p_prime_base_146[1]
   - 0.5 alpha[4] * z_p_prime_base_146[1]
   - 0.5 p[0] * z_p_prime_base_146[1] + 0.5 p[1] * z_p_prime_base_146[1]
   - 0.5 p[2] * z_p_prime_base_146[1] + 0.5 p[3] * z_p_prime_base_146[1]
   - 0.5 p[4] * z_p_prime_base_146[1] ] <= 0
 infeasibility_p_prime_base_146_row_2_leq_0: 5e-06 z_p_prime_base_146[2]
   + [ 0.5 alpha[0] * z_p_prime_base_146[2]
   - 0.5 alpha[1] * z_p_prime_base_146[2]
   - 0.5 alpha[2] * z_p_prime_base_146[2]
   + 0.5 alpha[3] * z_p_prime_base_146[2]
   + 0.5 alpha[4] * z_p_prime_base_146[2]
   + 0.5 p[0] * z_p_prime_base_146[2] + 0.5 p[1] * z_p_prime_base_146[2]
   + 0.5 p[2] * z_p_prime_base_146[2] - 0.5 p[3] * z_p_prime_base_146[2]
   + 0.5 p[4] * z_p_prime_base_146[2] ] <= 0
 infeasibility_p_prime_base_146_row_3_leq_0: 5e-06 z_p_prime_base_146[3]
   + [ - alpha[2] * z_p_prime_base_146[3] + p[2] * z_p_prime_base_146[3] ]
   <= 0
 infeasibility_p_prime_base_146_row_4_leq_0: 5e-06 z_p_prime_base_146[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_146[4]
   + 0.5 alpha[1] * z_p_prime_base_146[4]
   + 0.5 alpha[2] * z_p_prime_base_146[4]
   - 0.5 alpha[3] * z_p_prime_base_146[4]
   + 0.5 alpha[4] * z_p_prime_base_146[4]
   - 0.5 p[0] * z_p_prime_base_146[4] - 0.5 p[1] * z_p_prime_base_146[4]
   - 0.5 p[2] * z_p_prime_base_146[4] + 0.5 p[3] * z_p_prime_base_146[4]
   + 0.5 p[4] * z_p_prime_base_146[4] ] <= 0
 infeasibility_p_base_147_row_0_leq_0: 5e-06 z_p_base_147[0] + [
   p[0] * z_p_base_147[0] ] <= 0
 infeasibility_p_base_147_row_1_leq_0: 5e-06 z_p_base_147[1] + [
   - 0.5 p[0] * z_p_base_147[1] + 0.5 p[1] * z_p_base_147[1]
   + 0.5 p[2] * z_p_base_147[1] + 0.5 p[3] * z_p_base_147[1]
   - 0.5 p[4] * z_p_base_147[1] ] <= 0
 infeasibility_p_base_147_row_2_leq_0: 5e-06 z_p_base_147[2] + [
   0.5 p[0] * z_p_base_147[2] + 0.5 p[1] * z_p_base_147[2]
   - 0.5 p[2] * z_p_base_147[2] - 0.5 p[3] * z_p_base_147[2]
   + 0.5 p[4] * z_p_base_147[2] ] <= 0
 infeasibility_p_base_147_row_3_leq_0: 5e-06 z_p_base_147[3] + [
   p[2] * z_p_base_147[3] ] <= 0
 infeasibility_p_base_147_row_4_leq_0: 5e-06 z_p_base_147[4] + [
   - 0.5 p[0] * z_p_base_147[4] - 0.5 p[1] * z_p_base_147[4]
   - 0.5 p[2] * z_p_base_147[4] + 0.5 p[3] * z_p_base_147[4]
   + 0.5 p[4] * z_p_base_147[4] ] <= 0
 infeasibility_p_prime_base_147_row_0_leq_0: 5e-06 z_p_prime_base_147[0]
   + [ alpha[0] * z_p_prime_base_147[0] + p[0] * z_p_prime_base_147[0] ]
   <= 0
 infeasibility_p_prime_base_147_row_1_leq_0: 5e-06 z_p_prime_base_147[1]
   + [ - 0.5 alpha[0] * z_p_prime_base_147[1]
   - 0.5 alpha[1] * z_p_prime_base_147[1]
   - 0.5 alpha[2] * z_p_prime_base_147[1]
   - 0.5 alpha[3] * z_p_prime_base_147[1]
   - 0.5 alpha[4] * z_p_prime_base_147[1]
   - 0.5 p[0] * z_p_prime_base_147[1] + 0.5 p[1] * z_p_prime_base_147[1]
   + 0.5 p[2] * z_p_prime_base_147[1] + 0.5 p[3] * z_p_prime_base_147[1]
   - 0.5 p[4] * z_p_prime_base_147[1] ] <= 0
 infeasibility_p_prime_base_147_row_2_leq_0: 5e-06 z_p_prime_base_147[2]
   + [ 0.5 alpha[0] * z_p_prime_base_147[2]
   - 0.5 alpha[1] * z_p_prime_base_147[2]
   + 0.5 alpha[2] * z_p_prime_base_147[2]
   + 0.5 alpha[3] * z_p_prime_base_147[2]
   + 0.5 alpha[4] * z_p_prime_base_147[2]
   + 0.5 p[0] * z_p_prime_base_147[2] + 0.5 p[1] * z_p_prime_base_147[2]
   - 0.5 p[2] * z_p_prime_base_147[2] - 0.5 p[3] * z_p_prime_base_147[2]
   + 0.5 p[4] * z_p_prime_base_147[2] ] <= 0
 infeasibility_p_prime_base_147_row_3_leq_0: 5e-06 z_p_prime_base_147[3]
   + [ - alpha[2] * z_p_prime_base_147[3] + p[2] * z_p_prime_base_147[3] ]
   <= 0
 infeasibility_p_prime_base_147_row_4_leq_0: 5e-06 z_p_prime_base_147[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_147[4]
   + 0.5 alpha[1] * z_p_prime_base_147[4]
   + 0.5 alpha[2] * z_p_prime_base_147[4]
   - 0.5 alpha[3] * z_p_prime_base_147[4]
   + 0.5 alpha[4] * z_p_prime_base_147[4]
   - 0.5 p[0] * z_p_prime_base_147[4] - 0.5 p[1] * z_p_prime_base_147[4]
   - 0.5 p[2] * z_p_prime_base_147[4] + 0.5 p[3] * z_p_prime_base_147[4]
   + 0.5 p[4] * z_p_prime_base_147[4] ] <= 0
 infeasibility_p_base_148_row_0_leq_0: 5e-06 z_p_base_148[0] + [
   p[0] * z_p_base_148[0] ] <= 0
 infeasibility_p_base_148_row_1_leq_0: 5e-06 z_p_base_148[1] + [
   p[1] * z_p_base_148[1] ] <= 0
 infeasibility_p_base_148_row_2_leq_0: 5e-06 z_p_base_148[2] + [
   - 0.5 p[0] * z_p_base_148[2] - 0.5 p[1] * z_p_base_148[2]
   + 0.5 p[2] * z_p_base_148[2] + 0.5 p[3] * z_p_base_148[2]
   - 0.5 p[4] * z_p_base_148[2] ] <= 0
 infeasibility_p_base_148_row_3_leq_0: 5e-06 z_p_base_148[3] + [
   0.5 p[0] * z_p_base_148[3] + 0.5 p[1] * z_p_base_148[3]
   + 0.5 p[2] * z_p_base_148[3] - 0.5 p[3] * z_p_base_148[3]
   + 0.5 p[4] * z_p_base_148[3] ] <= 0
 infeasibility_p_base_148_row_4_leq_0: 5e-06 z_p_base_148[4] + [
   - 0.5 p[0] * z_p_base_148[4] - 0.5 p[1] * z_p_base_148[4]
   - 0.5 p[2] * z_p_base_148[4] + 0.5 p[3] * z_p_base_148[4]
   + 0.5 p[4] * z_p_base_148[4] ] <= 0
 infeasibility_p_prime_base_148_row_0_leq_0: 5e-06 z_p_prime_base_148[0]
   + [ alpha[0] * z_p_prime_base_148[0] + p[0] * z_p_prime_base_148[0] ]
   <= 0
 infeasibility_p_prime_base_148_row_1_leq_0: 5e-06 z_p_prime_base_148[1]
   + [ - alpha[1] * z_p_prime_base_148[1] + p[1] * z_p_prime_base_148[1] ]
   <= 0
 infeasibility_p_prime_base_148_row_2_leq_0: 5e-06 z_p_prime_base_148[2]
   + [ - 0.5 alpha[0] * z_p_prime_base_148[2]
   + 0.5 alpha[1] * z_p_prime_base_148[2]
   - 0.5 alpha[2] * z_p_prime_base_148[2]
   - 0.5 alpha[3] * z_p_prime_base_148[2]
   - 0.5 alpha[4] * z_p_prime_base_148[2]
   - 0.5 p[0] * z_p_prime_base_148[2] - 0.5 p[1] * z_p_prime_base_148[2]
   + 0.5 p[2] * z_p_prime_base_148[2] + 0.5 p[3] * z_p_prime_base_148[2]
   - 0.5 p[4] * z_p_prime_base_148[2] ] <= 0
 infeasibility_p_prime_base_148_row_3_leq_0: 5e-06 z_p_prime_base_148[3]
   + [ 0.5 alpha[0] * z_p_prime_base_148[3]
   - 0.5 alpha[1] * z_p_prime_base_148[3]
   - 0.5 alpha[2] * z_p_prime_base_148[3]
   + 0.5 alpha[3] * z_p_prime_base_148[3]
   + 0.5 alpha[4] * z_p_prime_base_148[3]
   + 0.5 p[0] * z_p_prime_base_148[3] + 0.5 p[1] * z_p_prime_base_148[3]
   + 0.5 p[2] * z_p_prime_base_148[3] - 0.5 p[3] * z_p_prime_base_148[3]
   + 0.5 p[4] * z_p_prime_base_148[3] ] <= 0
 infeasibility_p_prime_base_148_row_4_leq_0: 5e-06 z_p_prime_base_148[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_148[4]
   + 0.5 alpha[1] * z_p_prime_base_148[4]
   + 0.5 alpha[2] * z_p_prime_base_148[4]
   - 0.5 alpha[3] * z_p_prime_base_148[4]
   + 0.5 alpha[4] * z_p_prime_base_148[4]
   - 0.5 p[0] * z_p_prime_base_148[4] - 0.5 p[1] * z_p_prime_base_148[4]
   - 0.5 p[2] * z_p_prime_base_148[4] + 0.5 p[3] * z_p_prime_base_148[4]
   + 0.5 p[4] * z_p_prime_base_148[4] ] <= 0
 infeasibility_p_base_149_row_0_leq_0: 5e-06 z_p_base_149[0] + [
   p[0] * z_p_base_149[0] ] <= 0
 infeasibility_p_base_149_row_1_leq_0: 5e-06 z_p_base_149[1] + [
   p[1] * z_p_base_149[1] ] <= 0
 infeasibility_p_base_149_row_2_leq_0: 5e-06 z_p_base_149[2] + [
   - 0.5 p[0] * z_p_base_149[2] + 0.5 p[1] * z_p_base_149[2]
   + 0.5 p[2] * z_p_base_149[2] + 0.5 p[3] * z_p_base_149[2]
   - 0.5 p[4] * z_p_base_149[2] ] <= 0
 infeasibility_p_base_149_row_3_leq_0: 5e-06 z_p_base_149[3] + [
   0.5 p[0] * z_p_base_149[3] - 0.5 p[1] * z_p_base_149[3]
   + 0.5 p[2] * z_p_base_149[3] - 0.5 p[3] * z_p_base_149[3]
   + 0.5 p[4] * z_p_base_149[3] ] <= 0
 infeasibility_p_base_149_row_4_leq_0: 5e-06 z_p_base_149[4] + [
   - 0.5 p[0] * z_p_base_149[4] - 0.5 p[1] * z_p_base_149[4]
   - 0.5 p[2] * z_p_base_149[4] + 0.5 p[3] * z_p_base_149[4]
   + 0.5 p[4] * z_p_base_149[4] ] <= 0
 infeasibility_p_prime_base_149_row_0_leq_0: 5e-06 z_p_prime_base_149[0]
   + [ alpha[0] * z_p_prime_base_149[0] + p[0] * z_p_prime_base_149[0] ]
   <= 0
 infeasibility_p_prime_base_149_row_1_leq_0: 5e-06 z_p_prime_base_149[1]
   + [ - alpha[1] * z_p_prime_base_149[1] + p[1] * z_p_prime_base_149[1] ]
   <= 0
 infeasibility_p_prime_base_149_row_2_leq_0: 5e-06 z_p_prime_base_149[2]
   + [ - 0.5 alpha[0] * z_p_prime_base_149[2]
   - 0.5 alpha[1] * z_p_prime_base_149[2]
   - 0.5 alpha[2] * z_p_prime_base_149[2]
   - 0.5 alpha[3] * z_p_prime_base_149[2]
   - 0.5 alpha[4] * z_p_prime_base_149[2]
   - 0.5 p[0] * z_p_prime_base_149[2] + 0.5 p[1] * z_p_prime_base_149[2]
   + 0.5 p[2] * z_p_prime_base_149[2] + 0.5 p[3] * z_p_prime_base_149[2]
   - 0.5 p[4] * z_p_prime_base_149[2] ] <= 0
 infeasibility_p_prime_base_149_row_3_leq_0: 5e-06 z_p_prime_base_149[3]
   + [ 0.5 alpha[0] * z_p_prime_base_149[3]
   + 0.5 alpha[1] * z_p_prime_base_149[3]
   - 0.5 alpha[2] * z_p_prime_base_149[3]
   + 0.5 alpha[3] * z_p_prime_base_149[3]
   + 0.5 alpha[4] * z_p_prime_base_149[3]
   + 0.5 p[0] * z_p_prime_base_149[3] - 0.5 p[1] * z_p_prime_base_149[3]
   + 0.5 p[2] * z_p_prime_base_149[3] - 0.5 p[3] * z_p_prime_base_149[3]
   + 0.5 p[4] * z_p_prime_base_149[3] ] <= 0
 infeasibility_p_prime_base_149_row_4_leq_0: 5e-06 z_p_prime_base_149[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_149[4]
   + 0.5 alpha[1] * z_p_prime_base_149[4]
   + 0.5 alpha[2] * z_p_prime_base_149[4]
   - 0.5 alpha[3] * z_p_prime_base_149[4]
   + 0.5 alpha[4] * z_p_prime_base_149[4]
   - 0.5 p[0] * z_p_prime_base_149[4] - 0.5 p[1] * z_p_prime_base_149[4]
   - 0.5 p[2] * z_p_prime_base_149[4] + 0.5 p[3] * z_p_prime_base_149[4]
   + 0.5 p[4] * z_p_prime_base_149[4] ] <= 0
 infeasibility_p_base_150_row_0_leq_0: 5e-06 z_p_base_150[0] + [
   p[0] * z_p_base_150[0] ] <= 0
 infeasibility_p_base_150_row_1_leq_0: 5e-06 z_p_base_150[1] + [
   0.5 p[0] * z_p_base_150[1] + 0.5 p[1] * z_p_base_150[1]
   + 0.5 p[2] * z_p_base_150[1] - 0.5 p[3] * z_p_base_150[1]
   - 0.5 p[4] * z_p_base_150[1] ] <= 0
 infeasibility_p_base_150_row_2_leq_0: 5e-06 z_p_base_150[2] + [
   0.5 p[0] * z_p_base_150[2] + 0.5 p[1] * z_p_base_150[2]
   - 0.5 p[2] * z_p_base_150[2] + 0.5 p[3] * z_p_base_150[2]
   - 0.5 p[4] * z_p_base_150[2] ] <= 0
 infeasibility_p_base_150_row_3_leq_0: 5e-06 z_p_base_150[3] + [
   - p[0] * z_p_base_150[3] + p[4] * z_p_base_150[3] ] <= 0
 infeasibility_p_base_150_row_4_leq_0: 5e-06 z_p_base_150[4] + [
   - 0.5 p[0] * z_p_base_150[4] - 0.5 p[1] * z_p_base_150[4]
   + 0.5 p[2] * z_p_base_150[4] + 0.5 p[3] * z_p_base_150[4]
   + 0.5 p[4] * z_p_base_150[4] ] <= 0
 infeasibility_p_prime_base_150_row_0_leq_0: 5e-06 z_p_prime_base_150[0]
   + [ alpha[0] * z_p_prime_base_150[0] + p[0] * z_p_prime_base_150[0] ]
   <= 0
 infeasibility_p_prime_base_150_row_1_leq_0: 5e-06 z_p_prime_base_150[1]
   + [ 0.5 alpha[0] * z_p_prime_base_150[1]
   - 0.5 alpha[1] * z_p_prime_base_150[1]
   - 0.5 alpha[2] * z_p_prime_base_150[1]
   + 0.5 alpha[3] * z_p_prime_base_150[1]
   - 0.5 alpha[4] * z_p_prime_base_150[1]
   + 0.5 p[0] * z_p_prime_base_150[1] + 0.5 p[1] * z_p_prime_base_150[1]
   + 0.5 p[2] * z_p_prime_base_150[1] - 0.5 p[3] * z_p_prime_base_150[1]
   - 0.5 p[4] * z_p_prime_base_150[1] ] <= 0
 infeasibility_p_prime_base_150_row_2_leq_0: 5e-06 z_p_prime_base_150[2]
   + [ 0.5 alpha[0] * z_p_prime_base_150[2]
   - 0.5 alpha[1] * z_p_prime_base_150[2]
   + 0.5 alpha[2] * z_p_prime_base_150[2]
   - 0.5 alpha[3] * z_p_prime_base_150[2]
   - 0.5 alpha[4] * z_p_prime_base_150[2]
   + 0.5 p[0] * z_p_prime_base_150[2] + 0.5 p[1] * z_p_prime_base_150[2]
   - 0.5 p[2] * z_p_prime_base_150[2] + 0.5 p[3] * z_p_prime_base_150[2]
   - 0.5 p[4] * z_p_prime_base_150[2] ] <= 0
 infeasibility_p_prime_base_150_row_3_leq_0: 5e-06 z_p_prime_base_150[3]
   + [ - alpha[0] * z_p_prime_base_150[3]
   + alpha[4] * z_p_prime_base_150[3] - p[0] * z_p_prime_base_150[3]
   + p[4] * z_p_prime_base_150[3] ] <= 0
 infeasibility_p_prime_base_150_row_4_leq_0: 5e-06 z_p_prime_base_150[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_150[4]
   + 0.5 alpha[1] * z_p_prime_base_150[4]
   - 0.5 alpha[2] * z_p_prime_base_150[4]
   - 0.5 alpha[3] * z_p_prime_base_150[4]
   + 0.5 alpha[4] * z_p_prime_base_150[4]
   - 0.5 p[0] * z_p_prime_base_150[4] - 0.5 p[1] * z_p_prime_base_150[4]
   + 0.5 p[2] * z_p_prime_base_150[4] + 0.5 p[3] * z_p_prime_base_150[4]
   + 0.5 p[4] * z_p_prime_base_150[4] ] <= 0
 infeasibility_p_base_151_row_0_leq_0: 5e-06 z_p_base_151[0] + [
   p[0] * z_p_base_151[0] ] <= 0
 infeasibility_p_base_151_row_1_leq_0: 5e-06 z_p_base_151[1] + [
   0.5 p[0] * z_p_base_151[1] + 0.5 p[1] * z_p_base_151[1]
   + 0.5 p[2] * z_p_base_151[1] - 0.5 p[3] * z_p_base_151[1]
   - 0.5 p[4] * z_p_base_151[1] ] <= 0
 infeasibility_p_base_151_row_2_leq_0: 5e-06 z_p_base_151[2] + [
   p[3] * z_p_base_151[2] ] <= 0
 infeasibility_p_base_151_row_3_leq_0: 5e-06 z_p_base_151[3] + [
   - 0.5 p[0] * z_p_base_151[3] + 0.5 p[1] * z_p_base_151[3]
   - 0.5 p[2] * z_p_base_151[3] - 0.5 p[3] * z_p_base_151[3]
   + 0.5 p[4] * z_p_base_151[3] ] <= 0
 infeasibility_p_base_151_row_4_leq_0: 5e-06 z_p_base_151[4] + [
   - 0.5 p[0] * z_p_base_151[4] - 0.5 p[1] * z_p_base_151[4]
   + 0.5 p[2] * z_p_base_151[4] + 0.5 p[3] * z_p_base_151[4]
   + 0.5 p[4] * z_p_base_151[4] ] <= 0
 infeasibility_p_prime_base_151_row_0_leq_0: 5e-06 z_p_prime_base_151[0]
   + [ alpha[0] * z_p_prime_base_151[0] + p[0] * z_p_prime_base_151[0] ]
   <= 0
 infeasibility_p_prime_base_151_row_1_leq_0: 5e-06 z_p_prime_base_151[1]
   + [ 0.5 alpha[0] * z_p_prime_base_151[1]
   - 0.5 alpha[1] * z_p_prime_base_151[1]
   - 0.5 alpha[2] * z_p_prime_base_151[1]
   + 0.5 alpha[3] * z_p_prime_base_151[1]
   - 0.5 alpha[4] * z_p_prime_base_151[1]
   + 0.5 p[0] * z_p_prime_base_151[1] + 0.5 p[1] * z_p_prime_base_151[1]
   + 0.5 p[2] * z_p_prime_base_151[1] - 0.5 p[3] * z_p_prime_base_151[1]
   - 0.5 p[4] * z_p_prime_base_151[1] ] <= 0
 infeasibility_p_prime_base_151_row_2_leq_0: 5e-06 z_p_prime_base_151[2]
   + [ - alpha[3] * z_p_prime_base_151[2] + p[3] * z_p_prime_base_151[2] ]
   <= 0
 infeasibility_p_prime_base_151_row_3_leq_0: 5e-06 z_p_prime_base_151[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_151[3]
   - 0.5 alpha[1] * z_p_prime_base_151[3]
   + 0.5 alpha[2] * z_p_prime_base_151[3]
   + 0.5 alpha[3] * z_p_prime_base_151[3]
   + 0.5 alpha[4] * z_p_prime_base_151[3]
   - 0.5 p[0] * z_p_prime_base_151[3] + 0.5 p[1] * z_p_prime_base_151[3]
   - 0.5 p[2] * z_p_prime_base_151[3] - 0.5 p[3] * z_p_prime_base_151[3]
   + 0.5 p[4] * z_p_prime_base_151[3] ] <= 0
 infeasibility_p_prime_base_151_row_4_leq_0: 5e-06 z_p_prime_base_151[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_151[4]
   + 0.5 alpha[1] * z_p_prime_base_151[4]
   - 0.5 alpha[2] * z_p_prime_base_151[4]
   - 0.5 alpha[3] * z_p_prime_base_151[4]
   + 0.5 alpha[4] * z_p_prime_base_151[4]
   - 0.5 p[0] * z_p_prime_base_151[4] - 0.5 p[1] * z_p_prime_base_151[4]
   + 0.5 p[2] * z_p_prime_base_151[4] + 0.5 p[3] * z_p_prime_base_151[4]
   + 0.5 p[4] * z_p_prime_base_151[4] ] <= 0
 infeasibility_p_base_152_row_0_leq_0: 5e-06 z_p_base_152[0] + [
   p[0] * z_p_base_152[0] ] <= 0
 infeasibility_p_base_152_row_1_leq_0: 5e-06 z_p_base_152[1] + [
   p[2] * z_p_base_152[1] ] <= 0
 infeasibility_p_base_152_row_2_leq_0: 5e-06 z_p_base_152[2] + [
   0.5 p[0] * z_p_base_152[2] + 0.5 p[1] * z_p_base_152[2]
   - 0.5 p[2] * z_p_base_152[2] + 0.5 p[3] * z_p_base_152[2]
   - 0.5 p[4] * z_p_base_152[2] ] <= 0
 infeasibility_p_base_152_row_3_leq_0: 5e-06 z_p_base_152[3] + [
   - 0.5 p[0] * z_p_base_152[3] + 0.5 p[1] * z_p_base_152[3]
   - 0.5 p[2] * z_p_base_152[3] - 0.5 p[3] * z_p_base_152[3]
   + 0.5 p[4] * z_p_base_152[3] ] <= 0
 infeasibility_p_base_152_row_4_leq_0: 5e-06 z_p_base_152[4] + [
   - 0.5 p[0] * z_p_base_152[4] - 0.5 p[1] * z_p_base_152[4]
   + 0.5 p[2] * z_p_base_152[4] + 0.5 p[3] * z_p_base_152[4]
   + 0.5 p[4] * z_p_base_152[4] ] <= 0
 infeasibility_p_prime_base_152_row_0_leq_0: 5e-06 z_p_prime_base_152[0]
   + [ alpha[0] * z_p_prime_base_152[0] + p[0] * z_p_prime_base_152[0] ]
   <= 0
 infeasibility_p_prime_base_152_row_1_leq_0: 5e-06 z_p_prime_base_152[1]
   + [ - alpha[2] * z_p_prime_base_152[1] + p[2] * z_p_prime_base_152[1] ]
   <= 0
 infeasibility_p_prime_base_152_row_2_leq_0: 5e-06 z_p_prime_base_152[2]
   + [ 0.5 alpha[0] * z_p_prime_base_152[2]
   - 0.5 alpha[1] * z_p_prime_base_152[2]
   + 0.5 alpha[2] * z_p_prime_base_152[2]
   - 0.5 alpha[3] * z_p_prime_base_152[2]
   - 0.5 alpha[4] * z_p_prime_base_152[2]
   + 0.5 p[0] * z_p_prime_base_152[2] + 0.5 p[1] * z_p_prime_base_152[2]
   - 0.5 p[2] * z_p_prime_base_152[2] + 0.5 p[3] * z_p_prime_base_152[2]
   - 0.5 p[4] * z_p_prime_base_152[2] ] <= 0
 infeasibility_p_prime_base_152_row_3_leq_0: 5e-06 z_p_prime_base_152[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_152[3]
   - 0.5 alpha[1] * z_p_prime_base_152[3]
   + 0.5 alpha[2] * z_p_prime_base_152[3]
   + 0.5 alpha[3] * z_p_prime_base_152[3]
   + 0.5 alpha[4] * z_p_prime_base_152[3]
   - 0.5 p[0] * z_p_prime_base_152[3] + 0.5 p[1] * z_p_prime_base_152[3]
   - 0.5 p[2] * z_p_prime_base_152[3] - 0.5 p[3] * z_p_prime_base_152[3]
   + 0.5 p[4] * z_p_prime_base_152[3] ] <= 0
 infeasibility_p_prime_base_152_row_4_leq_0: 5e-06 z_p_prime_base_152[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_152[4]
   + 0.5 alpha[1] * z_p_prime_base_152[4]
   - 0.5 alpha[2] * z_p_prime_base_152[4]
   - 0.5 alpha[3] * z_p_prime_base_152[4]
   + 0.5 alpha[4] * z_p_prime_base_152[4]
   - 0.5 p[0] * z_p_prime_base_152[4] - 0.5 p[1] * z_p_prime_base_152[4]
   + 0.5 p[2] * z_p_prime_base_152[4] + 0.5 p[3] * z_p_prime_base_152[4]
   + 0.5 p[4] * z_p_prime_base_152[4] ] <= 0
 infeasibility_p_base_153_row_0_leq_0: 5e-06 z_p_base_153[0] + [
   p[0] * z_p_base_153[0] ] <= 0
 infeasibility_p_base_153_row_1_leq_0: 5e-06 z_p_base_153[1] + [
   0.5 p[0] * z_p_base_153[1] + 0.5 p[1] * z_p_base_153[1]
   + 0.5 p[2] * z_p_base_153[1] - 0.5 p[3] * z_p_base_153[1]
   - 0.5 p[4] * z_p_base_153[1] ] <= 0
 infeasibility_p_base_153_row_2_leq_0: 5e-06 z_p_base_153[2] + [
   - 0.5 p[0] * z_p_base_153[2] + 0.5 p[1] * z_p_base_153[2]
   - 0.5 p[2] * z_p_base_153[2] + 0.5 p[3] * z_p_base_153[2]
   + 0.5 p[4] * z_p_base_153[2] ] <= 0
 infeasibility_p_base_153_row_3_leq_0: 5e-06 z_p_base_153[3] + [
   0.5 p[0] * z_p_base_153[3] - 0.5 p[1] * z_p_base_153[3]
   + 0.5 p[2] * z_p_base_153[3] + 0.5 p[3] * z_p_base_153[3]
   - 0.5 p[4] * z_p_base_153[3] ] <= 0
 infeasibility_p_base_153_row_4_leq_0: 5e-06 z_p_base_153[4] + [
   - p[0] * z_p_base_153[4] + p[4] * z_p_base_153[4] ] <= 0
 infeasibility_p_prime_base_153_row_0_leq_0: 5e-06 z_p_prime_base_153[0]
   + [ alpha[0] * z_p_prime_base_153[0] + p[0] * z_p_prime_base_153[0] ]
   <= 0
 infeasibility_p_prime_base_153_row_1_leq_0: 5e-06 z_p_prime_base_153[1]
   + [ 0.5 alpha[0] * z_p_prime_base_153[1]
   - 0.5 alpha[1] * z_p_prime_base_153[1]
   - 0.5 alpha[2] * z_p_prime_base_153[1]
   + 0.5 alpha[3] * z_p_prime_base_153[1]
   - 0.5 alpha[4] * z_p_prime_base_153[1]
   + 0.5 p[0] * z_p_prime_base_153[1] + 0.5 p[1] * z_p_prime_base_153[1]
   + 0.5 p[2] * z_p_prime_base_153[1] - 0.5 p[3] * z_p_prime_base_153[1]
   - 0.5 p[4] * z_p_prime_base_153[1] ] <= 0
 infeasibility_p_prime_base_153_row_2_leq_0: 5e-06 z_p_prime_base_153[2]
   + [ - 0.5 alpha[0] * z_p_prime_base_153[2]
   - 0.5 alpha[1] * z_p_prime_base_153[2]
   + 0.5 alpha[2] * z_p_prime_base_153[2]
   - 0.5 alpha[3] * z_p_prime_base_153[2]
   + 0.5 alpha[4] * z_p_prime_base_153[2]
   - 0.5 p[0] * z_p_prime_base_153[2] + 0.5 p[1] * z_p_prime_base_153[2]
   - 0.5 p[2] * z_p_prime_base_153[2] + 0.5 p[3] * z_p_prime_base_153[2]
   + 0.5 p[4] * z_p_prime_base_153[2] ] <= 0
 infeasibility_p_prime_base_153_row_3_leq_0: 5e-06 z_p_prime_base_153[3]
   + [ 0.5 alpha[0] * z_p_prime_base_153[3]
   + 0.5 alpha[1] * z_p_prime_base_153[3]
   - 0.5 alpha[2] * z_p_prime_base_153[3]
   - 0.5 alpha[3] * z_p_prime_base_153[3]
   - 0.5 alpha[4] * z_p_prime_base_153[3]
   + 0.5 p[0] * z_p_prime_base_153[3] - 0.5 p[1] * z_p_prime_base_153[3]
   + 0.5 p[2] * z_p_prime_base_153[3] + 0.5 p[3] * z_p_prime_base_153[3]
   - 0.5 p[4] * z_p_prime_base_153[3] ] <= 0
 infeasibility_p_prime_base_153_row_4_leq_0: 5e-06 z_p_prime_base_153[4]
   + [ - alpha[0] * z_p_prime_base_153[4]
   + alpha[4] * z_p_prime_base_153[4] - p[0] * z_p_prime_base_153[4]
   + p[4] * z_p_prime_base_153[4] ] <= 0
 infeasibility_p_base_154_row_0_leq_0: 5e-06 z_p_base_154[0] + [
   p[0] * z_p_base_154[0] ] <= 0
 infeasibility_p_base_154_row_1_leq_0: 5e-06 z_p_base_154[1] + [
   - 0.5 p[0] * z_p_base_154[1] + 0.5 p[1] * z_p_base_154[1]
   + 0.5 p[2] * z_p_base_154[1] - 0.5 p[3] * z_p_base_154[1]
   + 0.5 p[4] * z_p_base_154[1] ] <= 0
 infeasibility_p_base_154_row_2_leq_0: 5e-06 z_p_base_154[2] + [
   0.5 p[0] * z_p_base_154[2] + 0.5 p[1] * z_p_base_154[2]
   - 0.5 p[2] * z_p_base_154[2] + 0.5 p[3] * z_p_base_154[2]
   - 0.5 p[4] * z_p_base_154[2] ] <= 0
 infeasibility_p_base_154_row_3_leq_0: 5e-06 z_p_base_154[3] + [
   0.5 p[0] * z_p_base_154[3] - 0.5 p[1] * z_p_base_154[3]
   + 0.5 p[2] * z_p_base_154[3] + 0.5 p[3] * z_p_base_154[3]
   - 0.5 p[4] * z_p_base_154[3] ] <= 0
 infeasibility_p_base_154_row_4_leq_0: 5e-06 z_p_base_154[4] + [
   - p[0] * z_p_base_154[4] + p[4] * z_p_base_154[4] ] <= 0
 infeasibility_p_prime_base_154_row_0_leq_0: 5e-06 z_p_prime_base_154[0]
   + [ alpha[0] * z_p_prime_base_154[0] + p[0] * z_p_prime_base_154[0] ]
   <= 0
 infeasibility_p_prime_base_154_row_1_leq_0: 5e-06 z_p_prime_base_154[1]
   + [ - 0.5 alpha[0] * z_p_prime_base_154[1]
   - 0.5 alpha[1] * z_p_prime_base_154[1]
   - 0.5 alpha[2] * z_p_prime_base_154[1]
   + 0.5 alpha[3] * z_p_prime_base_154[1]
   + 0.5 alpha[4] * z_p_prime_base_154[1]
   - 0.5 p[0] * z_p_prime_base_154[1] + 0.5 p[1] * z_p_prime_base_154[1]
   + 0.5 p[2] * z_p_prime_base_154[1] - 0.5 p[3] * z_p_prime_base_154[1]
   + 0.5 p[4] * z_p_prime_base_154[1] ] <= 0
 infeasibility_p_prime_base_154_row_2_leq_0: 5e-06 z_p_prime_base_154[2]
   + [ 0.5 alpha[0] * z_p_prime_base_154[2]
   - 0.5 alpha[1] * z_p_prime_base_154[2]
   + 0.5 alpha[2] * z_p_prime_base_154[2]
   - 0.5 alpha[3] * z_p_prime_base_154[2]
   - 0.5 alpha[4] * z_p_prime_base_154[2]
   + 0.5 p[0] * z_p_prime_base_154[2] + 0.5 p[1] * z_p_prime_base_154[2]
   - 0.5 p[2] * z_p_prime_base_154[2] + 0.5 p[3] * z_p_prime_base_154[2]
   - 0.5 p[4] * z_p_prime_base_154[2] ] <= 0
 infeasibility_p_prime_base_154_row_3_leq_0: 5e-06 z_p_prime_base_154[3]
   + [ 0.5 alpha[0] * z_p_prime_base_154[3]
   + 0.5 alpha[1] * z_p_prime_base_154[3]
   - 0.5 alpha[2] * z_p_prime_base_154[3]
   - 0.5 alpha[3] * z_p_prime_base_154[3]
   - 0.5 alpha[4] * z_p_prime_base_154[3]
   + 0.5 p[0] * z_p_prime_base_154[3] - 0.5 p[1] * z_p_prime_base_154[3]
   + 0.5 p[2] * z_p_prime_base_154[3] + 0.5 p[3] * z_p_prime_base_154[3]
   - 0.5 p[4] * z_p_prime_base_154[3] ] <= 0
 infeasibility_p_prime_base_154_row_4_leq_0: 5e-06 z_p_prime_base_154[4]
   + [ - alpha[0] * z_p_prime_base_154[4]
   + alpha[4] * z_p_prime_base_154[4] - p[0] * z_p_prime_base_154[4]
   + p[4] * z_p_prime_base_154[4] ] <= 0
 infeasibility_p_base_155_row_0_leq_0: 5e-06 z_p_base_155[0] + [
   p[0] * z_p_base_155[0] ] <= 0
 infeasibility_p_base_155_row_1_leq_0: 5e-06 z_p_base_155[1] + [
   0.5 p[0] * z_p_base_155[1] + 0.5 p[1] * z_p_base_155[1]
   + 0.5 p[2] * z_p_base_155[1] - 0.5 p[3] * z_p_base_155[1]
   - 0.5 p[4] * z_p_base_155[1] ] <= 0
 infeasibility_p_base_155_row_2_leq_0: 5e-06 z_p_base_155[2] + [
   - 0.5 p[0] * z_p_base_155[2] + 0.5 p[1] * z_p_base_155[2]
   - 0.5 p[2] * z_p_base_155[2] + 0.5 p[3] * z_p_base_155[2]
   + 0.5 p[4] * z_p_base_155[2] ] <= 0
 infeasibility_p_base_155_row_3_leq_0: 5e-06 z_p_base_155[3] + [
   p[3] * z_p_base_155[3] ] <= 0
 infeasibility_p_base_155_row_4_leq_0: 5e-06 z_p_base_155[4] + [
   - 0.5 p[0] * z_p_base_155[4] - 0.5 p[1] * z_p_base_155[4]
   + 0.5 p[2] * z_p_base_155[4] - 0.5 p[3] * z_p_base_155[4]
   + 0.5 p[4] * z_p_base_155[4] ] <= 0
 infeasibility_p_prime_base_155_row_0_leq_0: 5e-06 z_p_prime_base_155[0]
   + [ alpha[0] * z_p_prime_base_155[0] + p[0] * z_p_prime_base_155[0] ]
   <= 0
 infeasibility_p_prime_base_155_row_1_leq_0: 5e-06 z_p_prime_base_155[1]
   + [ 0.5 alpha[0] * z_p_prime_base_155[1]
   - 0.5 alpha[1] * z_p_prime_base_155[1]
   - 0.5 alpha[2] * z_p_prime_base_155[1]
   + 0.5 alpha[3] * z_p_prime_base_155[1]
   - 0.5 alpha[4] * z_p_prime_base_155[1]
   + 0.5 p[0] * z_p_prime_base_155[1] + 0.5 p[1] * z_p_prime_base_155[1]
   + 0.5 p[2] * z_p_prime_base_155[1] - 0.5 p[3] * z_p_prime_base_155[1]
   - 0.5 p[4] * z_p_prime_base_155[1] ] <= 0
 infeasibility_p_prime_base_155_row_2_leq_0: 5e-06 z_p_prime_base_155[2]
   + [ - 0.5 alpha[0] * z_p_prime_base_155[2]
   - 0.5 alpha[1] * z_p_prime_base_155[2]
   + 0.5 alpha[2] * z_p_prime_base_155[2]
   - 0.5 alpha[3] * z_p_prime_base_155[2]
   + 0.5 alpha[4] * z_p_prime_base_155[2]
   - 0.5 p[0] * z_p_prime_base_155[2] + 0.5 p[1] * z_p_prime_base_155[2]
   - 0.5 p[2] * z_p_prime_base_155[2] + 0.5 p[3] * z_p_prime_base_155[2]
   + 0.5 p[4] * z_p_prime_base_155[2] ] <= 0
 infeasibility_p_prime_base_155_row_3_leq_0: 5e-06 z_p_prime_base_155[3]
   + [ - alpha[3] * z_p_prime_base_155[3] + p[3] * z_p_prime_base_155[3] ]
   <= 0
 infeasibility_p_prime_base_155_row_4_leq_0: 5e-06 z_p_prime_base_155[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_155[4]
   + 0.5 alpha[1] * z_p_prime_base_155[4]
   - 0.5 alpha[2] * z_p_prime_base_155[4]
   + 0.5 alpha[3] * z_p_prime_base_155[4]
   + 0.5 alpha[4] * z_p_prime_base_155[4]
   - 0.5 p[0] * z_p_prime_base_155[4] - 0.5 p[1] * z_p_prime_base_155[4]
   + 0.5 p[2] * z_p_prime_base_155[4] - 0.5 p[3] * z_p_prime_base_155[4]
   + 0.5 p[4] * z_p_prime_base_155[4] ] <= 0
 infeasibility_p_base_156_row_0_leq_0: 5e-06 z_p_base_156[0] + [
   p[0] * z_p_base_156[0] ] <= 0
 infeasibility_p_base_156_row_1_leq_0: 5e-06 z_p_base_156[1] + [
   0.5 p[0] * z_p_base_156[1] + 0.5 p[1] * z_p_base_156[1]
   + 0.5 p[2] * z_p_base_156[1] + 0.5 p[3] * z_p_base_156[1]
   - 0.5 p[4] * z_p_base_156[1] ] <= 0
 infeasibility_p_base_156_row_2_leq_0: 5e-06 z_p_base_156[2] + [
   - 0.5 p[0] * z_p_base_156[2] + 0.5 p[1] * z_p_base_156[2]
   - 0.5 p[2] * z_p_base_156[2] - 0.5 p[3] * z_p_base_156[2]
   + 0.5 p[4] * z_p_base_156[2] ] <= 0
 infeasibility_p_base_156_row_3_leq_0: 5e-06 z_p_base_156[3] + [
   - 0.5 p[0] * z_p_base_156[3] - 0.5 p[1] * z_p_base_156[3]
   + 0.5 p[2] * z_p_base_156[3] - 0.5 p[3] * z_p_base_156[3]
   + 0.5 p[4] * z_p_base_156[3] ] <= 0
 infeasibility_p_base_156_row_4_leq_0: 5e-06 z_p_base_156[4] + [
   p[3] * z_p_base_156[4] ] <= 0
 infeasibility_p_prime_base_156_row_0_leq_0: 5e-06 z_p_prime_base_156[0]
   + [ alpha[0] * z_p_prime_base_156[0] + p[0] * z_p_prime_base_156[0] ]
   <= 0
 infeasibility_p_prime_base_156_row_1_leq_0: 5e-06 z_p_prime_base_156[1]
   + [ 0.5 alpha[0] * z_p_prime_base_156[1]
   - 0.5 alpha[1] * z_p_prime_base_156[1]
   - 0.5 alpha[2] * z_p_prime_base_156[1]
   - 0.5 alpha[3] * z_p_prime_base_156[1]
   - 0.5 alpha[4] * z_p_prime_base_156[1]
   + 0.5 p[0] * z_p_prime_base_156[1] + 0.5 p[1] * z_p_prime_base_156[1]
   + 0.5 p[2] * z_p_prime_base_156[1] + 0.5 p[3] * z_p_prime_base_156[1]
   - 0.5 p[4] * z_p_prime_base_156[1] ] <= 0
 infeasibility_p_prime_base_156_row_2_leq_0: 5e-06 z_p_prime_base_156[2]
   + [ - 0.5 alpha[0] * z_p_prime_base_156[2]
   - 0.5 alpha[1] * z_p_prime_base_156[2]
   + 0.5 alpha[2] * z_p_prime_base_156[2]
   + 0.5 alpha[3] * z_p_prime_base_156[2]
   + 0.5 alpha[4] * z_p_prime_base_156[2]
   - 0.5 p[0] * z_p_prime_base_156[2] + 0.5 p[1] * z_p_prime_base_156[2]
   - 0.5 p[2] * z_p_prime_base_156[2] - 0.5 p[3] * z_p_prime_base_156[2]
   + 0.5 p[4] * z_p_prime_base_156[2] ] <= 0
 infeasibility_p_prime_base_156_row_3_leq_0: 5e-06 z_p_prime_base_156[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_156[3]
   + 0.5 alpha[1] * z_p_prime_base_156[3]
   - 0.5 alpha[2] * z_p_prime_base_156[3]
   + 0.5 alpha[3] * z_p_prime_base_156[3]
   + 0.5 alpha[4] * z_p_prime_base_156[3]
   - 0.5 p[0] * z_p_prime_base_156[3] - 0.5 p[1] * z_p_prime_base_156[3]
   + 0.5 p[2] * z_p_prime_base_156[3] - 0.5 p[3] * z_p_prime_base_156[3]
   + 0.5 p[4] * z_p_prime_base_156[3] ] <= 0
 infeasibility_p_prime_base_156_row_4_leq_0: 5e-06 z_p_prime_base_156[4]
   + [ - alpha[3] * z_p_prime_base_156[4] + p[3] * z_p_prime_base_156[4] ]
   <= 0
 infeasibility_p_base_157_row_0_leq_0: 5e-06 z_p_base_157[0] + [
   p[0] * z_p_base_157[0] ] <= 0
 infeasibility_p_base_157_row_1_leq_0: 5e-06 z_p_base_157[1] + [
   p[1] * z_p_base_157[1] ] <= 0
 infeasibility_p_base_157_row_2_leq_0: 5e-06 z_p_base_157[2] + [
   0.5 p[0] * z_p_base_157[2] - 0.5 p[1] * z_p_base_157[2]
   + 0.5 p[2] * z_p_base_157[2] + 0.5 p[3] * z_p_base_157[2]
   - 0.5 p[4] * z_p_base_157[2] ] <= 0
 infeasibility_p_base_157_row_3_leq_0: 5e-06 z_p_base_157[3] + [
   - 0.5 p[0] * z_p_base_157[3] - 0.5 p[1] * z_p_base_157[3]
   + 0.5 p[2] * z_p_base_157[3] - 0.5 p[3] * z_p_base_157[3]
   + 0.5 p[4] * z_p_base_157[3] ] <= 0
 infeasibility_p_base_157_row_4_leq_0: 5e-06 z_p_base_157[4] + [
   - 0.5 p[0] * z_p_base_157[4] + 0.5 p[1] * z_p_base_157[4]
   - 0.5 p[2] * z_p_base_157[4] + 0.5 p[3] * z_p_base_157[4]
   + 0.5 p[4] * z_p_base_157[4] ] <= 0
 infeasibility_p_prime_base_157_row_0_leq_0: 5e-06 z_p_prime_base_157[0]
   + [ alpha[0] * z_p_prime_base_157[0] + p[0] * z_p_prime_base_157[0] ]
   <= 0
 infeasibility_p_prime_base_157_row_1_leq_0: 5e-06 z_p_prime_base_157[1]
   + [ - alpha[1] * z_p_prime_base_157[1] + p[1] * z_p_prime_base_157[1] ]
   <= 0
 infeasibility_p_prime_base_157_row_2_leq_0: 5e-06 z_p_prime_base_157[2]
   + [ 0.5 alpha[0] * z_p_prime_base_157[2]
   + 0.5 alpha[1] * z_p_prime_base_157[2]
   - 0.5 alpha[2] * z_p_prime_base_157[2]
   - 0.5 alpha[3] * z_p_prime_base_157[2]
   - 0.5 alpha[4] * z_p_prime_base_157[2]
   + 0.5 p[0] * z_p_prime_base_157[2] - 0.5 p[1] * z_p_prime_base_157[2]
   + 0.5 p[2] * z_p_prime_base_157[2] + 0.5 p[3] * z_p_prime_base_157[2]
   - 0.5 p[4] * z_p_prime_base_157[2] ] <= 0
 infeasibility_p_prime_base_157_row_3_leq_0: 5e-06 z_p_prime_base_157[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_157[3]
   + 0.5 alpha[1] * z_p_prime_base_157[3]
   - 0.5 alpha[2] * z_p_prime_base_157[3]
   + 0.5 alpha[3] * z_p_prime_base_157[3]
   + 0.5 alpha[4] * z_p_prime_base_157[3]
   - 0.5 p[0] * z_p_prime_base_157[3] - 0.5 p[1] * z_p_prime_base_157[3]
   + 0.5 p[2] * z_p_prime_base_157[3] - 0.5 p[3] * z_p_prime_base_157[3]
   + 0.5 p[4] * z_p_prime_base_157[3] ] <= 0
 infeasibility_p_prime_base_157_row_4_leq_0: 5e-06 z_p_prime_base_157[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_157[4]
   - 0.5 alpha[1] * z_p_prime_base_157[4]
   + 0.5 alpha[2] * z_p_prime_base_157[4]
   - 0.5 alpha[3] * z_p_prime_base_157[4]
   + 0.5 alpha[4] * z_p_prime_base_157[4]
   - 0.5 p[0] * z_p_prime_base_157[4] + 0.5 p[1] * z_p_prime_base_157[4]
   - 0.5 p[2] * z_p_prime_base_157[4] + 0.5 p[3] * z_p_prime_base_157[4]
   + 0.5 p[4] * z_p_prime_base_157[4] ] <= 0
 infeasibility_p_base_158_row_0_leq_0: 5e-06 z_p_base_158[0] + [
   p[0] * z_p_base_158[0] ] <= 0
 infeasibility_p_base_158_row_1_leq_0: 5e-06 z_p_base_158[1] + [
   0.5 p[0] * z_p_base_158[1] + 0.5 p[1] * z_p_base_158[1]
   - 0.5 p[2] * z_p_base_158[1] + 0.5 p[3] * z_p_base_158[1]
   - 0.5 p[4] * z_p_base_158[1] ] <= 0
 infeasibility_p_base_158_row_2_leq_0: 5e-06 z_p_base_158[2] + [
   - 0.5 p[0] * z_p_base_158[2] + 0.5 p[1] * z_p_base_158[2]
   + 0.5 p[2] * z_p_base_158[2] - 0.5 p[3] * z_p_base_158[2]
   + 0.5 p[4] * z_p_base_158[2] ] <= 0
 infeasibility_p_base_158_row_3_leq_0: 5e-06 z_p_base_158[3] + [
   p[2] * z_p_base_158[3] ] <= 0
 infeasibility_p_base_158_row_4_leq_0: 5e-06 z_p_base_158[4] + [
   - 0.5 p[0] * z_p_base_158[4] - 0.5 p[1] * z_p_base_158[4]
   - 0.5 p[2] * z_p_base_158[4] + 0.5 p[3] * z_p_base_158[4]
   + 0.5 p[4] * z_p_base_158[4] ] <= 0
 infeasibility_p_prime_base_158_row_0_leq_0: 5e-06 z_p_prime_base_158[0]
   + [ alpha[0] * z_p_prime_base_158[0] + p[0] * z_p_prime_base_158[0] ]
   <= 0
 infeasibility_p_prime_base_158_row_1_leq_0: 5e-06 z_p_prime_base_158[1]
   + [ 0.5 alpha[0] * z_p_prime_base_158[1]
   - 0.5 alpha[1] * z_p_prime_base_158[1]
   + 0.5 alpha[2] * z_p_prime_base_158[1]
   - 0.5 alpha[3] * z_p_prime_base_158[1]
   - 0.5 alpha[4] * z_p_prime_base_158[1]
   + 0.5 p[0] * z_p_prime_base_158[1] + 0.5 p[1] * z_p_prime_base_158[1]
   - 0.5 p[2] * z_p_prime_base_158[1] + 0.5 p[3] * z_p_prime_base_158[1]
   - 0.5 p[4] * z_p_prime_base_158[1] ] <= 0
 infeasibility_p_prime_base_158_row_2_leq_0: 5e-06 z_p_prime_base_158[2]
   + [ - 0.5 alpha[0] * z_p_prime_base_158[2]
   - 0.5 alpha[1] * z_p_prime_base_158[2]
   - 0.5 alpha[2] * z_p_prime_base_158[2]
   + 0.5 alpha[3] * z_p_prime_base_158[2]
   + 0.5 alpha[4] * z_p_prime_base_158[2]
   - 0.5 p[0] * z_p_prime_base_158[2] + 0.5 p[1] * z_p_prime_base_158[2]
   + 0.5 p[2] * z_p_prime_base_158[2] - 0.5 p[3] * z_p_prime_base_158[2]
   + 0.5 p[4] * z_p_prime_base_158[2] ] <= 0
 infeasibility_p_prime_base_158_row_3_leq_0: 5e-06 z_p_prime_base_158[3]
   + [ - alpha[2] * z_p_prime_base_158[3] + p[2] * z_p_prime_base_158[3] ]
   <= 0
 infeasibility_p_prime_base_158_row_4_leq_0: 5e-06 z_p_prime_base_158[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_158[4]
   + 0.5 alpha[1] * z_p_prime_base_158[4]
   + 0.5 alpha[2] * z_p_prime_base_158[4]
   - 0.5 alpha[3] * z_p_prime_base_158[4]
   + 0.5 alpha[4] * z_p_prime_base_158[4]
   - 0.5 p[0] * z_p_prime_base_158[4] - 0.5 p[1] * z_p_prime_base_158[4]
   - 0.5 p[2] * z_p_prime_base_158[4] + 0.5 p[3] * z_p_prime_base_158[4]
   + 0.5 p[4] * z_p_prime_base_158[4] ] <= 0
 infeasibility_p_base_159_row_0_leq_0: 5e-06 z_p_base_159[0] + [
   p[0] * z_p_base_159[0] ] <= 0
 infeasibility_p_base_159_row_1_leq_0: 5e-06 z_p_base_159[1] + [
   0.5 p[0] * z_p_base_159[1] + 0.5 p[1] * z_p_base_159[1]
   + 0.5 p[2] * z_p_base_159[1] + 0.5 p[3] * z_p_base_159[1]
   - 0.5 p[4] * z_p_base_159[1] ] <= 0
 infeasibility_p_base_159_row_2_leq_0: 5e-06 z_p_base_159[2] + [
   - 0.5 p[0] * z_p_base_159[2] + 0.5 p[1] * z_p_base_159[2]
   - 0.5 p[2] * z_p_base_159[2] - 0.5 p[3] * z_p_base_159[2]
   + 0.5 p[4] * z_p_base_159[2] ] <= 0
 infeasibility_p_base_159_row_3_leq_0: 5e-06 z_p_base_159[3] + [
   p[2] * z_p_base_159[3] ] <= 0
 infeasibility_p_base_159_row_4_leq_0: 5e-06 z_p_base_159[4] + [
   - 0.5 p[0] * z_p_base_159[4] - 0.5 p[1] * z_p_base_159[4]
   - 0.5 p[2] * z_p_base_159[4] + 0.5 p[3] * z_p_base_159[4]
   + 0.5 p[4] * z_p_base_159[4] ] <= 0
 infeasibility_p_prime_base_159_row_0_leq_0: 5e-06 z_p_prime_base_159[0]
   + [ alpha[0] * z_p_prime_base_159[0] + p[0] * z_p_prime_base_159[0] ]
   <= 0
 infeasibility_p_prime_base_159_row_1_leq_0: 5e-06 z_p_prime_base_159[1]
   + [ 0.5 alpha[0] * z_p_prime_base_159[1]
   - 0.5 alpha[1] * z_p_prime_base_159[1]
   - 0.5 alpha[2] * z_p_prime_base_159[1]
   - 0.5 alpha[3] * z_p_prime_base_159[1]
   - 0.5 alpha[4] * z_p_prime_base_159[1]
   + 0.5 p[0] * z_p_prime_base_159[1] + 0.5 p[1] * z_p_prime_base_159[1]
   + 0.5 p[2] * z_p_prime_base_159[1] + 0.5 p[3] * z_p_prime_base_159[1]
   - 0.5 p[4] * z_p_prime_base_159[1] ] <= 0
 infeasibility_p_prime_base_159_row_2_leq_0: 5e-06 z_p_prime_base_159[2]
   + [ - 0.5 alpha[0] * z_p_prime_base_159[2]
   - 0.5 alpha[1] * z_p_prime_base_159[2]
   + 0.5 alpha[2] * z_p_prime_base_159[2]
   + 0.5 alpha[3] * z_p_prime_base_159[2]
   + 0.5 alpha[4] * z_p_prime_base_159[2]
   - 0.5 p[0] * z_p_prime_base_159[2] + 0.5 p[1] * z_p_prime_base_159[2]
   - 0.5 p[2] * z_p_prime_base_159[2] - 0.5 p[3] * z_p_prime_base_159[2]
   + 0.5 p[4] * z_p_prime_base_159[2] ] <= 0
 infeasibility_p_prime_base_159_row_3_leq_0: 5e-06 z_p_prime_base_159[3]
   + [ - alpha[2] * z_p_prime_base_159[3] + p[2] * z_p_prime_base_159[3] ]
   <= 0
 infeasibility_p_prime_base_159_row_4_leq_0: 5e-06 z_p_prime_base_159[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_159[4]
   + 0.5 alpha[1] * z_p_prime_base_159[4]
   + 0.5 alpha[2] * z_p_prime_base_159[4]
   - 0.5 alpha[3] * z_p_prime_base_159[4]
   + 0.5 alpha[4] * z_p_prime_base_159[4]
   - 0.5 p[0] * z_p_prime_base_159[4] - 0.5 p[1] * z_p_prime_base_159[4]
   - 0.5 p[2] * z_p_prime_base_159[4] + 0.5 p[3] * z_p_prime_base_159[4]
   + 0.5 p[4] * z_p_prime_base_159[4] ] <= 0
 infeasibility_p_base_160_row_0_leq_0: 5e-06 z_p_base_160[0] + [
   p[0] * z_p_base_160[0] ] <= 0
 infeasibility_p_base_160_row_1_leq_0: 5e-06 z_p_base_160[1] + [
   p[1] * z_p_base_160[1] ] <= 0
 infeasibility_p_base_160_row_2_leq_0: 5e-06 z_p_base_160[2] + [
   0.5 p[0] * z_p_base_160[2] - 0.5 p[1] * z_p_base_160[2]
   + 0.5 p[2] * z_p_base_160[2] + 0.5 p[3] * z_p_base_160[2]
   - 0.5 p[4] * z_p_base_160[2] ] <= 0
 infeasibility_p_base_160_row_3_leq_0: 5e-06 z_p_base_160[3] + [
   - 0.5 p[0] * z_p_base_160[3] + 0.5 p[1] * z_p_base_160[3]
   + 0.5 p[2] * z_p_base_160[3] - 0.5 p[3] * z_p_base_160[3]
   + 0.5 p[4] * z_p_base_160[3] ] <= 0
 infeasibility_p_base_160_row_4_leq_0: 5e-06 z_p_base_160[4] + [
   - 0.5 p[0] * z_p_base_160[4] - 0.5 p[1] * z_p_base_160[4]
   - 0.5 p[2] * z_p_base_160[4] + 0.5 p[3] * z_p_base_160[4]
   + 0.5 p[4] * z_p_base_160[4] ] <= 0
 infeasibility_p_prime_base_160_row_0_leq_0: 5e-06 z_p_prime_base_160[0]
   + [ alpha[0] * z_p_prime_base_160[0] + p[0] * z_p_prime_base_160[0] ]
   <= 0
 infeasibility_p_prime_base_160_row_1_leq_0: 5e-06 z_p_prime_base_160[1]
   + [ - alpha[1] * z_p_prime_base_160[1] + p[1] * z_p_prime_base_160[1] ]
   <= 0
 infeasibility_p_prime_base_160_row_2_leq_0: 5e-06 z_p_prime_base_160[2]
   + [ 0.5 alpha[0] * z_p_prime_base_160[2]
   + 0.5 alpha[1] * z_p_prime_base_160[2]
   - 0.5 alpha[2] * z_p_prime_base_160[2]
   - 0.5 alpha[3] * z_p_prime_base_160[2]
   - 0.5 alpha[4] * z_p_prime_base_160[2]
   + 0.5 p[0] * z_p_prime_base_160[2] - 0.5 p[1] * z_p_prime_base_160[2]
   + 0.5 p[2] * z_p_prime_base_160[2] + 0.5 p[3] * z_p_prime_base_160[2]
   - 0.5 p[4] * z_p_prime_base_160[2] ] <= 0
 infeasibility_p_prime_base_160_row_3_leq_0: 5e-06 z_p_prime_base_160[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_160[3]
   - 0.5 alpha[1] * z_p_prime_base_160[3]
   - 0.5 alpha[2] * z_p_prime_base_160[3]
   + 0.5 alpha[3] * z_p_prime_base_160[3]
   + 0.5 alpha[4] * z_p_prime_base_160[3]
   - 0.5 p[0] * z_p_prime_base_160[3] + 0.5 p[1] * z_p_prime_base_160[3]
   + 0.5 p[2] * z_p_prime_base_160[3] - 0.5 p[3] * z_p_prime_base_160[3]
   + 0.5 p[4] * z_p_prime_base_160[3] ] <= 0
 infeasibility_p_prime_base_160_row_4_leq_0: 5e-06 z_p_prime_base_160[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_160[4]
   + 0.5 alpha[1] * z_p_prime_base_160[4]
   + 0.5 alpha[2] * z_p_prime_base_160[4]
   - 0.5 alpha[3] * z_p_prime_base_160[4]
   + 0.5 alpha[4] * z_p_prime_base_160[4]
   - 0.5 p[0] * z_p_prime_base_160[4] - 0.5 p[1] * z_p_prime_base_160[4]
   - 0.5 p[2] * z_p_prime_base_160[4] + 0.5 p[3] * z_p_prime_base_160[4]
   + 0.5 p[4] * z_p_prime_base_160[4] ] <= 0
 infeasibility_p_base_161_row_0_leq_0: 5e-06 z_p_base_161[0] + [
   p[0] * z_p_base_161[0] ] <= 0
 infeasibility_p_base_161_row_1_leq_0: 5e-06 z_p_base_161[1] + [
   p[1] * z_p_base_161[1] ] <= 0
 infeasibility_p_base_161_row_2_leq_0: 5e-06 z_p_base_161[2] + [
   0.5 p[0] * z_p_base_161[2] + 0.5 p[1] * z_p_base_161[2]
   + 0.5 p[2] * z_p_base_161[2] + 0.5 p[3] * z_p_base_161[2]
   - 0.5 p[4] * z_p_base_161[2] ] <= 0
 infeasibility_p_base_161_row_3_leq_0: 5e-06 z_p_base_161[3] + [
   - 0.5 p[0] * z_p_base_161[3] - 0.5 p[1] * z_p_base_161[3]
   + 0.5 p[2] * z_p_base_161[3] - 0.5 p[3] * z_p_base_161[3]
   + 0.5 p[4] * z_p_base_161[3] ] <= 0
 infeasibility_p_base_161_row_4_leq_0: 5e-06 z_p_base_161[4] + [
   - 0.5 p[0] * z_p_base_161[4] - 0.5 p[1] * z_p_base_161[4]
   - 0.5 p[2] * z_p_base_161[4] + 0.5 p[3] * z_p_base_161[4]
   + 0.5 p[4] * z_p_base_161[4] ] <= 0
 infeasibility_p_prime_base_161_row_0_leq_0: 5e-06 z_p_prime_base_161[0]
   + [ alpha[0] * z_p_prime_base_161[0] + p[0] * z_p_prime_base_161[0] ]
   <= 0
 infeasibility_p_prime_base_161_row_1_leq_0: 5e-06 z_p_prime_base_161[1]
   + [ - alpha[1] * z_p_prime_base_161[1] + p[1] * z_p_prime_base_161[1] ]
   <= 0
 infeasibility_p_prime_base_161_row_2_leq_0: 5e-06 z_p_prime_base_161[2]
   + [ 0.5 alpha[0] * z_p_prime_base_161[2]
   - 0.5 alpha[1] * z_p_prime_base_161[2]
   - 0.5 alpha[2] * z_p_prime_base_161[2]
   - 0.5 alpha[3] * z_p_prime_base_161[2]
   - 0.5 alpha[4] * z_p_prime_base_161[2]
   + 0.5 p[0] * z_p_prime_base_161[2] + 0.5 p[1] * z_p_prime_base_161[2]
   + 0.5 p[2] * z_p_prime_base_161[2] + 0.5 p[3] * z_p_prime_base_161[2]
   - 0.5 p[4] * z_p_prime_base_161[2] ] <= 0
 infeasibility_p_prime_base_161_row_3_leq_0: 5e-06 z_p_prime_base_161[3]
   + [ - 0.5 alpha[0] * z_p_prime_base_161[3]
   + 0.5 alpha[1] * z_p_prime_base_161[3]
   - 0.5 alpha[2] * z_p_prime_base_161[3]
   + 0.5 alpha[3] * z_p_prime_base_161[3]
   + 0.5 alpha[4] * z_p_prime_base_161[3]
   - 0.5 p[0] * z_p_prime_base_161[3] - 0.5 p[1] * z_p_prime_base_161[3]
   + 0.5 p[2] * z_p_prime_base_161[3] - 0.5 p[3] * z_p_prime_base_161[3]
   + 0.5 p[4] * z_p_prime_base_161[3] ] <= 0
 infeasibility_p_prime_base_161_row_4_leq_0: 5e-06 z_p_prime_base_161[4]
   + [ - 0.5 alpha[0] * z_p_prime_base_161[4]
   + 0.5 alpha[1] * z_p_prime_base_161[4]
   + 0.5 alpha[2] * z_p_prime_base_161[4]
   - 0.5 alpha[3] * z_p_prime_base_161[4]
   + 0.5 alpha[4] * z_p_prime_base_161[4]
   - 0.5 p[0] * z_p_prime_base_161[4] - 0.5 p[1] * z_p_prime_base_161[4]
   - 0.5 p[2] * z_p_prime_base_161[4] + 0.5 p[3] * z_p_prime_base_161[4]
   + 0.5 p[4] * z_p_prime_base_161[4] ] <= 0
 infeasibility_p_base_B_prime_row_0_leq_0: 5e-06 z_p_base_B_prime[0] + [
   0.5 p[0] * z_p_base_B_prime[0] + 0.5 p[1] * z_p_base_B_prime[0]
   + 0.5 p[2] * z_p_base_B_prime[0] - 0.5 p[3] * z_p_base_B_prime[0]
   - 0.5 p[4] * z_p_base_B_prime[0] ] <= 0
 infeasibility_p_base_B_prime_row_1_leq_0: 5e-06 z_p_base_B_prime[1] + [
   0.5 p[0] * z_p_base_B_prime[1] - 0.5 p[1] * z_p_base_B_prime[1]
   - 0.5 p[2] * z_p_base_B_prime[1] + 0.5 p[3] * z_p_base_B_prime[1]
   + 0.5 p[4] * z_p_base_B_prime[1] ] <= 0
 infeasibility_p_base_B_prime_row_2_leq_0: 5e-06 z_p_base_B_prime[2] + [
   - 0.5 p[0] * z_p_base_B_prime[2] + 0.5 p[1] * z_p_base_B_prime[2]
   - 0.5 p[2] * z_p_base_B_prime[2] + 0.5 p[3] * z_p_base_B_prime[2]
   + 0.5 p[4] * z_p_base_B_prime[2] ] <= 0
 infeasibility_p_base_B_prime_row_3_leq_0: 5e-06 z_p_base_B_prime[3] + [
   p[3] * z_p_base_B_prime[3] ] <= 0
 infeasibility_p_base_B_prime_row_4_leq_0: 5e-06 z_p_base_B_prime[4] + [
   p[2] * z_p_base_B_prime[4] - p[3] * z_p_base_B_prime[4] ] <= 0
 infeasibility_p_prime_base_B_row_0_leq_0: 5e-06 z_p_prime_base_B[0] + [
   0.5 alpha[0] * z_p_prime_base_B[0] - 0.5 alpha[1] * z_p_prime_base_B[0]
   + 0.5 alpha[2] * z_p_prime_base_B[0]
   + 0.5 alpha[3] * z_p_prime_base_B[0]
   - 0.5 alpha[4] * z_p_prime_base_B[0] + 0.5 p[0] * z_p_prime_base_B[0]
   + 0.5 p[1] * z_p_prime_base_B[0] - 0.5 p[2] * z_p_prime_base_B[0]
   - 0.5 p[3] * z_p_prime_base_B[0] - 0.5 p[4] * z_p_prime_base_B[0] ]
   <= 0
 infeasibility_p_prime_base_B_row_1_leq_0: 5e-06 z_p_prime_base_B[1] + [
   0.5 alpha[0] * z_p_prime_base_B[1] + 0.5 alpha[1] * z_p_prime_base_B[1]
   - 0.5 alpha[2] * z_p_prime_base_B[1]
   + 0.5 alpha[3] * z_p_prime_base_B[1]
   - 0.5 alpha[4] * z_p_prime_base_B[1] + 0.5 p[0] * z_p_prime_base_B[1]
   - 0.5 p[1] * z_p_prime_base_B[1] + 0.5 p[2] * z_p_prime_base_B[1]
   - 0.5 p[3] * z_p_prime_base_B[1] - 0.5 p[4] * z_p_prime_base_B[1] ]
   <= 0
 infeasibility_p_prime_base_B_row_2_leq_0: 5e-06 z_p_prime_base_B[2] + [
   - alpha[3] * z_p_prime_base_B[2] + p[3] * z_p_prime_base_B[2] ] <= 0
 infeasibility_p_prime_base_B_row_3_leq_0: 5e-06 z_p_prime_base_B[3] + [
   alpha[4] * z_p_prime_base_B[3] + p[4] * z_p_prime_base_B[3] ] <= 0
 infeasibility_p_prime_base_B_row_4_leq_0: 5e-06 z_p_prime_base_B[4] + [
   - 0.5 alpha[0] * z_p_prime_base_B[4]
   - 0.5 alpha[1] * z_p_prime_base_B[4]
   - 0.5 alpha[2] * z_p_prime_base_B[4]
   - 0.5 alpha[3] * z_p_prime_base_B[4]
   + 0.5 alpha[4] * z_p_prime_base_B[4] - 0.5 p[0] * z_p_prime_base_B[4]
   + 0.5 p[1] * z_p_prime_base_B[4] + 0.5 p[2] * z_p_prime_base_B[4]
   + 0.5 p[3] * z_p_prime_base_B[4] + 0.5 p[4] * z_p_prime_base_B[4] ]
   <= 0
Bounds
 alpha[0] <= 1
 alpha[1] <= 1
 alpha[2] <= 1
 alpha[3] <= 1
 alpha[4] <= 1
 0.1 <= p[0] <= 0.9
 0.1 <= p[1] <= 0.9
 0.1 <= p[2] <= 0.9
 0.1 <= p[3] <= 0.9
 0.1 <= p[4] <= 0.9
Binaries
 z_p_base_1[0] z_p_base_1[1] z_p_base_1[2] z_p_base_1[3] z_p_base_1[4]
 z_p_prime_base_1[0] z_p_prime_base_1[1] z_p_prime_base_1[2]
 z_p_prime_base_1[3] z_p_prime_base_1[4] z_p_base_2[0] z_p_base_2[1]
 z_p_base_2[2] z_p_base_2[3] z_p_base_2[4] z_p_prime_base_2[0]
 z_p_prime_base_2[1] z_p_prime_base_2[2] z_p_prime_base_2[3]
 z_p_prime_base_2[4] z_p_base_3[0] z_p_base_3[1] z_p_base_3[2]
 z_p_base_3[3] z_p_base_3[4] z_p_prime_base_3[0] z_p_prime_base_3[1]
 z_p_prime_base_3[2] z_p_prime_base_3[3] z_p_prime_base_3[4] z_p_base_4[0]
 z_p_base_4[1] z_p_base_4[2] z_p_base_4[3] z_p_base_4[4]
 z_p_prime_base_4[0] z_p_prime_base_4[1] z_p_prime_base_4[2]
 z_p_prime_base_4[3] z_p_prime_base_4[4] z_p_base_5[0] z_p_base_5[1]
 z_p_base_5[2] z_p_base_5[3] z_p_base_5[4] z_p_prime_base_5[0]
 z_p_prime_base_5[1] z_p_prime_base_5[2] z_p_prime_base_5[3]
 z_p_prime_base_5[4] z_p_base_6[0] z_p_base_6[1] z_p_base_6[2]
 z_p_base_6[3] z_p_base_6[4] z_p_prime_base_6[0] z_p_prime_base_6[1]
 z_p_prime_base_6[2] z_p_prime_base_6[3] z_p_prime_base_6[4] z_p_base_7[0]
 z_p_base_7[1] z_p_base_7[2] z_p_base_7[3] z_p_base_7[4]
 z_p_prime_base_7[0] z_p_prime_base_7[1] z_p_prime_base_7[2]
 z_p_prime_base_7[3] z_p_prime_base_7[4] z_p_base_8[0] z_p_base_8[1]
 z_p_base_8[2] z_p_base_8[3] z_p_base_8[4] z_p_prime_base_8[0]
 z_p_prime_base_8[1] z_p_prime_base_8[2] z_p_prime_base_8[3]
 z_p_prime_base_8[4] z_p_base_9[0] z_p_base_9[1] z_p_base_9[2]
 z_p_base_9[3] z_p_base_9[4] z_p_prime_base_9[0] z_p_prime_base_9[1]
 z_p_prime_base_9[2] z_p_prime_base_9[3] z_p_prime_base_9[4] z_p_base_10[0]
 z_p_base_10[1] z_p_base_10[2] z_p_base_10[3] z_p_base_10[4]
 z_p_prime_base_10[0] z_p_prime_base_10[1] z_p_prime_base_10[2]
 z_p_prime_base_10[3] z_p_prime_base_10[4] z_p_base_11[0] z_p_base_11[1]
 z_p_base_11[2] z_p_base_11[3] z_p_base_11[4] z_p_prime_base_11[0]
 z_p_prime_base_11[1] z_p_prime_base_11[2] z_p_prime_base_11[3]
 z_p_prime_base_11[4] z_p_base_12[0] z_p_base_12[1] z_p_base_12[2]
 z_p_base_12[3] z_p_base_12[4] z_p_prime_base_12[0] z_p_prime_base_12[1]
 z_p_prime_base_12[2] z_p_prime_base_12[3] z_p_prime_base_12[4]
 z_p_base_13[0] z_p_base_13[1] z_p_base_13[2] z_p_base_13[3] z_p_base_13[4]
 z_p_prime_base_13[0] z_p_prime_base_13[1] z_p_prime_base_13[2]
 z_p_prime_base_13[3] z_p_prime_base_13[4] z_p_base_14[0] z_p_base_14[1]
 z_p_base_14[2] z_p_base_14[3] z_p_base_14[4] z_p_prime_base_14[0]
 z_p_prime_base_14[1] z_p_prime_base_14[2] z_p_prime_base_14[3]
 z_p_prime_base_14[4] z_p_base_15[0] z_p_base_15[1] z_p_base_15[2]
 z_p_base_15[3] z_p_base_15[4] z_p_prime_base_15[0] z_p_prime_base_15[1]
 z_p_prime_base_15[2] z_p_prime_base_15[3] z_p_prime_base_15[4]
 z_p_base_16[0] z_p_base_16[1] z_p_base_16[2] z_p_base_16[3] z_p_base_16[4]
 z_p_prime_base_16[0] z_p_prime_base_16[1] z_p_prime_base_16[2]
 z_p_prime_base_16[3] z_p_prime_base_16[4] z_p_base_17[0] z_p_base_17[1]
 z_p_base_17[2] z_p_base_17[3] z_p_base_17[4] z_p_prime_base_17[0]
 z_p_prime_base_17[1] z_p_prime_base_17[2] z_p_prime_base_17[3]
 z_p_prime_base_17[4] z_p_base_18[0] z_p_base_18[1] z_p_base_18[2]
 z_p_base_18[3] z_p_base_18[4] z_p_prime_base_18[0] z_p_prime_base_18[1]
 z_p_prime_base_18[2] z_p_prime_base_18[3] z_p_prime_base_18[4]
 z_p_base_19[0] z_p_base_19[1] z_p_base_19[2] z_p_base_19[3] z_p_base_19[4]
 z_p_prime_base_19[0] z_p_prime_base_19[1] z_p_prime_base_19[2]
 z_p_prime_base_19[3] z_p_prime_base_19[4] z_p_base_20[0] z_p_base_20[1]
 z_p_base_20[2] z_p_base_20[3] z_p_base_20[4] z_p_prime_base_20[0]
 z_p_prime_base_20[1] z_p_prime_base_20[2] z_p_prime_base_20[3]
 z_p_prime_base_20[4] z_p_base_21[0] z_p_base_21[1] z_p_base_21[2]
 z_p_base_21[3] z_p_base_21[4] z_p_prime_base_21[0] z_p_prime_base_21[1]
 z_p_prime_base_21[2] z_p_prime_base_21[3] z_p_prime_base_21[4]
 z_p_base_22[0] z_p_base_22[1] z_p_base_22[2] z_p_base_22[3] z_p_base_22[4]
 z_p_prime_base_22[0] z_p_prime_base_22[1] z_p_prime_base_22[2]
 z_p_prime_base_22[3] z_p_prime_base_22[4] z_p_base_23[0] z_p_base_23[1]
 z_p_base_23[2] z_p_base_23[3] z_p_base_23[4] z_p_prime_base_23[0]
 z_p_prime_base_23[1] z_p_prime_base_23[2] z_p_prime_base_23[3]
 z_p_prime_base_23[4] z_p_base_24[0] z_p_base_24[1] z_p_base_24[2]
 z_p_base_24[3] z_p_base_24[4] z_p_prime_base_24[0] z_p_prime_base_24[1]
 z_p_prime_base_24[2] z_p_prime_base_24[3] z_p_prime_base_24[4]
 z_p_base_25[0] z_p_base_25[1] z_p_base_25[2] z_p_base_25[3] z_p_base_25[4]
 z_p_prime_base_25[0] z_p_prime_base_25[1] z_p_prime_base_25[2]
 z_p_prime_base_25[3] z_p_prime_base_25[4] z_p_base_26[0] z_p_base_26[1]
 z_p_base_26[2] z_p_base_26[3] z_p_base_26[4] z_p_prime_base_26[0]
 z_p_prime_base_26[1] z_p_prime_base_26[2] z_p_prime_base_26[3]
 z_p_prime_base_26[4] z_p_base_27[0] z_p_base_27[1] z_p_base_27[2]
 z_p_base_27[3] z_p_base_27[4] z_p_prime_base_27[0] z_p_prime_base_27[1]
 z_p_prime_base_27[2] z_p_prime_base_27[3] z_p_prime_base_27[4]
 z_p_base_28[0] z_p_base_28[1] z_p_base_28[2] z_p_base_28[3] z_p_base_28[4]
 z_p_prime_base_28[0] z_p_prime_base_28[1] z_p_prime_base_28[2]
 z_p_prime_base_28[3] z_p_prime_base_28[4] z_p_base_29[0] z_p_base_29[1]
 z_p_base_29[2] z_p_base_29[3] z_p_base_29[4] z_p_prime_base_29[0]
 z_p_prime_base_29[1] z_p_prime_base_29[2] z_p_prime_base_29[3]
 z_p_prime_base_29[4] z_p_base_30[0] z_p_base_30[1] z_p_base_30[2]
 z_p_base_30[3] z_p_base_30[4] z_p_prime_base_30[0] z_p_prime_base_30[1]
 z_p_prime_base_30[2] z_p_prime_base_30[3] z_p_prime_base_30[4]
 z_p_base_31[0] z_p_base_31[1] z_p_base_31[2] z_p_base_31[3] z_p_base_31[4]
 z_p_prime_base_31[0] z_p_prime_base_31[1] z_p_prime_base_31[2]
 z_p_prime_base_31[3] z_p_prime_base_31[4] z_p_base_32[0] z_p_base_32[1]
 z_p_base_32[2] z_p_base_32[3] z_p_base_32[4] z_p_prime_base_32[0]
 z_p_prime_base_32[1] z_p_prime_base_32[2] z_p_prime_base_32[3]
 z_p_prime_base_32[4] z_p_base_33[0] z_p_base_33[1] z_p_base_33[2]
 z_p_base_33[3] z_p_base_33[4] z_p_prime_base_33[0] z_p_prime_base_33[1]
 z_p_prime_base_33[2] z_p_prime_base_33[3] z_p_prime_base_33[4]
 z_p_base_34[0] z_p_base_34[1] z_p_base_34[2] z_p_base_34[3] z_p_base_34[4]
 z_p_prime_base_34[0] z_p_prime_base_34[1] z_p_prime_base_34[2]
 z_p_prime_base_34[3] z_p_prime_base_34[4] z_p_base_35[0] z_p_base_35[1]
 z_p_base_35[2] z_p_base_35[3] z_p_base_35[4] z_p_prime_base_35[0]
 z_p_prime_base_35[1] z_p_prime_base_35[2] z_p_prime_base_35[3]
 z_p_prime_base_35[4] z_p_base_36[0] z_p_base_36[1] z_p_base_36[2]
 z_p_base_36[3] z_p_base_36[4] z_p_prime_base_36[0] z_p_prime_base_36[1]
 z_p_prime_base_36[2] z_p_prime_base_36[3] z_p_prime_base_36[4]
 z_p_base_37[0] z_p_base_37[1] z_p_base_37[2] z_p_base_37[3] z_p_base_37[4]
 z_p_prime_base_37[0] z_p_prime_base_37[1] z_p_prime_base_37[2]
 z_p_prime_base_37[3] z_p_prime_base_37[4] z_p_base_38[0] z_p_base_38[1]
 z_p_base_38[2] z_p_base_38[3] z_p_base_38[4] z_p_prime_base_38[0]
 z_p_prime_base_38[1] z_p_prime_base_38[2] z_p_prime_base_38[3]
 z_p_prime_base_38[4] z_p_base_39[0] z_p_base_39[1] z_p_base_39[2]
 z_p_base_39[3] z_p_base_39[4] z_p_prime_base_39[0] z_p_prime_base_39[1]
 z_p_prime_base_39[2] z_p_prime_base_39[3] z_p_prime_base_39[4]
 z_p_base_40[0] z_p_base_40[1] z_p_base_40[2] z_p_base_40[3] z_p_base_40[4]
 z_p_prime_base_40[0] z_p_prime_base_40[1] z_p_prime_base_40[2]
 z_p_prime_base_40[3] z_p_prime_base_40[4] z_p_base_41[0] z_p_base_41[1]
 z_p_base_41[2] z_p_base_41[3] z_p_base_41[4] z_p_prime_base_41[0]
 z_p_prime_base_41[1] z_p_prime_base_41[2] z_p_prime_base_41[3]
 z_p_prime_base_41[4] z_p_base_42[0] z_p_base_42[1] z_p_base_42[2]
 z_p_base_42[3] z_p_base_42[4] z_p_prime_base_42[0] z_p_prime_base_42[1]
 z_p_prime_base_42[2] z_p_prime_base_42[3] z_p_prime_base_42[4]
 z_p_base_43[0] z_p_base_43[1] z_p_base_43[2] z_p_base_43[3] z_p_base_43[4]
 z_p_prime_base_43[0] z_p_prime_base_43[1] z_p_prime_base_43[2]
 z_p_prime_base_43[3] z_p_prime_base_43[4] z_p_base_44[0] z_p_base_44[1]
 z_p_base_44[2] z_p_base_44[3] z_p_base_44[4] z_p_prime_base_44[0]
 z_p_prime_base_44[1] z_p_prime_base_44[2] z_p_prime_base_44[3]
 z_p_prime_base_44[4] z_p_base_45[0] z_p_base_45[1] z_p_base_45[2]
 z_p_base_45[3] z_p_base_45[4] z_p_prime_base_45[0] z_p_prime_base_45[1]
 z_p_prime_base_45[2] z_p_prime_base_45[3] z_p_prime_base_45[4]
 z_p_base_46[0] z_p_base_46[1] z_p_base_46[2] z_p_base_46[3] z_p_base_46[4]
 z_p_prime_base_46[0] z_p_prime_base_46[1] z_p_prime_base_46[2]
 z_p_prime_base_46[3] z_p_prime_base_46[4] z_p_base_47[0] z_p_base_47[1]
 z_p_base_47[2] z_p_base_47[3] z_p_base_47[4] z_p_prime_base_47[0]
 z_p_prime_base_47[1] z_p_prime_base_47[2] z_p_prime_base_47[3]
 z_p_prime_base_47[4] z_p_base_48[0] z_p_base_48[1] z_p_base_48[2]
 z_p_base_48[3] z_p_base_48[4] z_p_prime_base_48[0] z_p_prime_base_48[1]
 z_p_prime_base_48[2] z_p_prime_base_48[3] z_p_prime_base_48[4]
 z_p_base_49[0] z_p_base_49[1] z_p_base_49[2] z_p_base_49[3] z_p_base_49[4]
 z_p_prime_base_49[0] z_p_prime_base_49[1] z_p_prime_base_49[2]
 z_p_prime_base_49[3] z_p_prime_base_49[4] z_p_base_50[0] z_p_base_50[1]
 z_p_base_50[2] z_p_base_50[3] z_p_base_50[4] z_p_prime_base_50[0]
 z_p_prime_base_50[1] z_p_prime_base_50[2] z_p_prime_base_50[3]
 z_p_prime_base_50[4] z_p_base_51[0] z_p_base_51[1] z_p_base_51[2]
 z_p_base_51[3] z_p_base_51[4] z_p_prime_base_51[0] z_p_prime_base_51[1]
 z_p_prime_base_51[2] z_p_prime_base_51[3] z_p_prime_base_51[4]
 z_p_base_52[0] z_p_base_52[1] z_p_base_52[2] z_p_base_52[3] z_p_base_52[4]
 z_p_prime_base_52[0] z_p_prime_base_52[1] z_p_prime_base_52[2]
 z_p_prime_base_52[3] z_p_prime_base_52[4] z_p_base_53[0] z_p_base_53[1]
 z_p_base_53[2] z_p_base_53[3] z_p_base_53[4] z_p_prime_base_53[0]
 z_p_prime_base_53[1] z_p_prime_base_53[2] z_p_prime_base_53[3]
 z_p_prime_base_53[4] z_p_base_54[0] z_p_base_54[1] z_p_base_54[2]
 z_p_base_54[3] z_p_base_54[4] z_p_prime_base_54[0] z_p_prime_base_54[1]
 z_p_prime_base_54[2] z_p_prime_base_54[3] z_p_prime_base_54[4]
 z_p_base_55[0] z_p_base_55[1] z_p_base_55[2] z_p_base_55[3] z_p_base_55[4]
 z_p_prime_base_55[0] z_p_prime_base_55[1] z_p_prime_base_55[2]
 z_p_prime_base_55[3] z_p_prime_base_55[4] z_p_base_56[0] z_p_base_56[1]
 z_p_base_56[2] z_p_base_56[3] z_p_base_56[4] z_p_prime_base_56[0]
 z_p_prime_base_56[1] z_p_prime_base_56[2] z_p_prime_base_56[3]
 z_p_prime_base_56[4] z_p_base_57[0] z_p_base_57[1] z_p_base_57[2]
 z_p_base_57[3] z_p_base_57[4] z_p_prime_base_57[0] z_p_prime_base_57[1]
 z_p_prime_base_57[2] z_p_prime_base_57[3] z_p_prime_base_57[4]
 z_p_base_58[0] z_p_base_58[1] z_p_base_58[2] z_p_base_58[3] z_p_base_58[4]
 z_p_prime_base_58[0] z_p_prime_base_58[1] z_p_prime_base_58[2]
 z_p_prime_base_58[3] z_p_prime_base_58[4] z_p_base_59[0] z_p_base_59[1]
 z_p_base_59[2] z_p_base_59[3] z_p_base_59[4] z_p_prime_base_59[0]
 z_p_prime_base_59[1] z_p_prime_base_59[2] z_p_prime_base_59[3]
 z_p_prime_base_59[4] z_p_base_60[0] z_p_base_60[1] z_p_base_60[2]
 z_p_base_60[3] z_p_base_60[4] z_p_prime_base_60[0] z_p_prime_base_60[1]
 z_p_prime_base_60[2] z_p_prime_base_60[3] z_p_prime_base_60[4]
 z_p_base_61[0] z_p_base_61[1] z_p_base_61[2] z_p_base_61[3] z_p_base_61[4]
 z_p_prime_base_61[0] z_p_prime_base_61[1] z_p_prime_base_61[2]
 z_p_prime_base_61[3] z_p_prime_base_61[4] z_p_base_62[0] z_p_base_62[1]
 z_p_base_62[2] z_p_base_62[3] z_p_base_62[4] z_p_prime_base_62[0]
 z_p_prime_base_62[1] z_p_prime_base_62[2] z_p_prime_base_62[3]
 z_p_prime_base_62[4] z_p_base_63[0] z_p_base_63[1] z_p_base_63[2]
 z_p_base_63[3] z_p_base_63[4] z_p_prime_base_63[0] z_p_prime_base_63[1]
 z_p_prime_base_63[2] z_p_prime_base_63[3] z_p_prime_base_63[4]
 z_p_base_64[0] z_p_base_64[1] z_p_base_64[2] z_p_base_64[3] z_p_base_64[4]
 z_p_prime_base_64[0] z_p_prime_base_64[1] z_p_prime_base_64[2]
 z_p_prime_base_64[3] z_p_prime_base_64[4] z_p_base_66[0] z_p_base_66[1]
 z_p_base_66[2] z_p_base_66[3] z_p_base_66[4] z_p_prime_base_66[0]
 z_p_prime_base_66[1] z_p_prime_base_66[2] z_p_prime_base_66[3]
 z_p_prime_base_66[4] z_p_base_67[0] z_p_base_67[1] z_p_base_67[2]
 z_p_base_67[3] z_p_base_67[4] z_p_prime_base_67[0] z_p_prime_base_67[1]
 z_p_prime_base_67[2] z_p_prime_base_67[3] z_p_prime_base_67[4]
 z_p_base_68[0] z_p_base_68[1] z_p_base_68[2] z_p_base_68[3] z_p_base_68[4]
 z_p_prime_base_68[0] z_p_prime_base_68[1] z_p_prime_base_68[2]
 z_p_prime_base_68[3] z_p_prime_base_68[4] z_p_base_69[0] z_p_base_69[1]
 z_p_base_69[2] z_p_base_69[3] z_p_base_69[4] z_p_prime_base_69[0]
 z_p_prime_base_69[1] z_p_prime_base_69[2] z_p_prime_base_69[3]
 z_p_prime_base_69[4] z_p_base_70[0] z_p_base_70[1] z_p_base_70[2]
 z_p_base_70[3] z_p_base_70[4] z_p_prime_base_70[0] z_p_prime_base_70[1]
 z_p_prime_base_70[2] z_p_prime_base_70[3] z_p_prime_base_70[4]
 z_p_base_71[0] z_p_base_71[1] z_p_base_71[2] z_p_base_71[3] z_p_base_71[4]
 z_p_prime_base_71[0] z_p_prime_base_71[1] z_p_prime_base_71[2]
 z_p_prime_base_71[3] z_p_prime_base_71[4] z_p_base_72[0] z_p_base_72[1]
 z_p_base_72[2] z_p_base_72[3] z_p_base_72[4] z_p_prime_base_72[0]
 z_p_prime_base_72[1] z_p_prime_base_72[2] z_p_prime_base_72[3]
 z_p_prime_base_72[4] z_p_base_73[0] z_p_base_73[1] z_p_base_73[2]
 z_p_base_73[3] z_p_base_73[4] z_p_prime_base_73[0] z_p_prime_base_73[1]
 z_p_prime_base_73[2] z_p_prime_base_73[3] z_p_prime_base_73[4]
 z_p_base_74[0] z_p_base_74[1] z_p_base_74[2] z_p_base_74[3] z_p_base_74[4]
 z_p_prime_base_74[0] z_p_prime_base_74[1] z_p_prime_base_74[2]
 z_p_prime_base_74[3] z_p_prime_base_74[4] z_p_base_75[0] z_p_base_75[1]
 z_p_base_75[2] z_p_base_75[3] z_p_base_75[4] z_p_prime_base_75[0]
 z_p_prime_base_75[1] z_p_prime_base_75[2] z_p_prime_base_75[3]
 z_p_prime_base_75[4] z_p_base_76[0] z_p_base_76[1] z_p_base_76[2]
 z_p_base_76[3] z_p_base_76[4] z_p_prime_base_76[0] z_p_prime_base_76[1]
 z_p_prime_base_76[2] z_p_prime_base_76[3] z_p_prime_base_76[4]
 z_p_base_77[0] z_p_base_77[1] z_p_base_77[2] z_p_base_77[3] z_p_base_77[4]
 z_p_prime_base_77[0] z_p_prime_base_77[1] z_p_prime_base_77[2]
 z_p_prime_base_77[3] z_p_prime_base_77[4] z_p_base_78[0] z_p_base_78[1]
 z_p_base_78[2] z_p_base_78[3] z_p_base_78[4] z_p_prime_base_78[0]
 z_p_prime_base_78[1] z_p_prime_base_78[2] z_p_prime_base_78[3]
 z_p_prime_base_78[4] z_p_base_79[0] z_p_base_79[1] z_p_base_79[2]
 z_p_base_79[3] z_p_base_79[4] z_p_prime_base_79[0] z_p_prime_base_79[1]
 z_p_prime_base_79[2] z_p_prime_base_79[3] z_p_prime_base_79[4]
 z_p_base_80[0] z_p_base_80[1] z_p_base_80[2] z_p_base_80[3] z_p_base_80[4]
 z_p_prime_base_80[0] z_p_prime_base_80[1] z_p_prime_base_80[2]
 z_p_prime_base_80[3] z_p_prime_base_80[4] z_p_base_81[0] z_p_base_81[1]
 z_p_base_81[2] z_p_base_81[3] z_p_base_81[4] z_p_prime_base_81[0]
 z_p_prime_base_81[1] z_p_prime_base_81[2] z_p_prime_base_81[3]
 z_p_prime_base_81[4] z_p_base_82[0] z_p_base_82[1] z_p_base_82[2]
 z_p_base_82[3] z_p_base_82[4] z_p_prime_base_82[0] z_p_prime_base_82[1]
 z_p_prime_base_82[2] z_p_prime_base_82[3] z_p_prime_base_82[4]
 z_p_base_83[0] z_p_base_83[1] z_p_base_83[2] z_p_base_83[3] z_p_base_83[4]
 z_p_prime_base_83[0] z_p_prime_base_83[1] z_p_prime_base_83[2]
 z_p_prime_base_83[3] z_p_prime_base_83[4] z_p_base_84[0] z_p_base_84[1]
 z_p_base_84[2] z_p_base_84[3] z_p_base_84[4] z_p_prime_base_84[0]
 z_p_prime_base_84[1] z_p_prime_base_84[2] z_p_prime_base_84[3]
 z_p_prime_base_84[4] z_p_base_85[0] z_p_base_85[1] z_p_base_85[2]
 z_p_base_85[3] z_p_base_85[4] z_p_prime_base_85[0] z_p_prime_base_85[1]
 z_p_prime_base_85[2] z_p_prime_base_85[3] z_p_prime_base_85[4]
 z_p_base_86[0] z_p_base_86[1] z_p_base_86[2] z_p_base_86[3] z_p_base_86[4]
 z_p_prime_base_86[0] z_p_prime_base_86[1] z_p_prime_base_86[2]
 z_p_prime_base_86[3] z_p_prime_base_86[4] z_p_base_87[0] z_p_base_87[1]
 z_p_base_87[2] z_p_base_87[3] z_p_base_87[4] z_p_prime_base_87[0]
 z_p_prime_base_87[1] z_p_prime_base_87[2] z_p_prime_base_87[3]
 z_p_prime_base_87[4] z_p_base_88[0] z_p_base_88[1] z_p_base_88[2]
 z_p_base_88[3] z_p_base_88[4] z_p_prime_base_88[0] z_p_prime_base_88[1]
 z_p_prime_base_88[2] z_p_prime_base_88[3] z_p_prime_base_88[4]
 z_p_base_89[0] z_p_base_89[1] z_p_base_89[2] z_p_base_89[3] z_p_base_89[4]
 z_p_prime_base_89[0] z_p_prime_base_89[1] z_p_prime_base_89[2]
 z_p_prime_base_89[3] z_p_prime_base_89[4] z_p_base_90[0] z_p_base_90[1]
 z_p_base_90[2] z_p_base_90[3] z_p_base_90[4] z_p_prime_base_90[0]
 z_p_prime_base_90[1] z_p_prime_base_90[2] z_p_prime_base_90[3]
 z_p_prime_base_90[4] z_p_base_91[0] z_p_base_91[1] z_p_base_91[2]
 z_p_base_91[3] z_p_base_91[4] z_p_prime_base_91[0] z_p_prime_base_91[1]
 z_p_prime_base_91[2] z_p_prime_base_91[3] z_p_prime_base_91[4]
 z_p_base_92[0] z_p_base_92[1] z_p_base_92[2] z_p_base_92[3] z_p_base_92[4]
 z_p_prime_base_92[0] z_p_prime_base_92[1] z_p_prime_base_92[2]
 z_p_prime_base_92[3] z_p_prime_base_92[4] z_p_base_93[0] z_p_base_93[1]
 z_p_base_93[2] z_p_base_93[3] z_p_base_93[4] z_p_prime_base_93[0]
 z_p_prime_base_93[1] z_p_prime_base_93[2] z_p_prime_base_93[3]
 z_p_prime_base_93[4] z_p_base_94[0] z_p_base_94[1] z_p_base_94[2]
 z_p_base_94[3] z_p_base_94[4] z_p_prime_base_94[0] z_p_prime_base_94[1]
 z_p_prime_base_94[2] z_p_prime_base_94[3] z_p_prime_base_94[4]
 z_p_base_95[0] z_p_base_95[1] z_p_base_95[2] z_p_base_95[3] z_p_base_95[4]
 z_p_prime_base_95[0] z_p_prime_base_95[1] z_p_prime_base_95[2]
 z_p_prime_base_95[3] z_p_prime_base_95[4] z_p_base_96[0] z_p_base_96[1]
 z_p_base_96[2] z_p_base_96[3] z_p_base_96[4] z_p_prime_base_96[0]
 z_p_prime_base_96[1] z_p_prime_base_96[2] z_p_prime_base_96[3]
 z_p_prime_base_96[4] z_p_base_97[0] z_p_base_97[1] z_p_base_97[2]
 z_p_base_97[3] z_p_base_97[4] z_p_prime_base_97[0] z_p_prime_base_97[1]
 z_p_prime_base_97[2] z_p_prime_base_97[3] z_p_prime_base_97[4]
 z_p_base_98[0] z_p_base_98[1] z_p_base_98[2] z_p_base_98[3] z_p_base_98[4]
 z_p_prime_base_98[0] z_p_prime_base_98[1] z_p_prime_base_98[2]
 z_p_prime_base_98[3] z_p_prime_base_98[4] z_p_base_99[0] z_p_base_99[1]
 z_p_base_99[2] z_p_base_99[3] z_p_base_99[4] z_p_prime_base_99[0]
 z_p_prime_base_99[1] z_p_prime_base_99[2] z_p_prime_base_99[3]
 z_p_prime_base_99[4] z_p_base_100[0] z_p_base_100[1] z_p_base_100[2]
 z_p_base_100[3] z_p_base_100[4] z_p_prime_base_100[0]
 z_p_prime_base_100[1] z_p_prime_base_100[2] z_p_prime_base_100[3]
 z_p_prime_base_100[4] z_p_base_101[0] z_p_base_101[1] z_p_base_101[2]
 z_p_base_101[3] z_p_base_101[4] z_p_prime_base_101[0]
 z_p_prime_base_101[1] z_p_prime_base_101[2] z_p_prime_base_101[3]
 z_p_prime_base_101[4] z_p_base_102[0] z_p_base_102[1] z_p_base_102[2]
 z_p_base_102[3] z_p_base_102[4] z_p_prime_base_102[0]
 z_p_prime_base_102[1] z_p_prime_base_102[2] z_p_prime_base_102[3]
 z_p_prime_base_102[4] z_p_base_103[0] z_p_base_103[1] z_p_base_103[2]
 z_p_base_103[3] z_p_base_103[4] z_p_prime_base_103[0]
 z_p_prime_base_103[1] z_p_prime_base_103[2] z_p_prime_base_103[3]
 z_p_prime_base_103[4] z_p_base_104[0] z_p_base_104[1] z_p_base_104[2]
 z_p_base_104[3] z_p_base_104[4] z_p_prime_base_104[0]
 z_p_prime_base_104[1] z_p_prime_base_104[2] z_p_prime_base_104[3]
 z_p_prime_base_104[4] z_p_base_105[0] z_p_base_105[1] z_p_base_105[2]
 z_p_base_105[3] z_p_base_105[4] z_p_prime_base_105[0]
 z_p_prime_base_105[1] z_p_prime_base_105[2] z_p_prime_base_105[3]
 z_p_prime_base_105[4] z_p_base_106[0] z_p_base_106[1] z_p_base_106[2]
 z_p_base_106[3] z_p_base_106[4] z_p_prime_base_106[0]
 z_p_prime_base_106[1] z_p_prime_base_106[2] z_p_prime_base_106[3]
 z_p_prime_base_106[4] z_p_base_107[0] z_p_base_107[1] z_p_base_107[2]
 z_p_base_107[3] z_p_base_107[4] z_p_prime_base_107[0]
 z_p_prime_base_107[1] z_p_prime_base_107[2] z_p_prime_base_107[3]
 z_p_prime_base_107[4] z_p_base_108[0] z_p_base_108[1] z_p_base_108[2]
 z_p_base_108[3] z_p_base_108[4] z_p_prime_base_108[0]
 z_p_prime_base_108[1] z_p_prime_base_108[2] z_p_prime_base_108[3]
 z_p_prime_base_108[4] z_p_base_109[0] z_p_base_109[1] z_p_base_109[2]
 z_p_base_109[3] z_p_base_109[4] z_p_prime_base_109[0]
 z_p_prime_base_109[1] z_p_prime_base_109[2] z_p_prime_base_109[3]
 z_p_prime_base_109[4] z_p_base_110[0] z_p_base_110[1] z_p_base_110[2]
 z_p_base_110[3] z_p_base_110[4] z_p_prime_base_110[0]
 z_p_prime_base_110[1] z_p_prime_base_110[2] z_p_prime_base_110[3]
 z_p_prime_base_110[4] z_p_base_111[0] z_p_base_111[1] z_p_base_111[2]
 z_p_base_111[3] z_p_base_111[4] z_p_prime_base_111[0]
 z_p_prime_base_111[1] z_p_prime_base_111[2] z_p_prime_base_111[3]
 z_p_prime_base_111[4] z_p_base_112[0] z_p_base_112[1] z_p_base_112[2]
 z_p_base_112[3] z_p_base_112[4] z_p_prime_base_112[0]
 z_p_prime_base_112[1] z_p_prime_base_112[2] z_p_prime_base_112[3]
 z_p_prime_base_112[4] z_p_base_113[0] z_p_base_113[1] z_p_base_113[2]
 z_p_base_113[3] z_p_base_113[4] z_p_prime_base_113[0]
 z_p_prime_base_113[1] z_p_prime_base_113[2] z_p_prime_base_113[3]
 z_p_prime_base_113[4] z_p_base_114[0] z_p_base_114[1] z_p_base_114[2]
 z_p_base_114[3] z_p_base_114[4] z_p_prime_base_114[0]
 z_p_prime_base_114[1] z_p_prime_base_114[2] z_p_prime_base_114[3]
 z_p_prime_base_114[4] z_p_base_115[0] z_p_base_115[1] z_p_base_115[2]
 z_p_base_115[3] z_p_base_115[4] z_p_prime_base_115[0]
 z_p_prime_base_115[1] z_p_prime_base_115[2] z_p_prime_base_115[3]
 z_p_prime_base_115[4] z_p_base_116[0] z_p_base_116[1] z_p_base_116[2]
 z_p_base_116[3] z_p_base_116[4] z_p_prime_base_116[0]
 z_p_prime_base_116[1] z_p_prime_base_116[2] z_p_prime_base_116[3]
 z_p_prime_base_116[4] z_p_base_117[0] z_p_base_117[1] z_p_base_117[2]
 z_p_base_117[3] z_p_base_117[4] z_p_prime_base_117[0]
 z_p_prime_base_117[1] z_p_prime_base_117[2] z_p_prime_base_117[3]
 z_p_prime_base_117[4] z_p_base_118[0] z_p_base_118[1] z_p_base_118[2]
 z_p_base_118[3] z_p_base_118[4] z_p_prime_base_118[0]
 z_p_prime_base_118[1] z_p_prime_base_118[2] z_p_prime_base_118[3]
 z_p_prime_base_118[4] z_p_base_119[0] z_p_base_119[1] z_p_base_119[2]
 z_p_base_119[3] z_p_base_119[4] z_p_prime_base_119[0]
 z_p_prime_base_119[1] z_p_prime_base_119[2] z_p_prime_base_119[3]
 z_p_prime_base_119[4] z_p_base_120[0] z_p_base_120[1] z_p_base_120[2]
 z_p_base_120[3] z_p_base_120[4] z_p_prime_base_120[0]
 z_p_prime_base_120[1] z_p_prime_base_120[2] z_p_prime_base_120[3]
 z_p_prime_base_120[4] z_p_base_121[0] z_p_base_121[1] z_p_base_121[2]
 z_p_base_121[3] z_p_base_121[4] z_p_prime_base_121[0]
 z_p_prime_base_121[1] z_p_prime_base_121[2] z_p_prime_base_121[3]
 z_p_prime_base_121[4] z_p_base_122[0] z_p_base_122[1] z_p_base_122[2]
 z_p_base_122[3] z_p_base_122[4] z_p_prime_base_122[0]
 z_p_prime_base_122[1] z_p_prime_base_122[2] z_p_prime_base_122[3]
 z_p_prime_base_122[4] z_p_base_123[0] z_p_base_123[1] z_p_base_123[2]
 z_p_base_123[3] z_p_base_123[4] z_p_prime_base_123[0]
 z_p_prime_base_123[1] z_p_prime_base_123[2] z_p_prime_base_123[3]
 z_p_prime_base_123[4] z_p_base_124[0] z_p_base_124[1] z_p_base_124[2]
 z_p_base_124[3] z_p_base_124[4] z_p_prime_base_124[0]
 z_p_prime_base_124[1] z_p_prime_base_124[2] z_p_prime_base_124[3]
 z_p_prime_base_124[4] z_p_base_125[0] z_p_base_125[1] z_p_base_125[2]
 z_p_base_125[3] z_p_base_125[4] z_p_prime_base_125[0]
 z_p_prime_base_125[1] z_p_prime_base_125[2] z_p_prime_base_125[3]
 z_p_prime_base_125[4] z_p_base_126[0] z_p_base_126[1] z_p_base_126[2]
 z_p_base_126[3] z_p_base_126[4] z_p_prime_base_126[0]
 z_p_prime_base_126[1] z_p_prime_base_126[2] z_p_prime_base_126[3]
 z_p_prime_base_126[4] z_p_base_127[0] z_p_base_127[1] z_p_base_127[2]
 z_p_base_127[3] z_p_base_127[4] z_p_prime_base_127[0]
 z_p_prime_base_127[1] z_p_prime_base_127[2] z_p_prime_base_127[3]
 z_p_prime_base_127[4] z_p_base_128[0] z_p_base_128[1] z_p_base_128[2]
 z_p_base_128[3] z_p_base_128[4] z_p_prime_base_128[0]
 z_p_prime_base_128[1] z_p_prime_base_128[2] z_p_prime_base_128[3]
 z_p_prime_base_128[4] z_p_base_129[0] z_p_base_129[1] z_p_base_129[2]
 z_p_base_129[3] z_p_base_129[4] z_p_prime_base_129[0]
 z_p_prime_base_129[1] z_p_prime_base_129[2] z_p_prime_base_129[3]
 z_p_prime_base_129[4] z_p_base_130[0] z_p_base_130[1] z_p_base_130[2]
 z_p_base_130[3] z_p_base_130[4] z_p_prime_base_130[0]
 z_p_prime_base_130[1] z_p_prime_base_130[2] z_p_prime_base_130[3]
 z_p_prime_base_130[4] z_p_base_131[0] z_p_base_131[1] z_p_base_131[2]
 z_p_base_131[3] z_p_base_131[4] z_p_prime_base_131[0]
 z_p_prime_base_131[1] z_p_prime_base_131[2] z_p_prime_base_131[3]
 z_p_prime_base_131[4] z_p_base_132[0] z_p_base_132[1] z_p_base_132[2]
 z_p_base_132[3] z_p_base_132[4] z_p_prime_base_132[0]
 z_p_prime_base_132[1] z_p_prime_base_132[2] z_p_prime_base_132[3]
 z_p_prime_base_132[4] z_p_base_133[0] z_p_base_133[1] z_p_base_133[2]
 z_p_base_133[3] z_p_base_133[4] z_p_prime_base_133[0]
 z_p_prime_base_133[1] z_p_prime_base_133[2] z_p_prime_base_133[3]
 z_p_prime_base_133[4] z_p_base_134[0] z_p_base_134[1] z_p_base_134[2]
 z_p_base_134[3] z_p_base_134[4] z_p_prime_base_134[0]
 z_p_prime_base_134[1] z_p_prime_base_134[2] z_p_prime_base_134[3]
 z_p_prime_base_134[4] z_p_base_135[0] z_p_base_135[1] z_p_base_135[2]
 z_p_base_135[3] z_p_base_135[4] z_p_prime_base_135[0]
 z_p_prime_base_135[1] z_p_prime_base_135[2] z_p_prime_base_135[3]
 z_p_prime_base_135[4] z_p_base_136[0] z_p_base_136[1] z_p_base_136[2]
 z_p_base_136[3] z_p_base_136[4] z_p_prime_base_136[0]
 z_p_prime_base_136[1] z_p_prime_base_136[2] z_p_prime_base_136[3]
 z_p_prime_base_136[4] z_p_base_137[0] z_p_base_137[1] z_p_base_137[2]
 z_p_base_137[3] z_p_base_137[4] z_p_prime_base_137[0]
 z_p_prime_base_137[1] z_p_prime_base_137[2] z_p_prime_base_137[3]
 z_p_prime_base_137[4] z_p_base_138[0] z_p_base_138[1] z_p_base_138[2]
 z_p_base_138[3] z_p_base_138[4] z_p_prime_base_138[0]
 z_p_prime_base_138[1] z_p_prime_base_138[2] z_p_prime_base_138[3]
 z_p_prime_base_138[4] z_p_base_139[0] z_p_base_139[1] z_p_base_139[2]
 z_p_base_139[3] z_p_base_139[4] z_p_prime_base_139[0]
 z_p_prime_base_139[1] z_p_prime_base_139[2] z_p_prime_base_139[3]
 z_p_prime_base_139[4] z_p_base_140[0] z_p_base_140[1] z_p_base_140[2]
 z_p_base_140[3] z_p_base_140[4] z_p_prime_base_140[0]
 z_p_prime_base_140[1] z_p_prime_base_140[2] z_p_prime_base_140[3]
 z_p_prime_base_140[4] z_p_base_141[0] z_p_base_141[1] z_p_base_141[2]
 z_p_base_141[3] z_p_base_141[4] z_p_prime_base_141[0]
 z_p_prime_base_141[1] z_p_prime_base_141[2] z_p_prime_base_141[3]
 z_p_prime_base_141[4] z_p_base_142[0] z_p_base_142[1] z_p_base_142[2]
 z_p_base_142[3] z_p_base_142[4] z_p_prime_base_142[0]
 z_p_prime_base_142[1] z_p_prime_base_142[2] z_p_prime_base_142[3]
 z_p_prime_base_142[4] z_p_base_143[0] z_p_base_143[1] z_p_base_143[2]
 z_p_base_143[3] z_p_base_143[4] z_p_prime_base_143[0]
 z_p_prime_base_143[1] z_p_prime_base_143[2] z_p_prime_base_143[3]
 z_p_prime_base_143[4] z_p_base_144[0] z_p_base_144[1] z_p_base_144[2]
 z_p_base_144[3] z_p_base_144[4] z_p_prime_base_144[0]
 z_p_prime_base_144[1] z_p_prime_base_144[2] z_p_prime_base_144[3]
 z_p_prime_base_144[4] z_p_base_145[0] z_p_base_145[1] z_p_base_145[2]
 z_p_base_145[3] z_p_base_145[4] z_p_prime_base_145[0]
 z_p_prime_base_145[1] z_p_prime_base_145[2] z_p_prime_base_145[3]
 z_p_prime_base_145[4] z_p_base_146[0] z_p_base_146[1] z_p_base_146[2]
 z_p_base_146[3] z_p_base_146[4] z_p_prime_base_146[0]
 z_p_prime_base_146[1] z_p_prime_base_146[2] z_p_prime_base_146[3]
 z_p_prime_base_146[4] z_p_base_147[0] z_p_base_147[1] z_p_base_147[2]
 z_p_base_147[3] z_p_base_147[4] z_p_prime_base_147[0]
 z_p_prime_base_147[1] z_p_prime_base_147[2] z_p_prime_base_147[3]
 z_p_prime_base_147[4] z_p_base_148[0] z_p_base_148[1] z_p_base_148[2]
 z_p_base_148[3] z_p_base_148[4] z_p_prime_base_148[0]
 z_p_prime_base_148[1] z_p_prime_base_148[2] z_p_prime_base_148[3]
 z_p_prime_base_148[4] z_p_base_149[0] z_p_base_149[1] z_p_base_149[2]
 z_p_base_149[3] z_p_base_149[4] z_p_prime_base_149[0]
 z_p_prime_base_149[1] z_p_prime_base_149[2] z_p_prime_base_149[3]
 z_p_prime_base_149[4] z_p_base_150[0] z_p_base_150[1] z_p_base_150[2]
 z_p_base_150[3] z_p_base_150[4] z_p_prime_base_150[0]
 z_p_prime_base_150[1] z_p_prime_base_150[2] z_p_prime_base_150[3]
 z_p_prime_base_150[4] z_p_base_151[0] z_p_base_151[1] z_p_base_151[2]
 z_p_base_151[3] z_p_base_151[4] z_p_prime_base_151[0]
 z_p_prime_base_151[1] z_p_prime_base_151[2] z_p_prime_base_151[3]
 z_p_prime_base_151[4] z_p_base_152[0] z_p_base_152[1] z_p_base_152[2]
 z_p_base_152[3] z_p_base_152[4] z_p_prime_base_152[0]
 z_p_prime_base_152[1] z_p_prime_base_152[2] z_p_prime_base_152[3]
 z_p_prime_base_152[4] z_p_base_153[0] z_p_base_153[1] z_p_base_153[2]
 z_p_base_153[3] z_p_base_153[4] z_p_prime_base_153[0]
 z_p_prime_base_153[1] z_p_prime_base_153[2] z_p_prime_base_153[3]
 z_p_prime_base_153[4] z_p_base_154[0] z_p_base_154[1] z_p_base_154[2]
 z_p_base_154[3] z_p_base_154[4] z_p_prime_base_154[0]
 z_p_prime_base_154[1] z_p_prime_base_154[2] z_p_prime_base_154[3]
 z_p_prime_base_154[4] z_p_base_155[0] z_p_base_155[1] z_p_base_155[2]
 z_p_base_155[3] z_p_base_155[4] z_p_prime_base_155[0]
 z_p_prime_base_155[1] z_p_prime_base_155[2] z_p_prime_base_155[3]
 z_p_prime_base_155[4] z_p_base_156[0] z_p_base_156[1] z_p_base_156[2]
 z_p_base_156[3] z_p_base_156[4] z_p_prime_base_156[0]
 z_p_prime_base_156[1] z_p_prime_base_156[2] z_p_prime_base_156[3]
 z_p_prime_base_156[4] z_p_base_157[0] z_p_base_157[1] z_p_base_157[2]
 z_p_base_157[3] z_p_base_157[4] z_p_prime_base_157[0]
 z_p_prime_base_157[1] z_p_prime_base_157[2] z_p_prime_base_157[3]
 z_p_prime_base_157[4] z_p_base_158[0] z_p_base_158[1] z_p_base_158[2]
 z_p_base_158[3] z_p_base_158[4] z_p_prime_base_158[0]
 z_p_prime_base_158[1] z_p_prime_base_158[2] z_p_prime_base_158[3]
 z_p_prime_base_158[4] z_p_base_159[0] z_p_base_159[1] z_p_base_159[2]
 z_p_base_159[3] z_p_base_159[4] z_p_prime_base_159[0]
 z_p_prime_base_159[1] z_p_prime_base_159[2] z_p_prime_base_159[3]
 z_p_prime_base_159[4] z_p_base_160[0] z_p_base_160[1] z_p_base_160[2]
 z_p_base_160[3] z_p_base_160[4] z_p_prime_base_160[0]
 z_p_prime_base_160[1] z_p_prime_base_160[2] z_p_prime_base_160[3]
 z_p_prime_base_160[4] z_p_base_161[0] z_p_base_161[1] z_p_base_161[2]
 z_p_base_161[3] z_p_base_161[4] z_p_prime_base_161[0]
 z_p_prime_base_161[1] z_p_prime_base_161[2] z_p_prime_base_161[3]
 z_p_prime_base_161[4] z_p_base_B_prime[0] z_p_base_B_prime[1]
 z_p_base_B_prime[2] z_p_base_B_prime[3] z_p_base_B_prime[4]
 z_p_prime_base_B[0] z_p_prime_base_B[1] z_p_prime_base_B[2]
 z_p_prime_base_B[3] z_p_prime_base_B[4]
End
