function [Cells_AB_sum_t_index, Cells_AB_sum_index] = fxn_MPPCA_ensemble_connect_cal(temp_ca_data, ensemble_set_A_IDs, ensemble_set_B_IDs, binarize_thr)
%%
temp_ca_A = temp_ca_data(:, ensemble_set_A_IDs);
temp_ca_B = temp_ca_data(:, ensemble_set_B_IDs);

temp_ca_A_logi = temp_ca_A > binarize_thr;
temp_ca_B_logi = temp_ca_B > binarize_thr;

Cells_A_sum = sum(temp_ca_A_logi, 2);
Cells_B_sum = sum(temp_ca_B_logi, 2);

for i = 1:size(Cells_A_sum,1)
Cells_AB_sum_t(i,1) = Cells_A_sum(i)*Cells_B_sum(i);
end

Cells_AB_sum_t_index = Cells_AB_sum_t / (size(temp_ca_A,2)*size(temp_ca_B,2));
Cells_AB_sum_index   = sum(Cells_AB_sum_t_index)/size(temp_ca_A,1);
%%
end