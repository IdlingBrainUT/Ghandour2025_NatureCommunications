%% MPPCA cell-cell pairwise interation 
function [result_connect] = fxn_MPPCA_ensemble_connectivity(result_data_cell, result_MPPCA, session_num_assess, ensemble_set_A, ensemble_set_B, shuffle_mode, binarize_thr)
%% comment
% 211118: 1st ver
%% debug info
% load('pairwise_demo.mat');
% session_num_assess        = 1 ; % Caution: Don't select same sassion to reference session, logically failed.
% binarize_thr              = prms_reactivation_SD_thr ; % same to red-highlighted thresholding.
% ensemble_set_A            = [ 1 ] ; % input pattern num
% ensemble_set_B            = [ 5 ] ; % input pattern num
% shuffle_mode              = 1;      % 1:circshift, 2:randperm
%% system parameter
shift = 1; % system parameter, don't change.

%% cell sorting
temp_ca_data = result_data_cell{session_num_assess+shift,7};

ensemble_set_A_IDs = result_MPPCA.neuron_sig_IDs{1, ensemble_set_A};
ensemble_set_B_IDs = result_MPPCA.neuron_sig_IDs{1, ensemble_set_B};

% for venn diagram
IDs_both_AB = intersect(ensemble_set_A_IDs,ensemble_set_B_IDs);
IDs_spec_A  = setdiff(ensemble_set_A_IDs,ensemble_set_B_IDs);
IDs_spec_B  = setdiff(ensemble_set_B_IDs,ensemble_set_A_IDs);

%   A = [300 200]; I = 150;
figure('Position', [1300 800 280 150]);
fxn_venn([numel(IDs_spec_A), numel(IDs_spec_B)], numel(IDs_both_AB))
xticklabels([]); yticklabels([]); title('Venn diagram')
%% Calculate connectivity
[Cells_AB_sum_t_index, Cells_AB_sum_index] = ...
  fxn_MPPCA_ensemble_connect_cal (temp_ca_data, ensemble_set_A_IDs, ensemble_set_B_IDs, binarize_thr);
%% Calculate shuffled connetctivity
% shuffle data generation
iteration = 1000;
% shuffle_mode = 2;
% shuffled_data = fxn_circ_shuffling(data_input,iteration);
temp_shuffled_data = fxn_data_shuffling(temp_ca_data, iteration, shuffle_mode);
%% shuffle data cal 
for s_i = 1:iteration
[Cells_AB_shuffle_t_index{s_i,1}, Cells_AB_shuffle_index{s_i,1}] = ... 
    fxn_MPPCA_ensemble_connect_cal(temp_shuffled_data{s_i,1} , ensemble_set_A_IDs, ensemble_set_B_IDs, binarize_thr);
end
disp(['Shuffled data iteration is now ', num2str(s_i)]);
%%
figure('Position', [50 50 600 200]) ;
subplot(121);
histogram(cell2mat(Cells_AB_shuffle_index))
hold on
xline(Cells_AB_sum_index,'r','Actual value','LineWidth',2)
%%
result_connect.connect_index = Cells_AB_sum_index;
result_connect.shuffle_index = Cells_AB_shuffle_index;
%%
end