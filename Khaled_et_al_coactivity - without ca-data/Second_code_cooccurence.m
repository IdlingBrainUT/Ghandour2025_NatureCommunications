%% less-fire -> bin -> ocurrence cal -> histogram

% function [result_sum_dur_range, data_for_hist1, data_for_hist2, i, ii, iii] = ...
%     func_bin_ocurrence_hist_triple(ca_filt_data, cut_frame_num, bin_frame_num, Cell_type_A_ID, Cell_type_B_ID, Cell_type_C_ID, cal_range_start, cal_range_end)

% ca_filt_data -> input data
% cut_frame_num = 1; % 0 = no cut, 1 = 50ms-frame-firing cut
% bin_frame_num = 10; % 500ms binning
% Cell_type_A_ID = [1,2,3,4,5];
% Cell_type_B_ID = [6,7,8,9,10,11,12,13,14,15];
% cal_range_start = 100; % input range in second scale
% cal_range_end   = 200; % input range in second scale


%% calculation double positive

ca_temp       = ca_filt_data; % input data
cut_frame_num = 0; % 0 = no cut, 1 = 50ms-frame-firing cut
bin_frame_num = 10; % 500ms binning
%% 
    cal_range_start = 2166; % input range in second scale
    cal_range_end   = 2226; % input range in second scale
    
%%  sample for THREE factors

        Cell_type_A_ID = [4,7,22,58,76,92]; % input Event1-Ai12 cell ID
        Cell_type_B_ID = [2,3,5,6,12,15,17,19,21,24,26,30,32,34,35,38,40,42,45,46,52,54,57,60,62,64,68,74,78,80,81,85,90,91,96,105,106,109,110,111,112,116,119,128,129,130,144]; % input E2-SQARE-Ai12 cell ID
        Cell_type_C_ID = [8,13,18,23,27,28,37,41,48,55,59,70,71,73,75,97]; % input E2-Shock-Ai12 cell ID


%% calculate occurence between SINGLE factors

[sinA_Res_occ, sinA_Res_hist_all, sinA_Res_hist_select, sinA_Res_num] = ...
     fxn_bin_ocurrence_hist_single(ca_temp, cut_frame_num, bin_frame_num, ...
     Cell_type_A_ID, cal_range_start, cal_range_end);    
 
 [sinB_Res_occ, sinB_Res_hist_all, sinB_Res_hist_select, sinB_Res_num] = ...
     fxn_bin_ocurrence_hist_single(ca_temp, cut_frame_num, bin_frame_num, ...
     Cell_type_B_ID, cal_range_start, cal_range_end);    
 
 [sinC_Res_occ, sinC_Res_hist_all, sinC_Res_hist_select, sinC_Res_num] = ...
     fxn_bin_ocurrence_hist_single(ca_temp, cut_frame_num, bin_frame_num, ...
     Cell_type_C_ID, cal_range_start, cal_range_end);    
    
%% calculate occurence between TWO factors

[dob_AB_Res_occ, dob_AB_Res_hist_all, dob_AB_Res_hist_select, dob_AB_Res_A_num, dob_AB_Res_B_num] = ...
     fxn_bin_ocurrence_hist_double(ca_temp, cut_frame_num, bin_frame_num, ...
     Cell_type_A_ID, Cell_type_B_ID,  cal_range_start, cal_range_end);    
 
[dob_AC_Res_occ, dob_AC_Res_hist_all, dob_AC_Res_hist_select, dob_AC_Res_A_num, dob_AC_Res_C_num] = ...
     fxn_bin_ocurrence_hist_double(ca_temp, cut_frame_num, bin_frame_num, ...
     Cell_type_A_ID, Cell_type_C_ID,  cal_range_start, cal_range_end);    
 
%%  calculate occurence between THREE factors

[tri_ABC_Res_occ, tri_ABC_Res_hist_all, tri_ABC_Res_hist_select, tri_ABC_Res_A_num, tri_ABC_Res_B_num, tri_ABC_Res_C_num]  = ...
     fxn_bin_ocurrence_hist_triple(ca_temp, cut_frame_num, bin_frame_num, ...
     Cell_type_A_ID, Cell_type_B_ID, Cell_type_C_ID, cal_range_start, cal_range_end);

 %%