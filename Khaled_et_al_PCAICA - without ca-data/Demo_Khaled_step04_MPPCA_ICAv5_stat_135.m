%% MPPCA-ICA calculation code
clc; clear; close all; disp('Previous data is cleared')
data_cell = {}; shift=1; % Don't modify this line
addpath('homemade_functions\') 
%% load data set
load('calcium_135.mat') % Ca data

%% load ca data 
ca_raw_data      = calcium; % input time x neuron matrix
bin_frame_num    = 4;   % 20   -> 1s, 1s binning, For PCA-ICA, I reccomend to not change here, first.

%% %%  mouse sod 5 session frame information

% input frame info
data_cell{1+shift,1} = ('Day1pre');   data_cell{1+shift,2} = [1241	:8049] ;  
data_cell{2+shift,1} = ('Day2pre');    data_cell{2+shift,2} = [9291	:10639], [12681 :14029] ;  
data_cell{3+shift,1} = ('A');   data_cell{3+shift,2} = [15614	:22986] ;
data_cell{4+shift,1} = ('Day2post');   data_cell{4+shift,2} = [24491	:25839], [29030 :30378] ;
data_cell{5+shift,1} = ('Day3pre');   data_cell{5+shift,2} = [32779	:41565] ;
data_cell{6+shift,1} = ('A2');    data_cell{6+shift,2} = [41566	:42914] ;
data_cell{7+shift,1} = ('B');    data_cell{7+shift,2} = [45330	:46678] ;
data_cell{8+shift,1} = ('Day3post');    data_cell{8+shift,2} = [50225	:56683] ; 

%% data processing
% This code includes z-score processing, even you select either "yes" or "no" to filter fluctuation.

% Plseae choose 'y' (yes.)
[result_ca_filt_data_ct, result_data_cell] = fxn_MPPCA_ca_data_process(data_cell, ca_raw_data, bin_frame_num);
%% # Assign reference and target section. Select either this section or below, (Disable eiteher here or below)
% # reference vs. all session analysis mode. I recomend this first to check data aspect.
reference_session_num = 7;

data_bin_z_reference = result_data_cell{reference_session_num+shift,7}; % input data for MPPCA-ICA
data_bin_z_target    = result_ca_filt_data_ct; % input data for MPPCA-ICA
%% # Assign reference and target section. Select either this section or below, (Disable eiteher here or above)
% # session by session analysis mode. In this mode, disable the color-highlight funciton.
% reference_session_num = 4;
% target_session_num = 2;
% 
% data_bin_z_reference = zscore(result_data_cell{reference_session_num+shift,7}); % input data for MPPCA-ICA
% data_bin_z_target = zscore(result_data_cell{target_session_num+shift,7}); % input data for MPPCA-ICA

%% Parameters
prms_forced_IC              = 0;        % def=0;    0:MP-PCA, 1-x:forced_IC num input;  
prms_ICA_mode               = 3;        % def=3 (MATLAB ICA code);    1:fastica(Disabled), 2:fastICA(Disabled), 3:Reconstruction ICA (RICA)
prms_SD_thr                 = 2.0;      % def=2.5   (SD) ensemble weight threshold
prms_reactivation_SD_thr    = 2.0;        % def=2;    (SD) filled coloring threshold for reactivation strength 
prms_target_region1         = [0: 0] ;   % yellow, input second scale. like as [a:b]. disable: [0:0] 
prms_target_region2         = [0: 0] ;   % green, input second scale. like as [c:d]. disable: [0:0]
prms_ticklabel_mode         = 1 ;       % def=0;    0:off mode. 1:session name on mode. when you select whole data as target reference
%% Calculate MPPCA-ICA
[result_MPPCA] = fxn_MPPCA_ICAv05(data_bin_z_reference, data_bin_z_target, prms_forced_IC, prms_ICA_mode, prms_SD_thr,...
                                   prms_reactivation_SD_thr, prms_target_region1, prms_target_region2, bin_frame_num, ...
                                   result_data_cell, prms_ticklabel_mode); % calculation
%% Data handling
disp(result_MPPCA);
disp('Calculaton finished!');

% See the struct-variable result_PCAICA.
% In this struct, there are the data for quantification and cell-ID for
% tracking.

% ith_decompose         is significant patterns embedded within calcium data
% r_strength_ref        is r_strength(z, zscored) values in reference session.
% r_strength_target     is r_strength(z, zscored) values in target session.
% neuorn_weight         is neuron_weight for each pattern. This vector is visualized by stem graph.
% neuron_sig_IDs        is neuron ID contributing each pattern by thresholding.
%% Statistical analysis
thr_freq_search = 2; % def=2 SD. Recommended to be same to prms_reactivation_SD_thr. Threshold for frequency cal. McHugh uses 5.
stat_ref = 4; % reference section
stat_tar = 2; % target section
thr_stat_fold = 0; % def=0, I reccomend zero.

result_MPPCA_stat = fxn_MPPCA_ICA_stat(result_MPPCA, thr_freq_search, stat_ref, stat_tar, thr_stat_fold);

%%