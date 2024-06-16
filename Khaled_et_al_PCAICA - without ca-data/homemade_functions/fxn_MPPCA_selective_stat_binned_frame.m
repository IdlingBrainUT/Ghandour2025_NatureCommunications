%% Statistical comparison among MPPCA ensembles across behavioral sesisons
function [result_MPPCA_stat] = fxn_MPPCA_selective_stat_binned_frame(select_cell, result_MPPCA, prms_MPPCA)
%% update note
% 211104: 1st code
% 220804: updated by USM request
%% debug mode
% thr_freq_search = 2; % def=2 SD, SD threshold for frequency cal. McHugh uses Raw RS 5.
% stat_ref = 16; % 
% stat_tar = 15; % 
% thr_stat_fold = 4; % 

%% Decompose parameters
global shift
prms_thr_freq_search   = prms_MPPCA.prms_thr_freq_search ; % def=2 SD, SD threshold for frequency cal. McHugh uses Raw RS 5.
prms_stat_ref          = prms_MPPCA.prms_stat_ref ; % 
prms_stat_tar          = prms_MPPCA.prms_stat_tar ; % 
prms_thr_stat_fold     = prms_MPPCA.prms_thr_stat_fold ; % 
prms_frame_num         = 1 ; % already binned
prms_sample_fps        = 1 ; % already binned
%% data processing

%new part
select_cell{1,1} = ('session'); select_cell{1,2} = ('frame range'); 
select_cell{1,3} = ('butter res'); select_cell{1,4} = ('Ca_bin_time');
select_cell{1,5} = ('Ca_bin_time'); select_cell{1,6} = ('Ca_cell_num'); select_cell{1,7} = ('Ca_select_binned');
select_cell{1,8} = ('Cumulative frames'); select_cell{1,9} = ('Cumulative time'); 
select_cell{1,10} = ('frame-t onset'); select_cell{1,11} = ('frame-t offset'); select_cell{1,12} = ('frame-t range'); 

total_session_num = size(select_cell(2:end,1),1);
ca_filt_data_ct = [];
cumulative_frames = 0;

ca_filt_data = result_MPPCA.data_bin_all_ct;
for i = 1:total_session_num
    select_cell{i+shift,3} = ca_filt_data(select_cell{i+shift,2},:); % butter filt

        select_cell{i+shift,4} = select_cell{i+shift,3}; % disable z-score cutoff 
    
        [select_cell{i+shift,5}, select_cell{i+shift,6}, select_cell{i+shift,7}] = ...
        fxn_mod_round_binning_time_v2(select_cell{i+shift,4}, prms_frame_num, prms_sample_fps);
    cumulative_frames = cumulative_frames + length(select_cell{i+shift,5});
    select_cell{i+shift,8} =  cumulative_frames;
    select_cell{i+shift,9} =  select_cell{i+shift,8}/(prms_sample_fps/prms_frame_num);    
    ca_filt_select_ct = cat(1, ca_filt_data_ct, select_cell{i+shift,7});    
end 

for i = 1:total_session_num
    if i == 1
    select_cell{i+shift,10} = 1;                 select_cell{i+shift,11} = select_cell{i+shift,8};
    select_cell{i+shift,12} = [select_cell{i+shift,10}: select_cell{i+shift,11}];
    else 
    select_cell{i+shift,10} = select_cell{i,8}+1;  select_cell{i+shift,11} = select_cell{i+shift,8};
    select_cell{i+shift,12} = [select_cell{i+shift,10}: select_cell{i+shift,11}];
    end
end

select_cell{1,13} = ('r_kt_raw'); select_cell{1,14} = ('r_kt_z'); 

for i = 1:total_session_num
select_cell{i+shift,13} =  result_MPPCA.r_strength_target(select_cell{i+shift,2},:);
select_cell{i+shift,14} =  result_MPPCA.r_strength_targetz(select_cell{i+shift,2},:);
end
%%
result_select_cell_stat = select_cell(:,[1,7,12,13,14]);

%%
num_of_session = size(result_select_cell_stat,1)-1;
num_of_ensemble = numel(result_MPPCA.neuron_sig_IDs);

result_select_cell_stat{1,6} = ('Mean'); result_select_cell_stat{1,7} = ('Frq search');
result_select_cell_stat{1,8} = ('Cut off RS'); result_select_cell_stat{1,9} = ('Cut off zRS');
result_select_cell_stat{1,9} = ('Cut off Mean'); result_select_cell_stat{1,10} = ('Mean fold change');

result_select_cell_stat{1,11} = ('Max'); result_select_cell_stat{1,12} = ('zero pad'); 
result_select_cell_stat{1,13} = ('Local maxima'); result_select_cell_stat{1,14} = ('Pad off local maxima');
result_select_cell_stat{1,15} = ('Maxima address'); result_select_cell_stat{1,16} = ('Maxima val');
result_select_cell_stat{1,17} = ('Num of maxima'); result_select_cell_stat{1,18} = ('Nomarized maxima');
result_select_cell_stat{1,19} = ('Binalized above thr'); result_select_cell_stat{1,20} = ('Num of raster');
result_select_cell_stat{1,21} = ('Frequency'); result_select_cell_stat{1,21} = ('Freq mean fold');
%% Mean cal
% 4 = raw_r_kt, 5 = zscored r_kt.

for i = 1:num_of_session
    result_select_cell_stat{shift+i,6} = mean(result_select_cell_stat{shift+i,5},1); % mean
    
    result_select_cell_stat{shift+i,7} = find(result_select_cell_stat{shift+i,5} < prms_thr_freq_search );
    result_select_cell_stat{shift+i,8} = result_select_cell_stat{shift+i,5};
    result_select_cell_stat{shift+i,8}(result_select_cell_stat{shift+i,7}) = 0; % replace to zero

    result_select_cell_stat{shift+i,9} = mean(result_select_cell_stat{shift+i,8},1); % mean

    result_select_cell_stat{shift+i,11} = max(result_select_cell_stat{shift+i,5},[],1); % Max
    
    for ii = 1:num_of_ensemble
    result_select_cell_stat{shift+i,12}(:,ii) = [0; (result_select_cell_stat{shift+i,8}(:,ii));0] ;
    result_select_cell_stat{shift+i,13}(:,ii) = double(islocalmax(result_select_cell_stat{shift+i,12}(:,ii)));
    result_select_cell_stat{shift+i,14}(:,ii) = result_select_cell_stat{shift+i,13}(2:end-1,ii);
    result_select_cell_stat{shift+i,15}{1,ii} = find(result_select_cell_stat{shift+i,14}(:,ii) == 1);
    result_select_cell_stat{shift+i,16}{1,ii} = result_select_cell_stat{shift+i,8}( (result_select_cell_stat{shift+i,15}{1,ii}),ii);
    result_select_cell_stat{shift+i,17}(:,ii) = numel(result_select_cell_stat{shift+i,16}{1,ii}); 
    result_select_cell_stat{shift+i,18}(:,ii) = sum(result_select_cell_stat{shift+i,16}{1,ii})./ size(result_select_cell_stat{shift+i,8},1);
    
    % binarize for frequency cal
    result_select_cell_stat{shift+i,19} = (result_select_cell_stat{shift+i,8} > 0);
    result_select_cell_stat{shift+i,20} = sum(result_select_cell_stat{shift+i,19});
    result_select_cell_stat{shift+i,21} = (result_select_cell_stat{shift+i,20}) ./ size(result_select_cell_stat{shift+i,8},1);
    end
end
%% stat for reactivation strength
for i = 1:num_of_session   
    result_select_cell_stat{shift+i,10} = result_select_cell_stat{shift+i,9} ./ result_select_cell_stat{shift+prms_stat_ref,9}; % mean   
end

temp_fold_change_stat_RS = {};

for i = 1:num_of_ensemble
temp_fold_change_stat_RS{1,i} = result_select_cell_stat{shift+prms_stat_ref,9}(i) ./ result_select_cell_stat{shift+prms_stat_tar,9}(i);
[temp_fold_change_stat_RS{2,i}, temp_fold_change_stat_RS{3,i}, temp_fold_change_stat_RS{4,i}] = ...
    ranksum(result_select_cell_stat{shift+prms_stat_ref,8}(:,i), result_select_cell_stat{shift+prms_stat_tar,8}(:,i));  
temp_fold_change_stat_RS{5,i} =  temp_fold_change_stat_RS{1,i} >= prms_thr_stat_fold && temp_fold_change_stat_RS{3,i} == 1;
end
%% stat for frequency
for i = 1:num_of_session   
    result_select_cell_stat{shift+i,22} = result_select_cell_stat{shift+i,21} ./ result_select_cell_stat{shift+prms_stat_ref,21}; % mean   
end

temp_fold_change_stat_Frq = {};

for i = 1:num_of_ensemble
temp_fold_change_stat_Frq{1,i} = result_select_cell_stat{shift+prms_stat_ref,21}(i) ./ result_select_cell_stat{shift+prms_stat_tar,21}(i);
[temp_fold_change_stat_Frq{2,i}, temp_fold_change_stat_Frq{3,i}, temp_fold_change_stat_Frq{4,i}] = ...
    ranksum(result_select_cell_stat{shift+prms_stat_ref,19}(:,i), result_select_cell_stat{shift+prms_stat_tar,19}(:,i));  
temp_fold_change_stat_Frq{5,i} =  temp_fold_change_stat_Frq{1,i} >= prms_thr_stat_fold && temp_fold_change_stat_Frq{3,i} == 1;
end
%% Sorting by stat results
posi_ID = find(cell2mat(temp_fold_change_stat_RS(5,:)) == 1);
nega_ID = 1:num_of_ensemble;
nega_ID(posi_ID) = [];
%% result stat mean
result_stat_mean = cell2mat(result_select_cell_stat(2:end,9)); % 9: raw RS, 10: fold change of RS to target
result_stat_mean_posi = mean(result_stat_mean(:,posi_ID),2);
result_stat_mean_nega = mean(result_stat_mean(:,nega_ID),2);

figure('Position',[400 100 1400 800]); % imagesc(result_stat_mean);
subplot(221)
plot(result_stat_mean(:,posi_ID),'m','LineWidth',1)
hold on
plot(result_stat_mean(:,nega_ID),'c','LineWidth',1)
hold on
plot(result_stat_mean_posi,'r','LineWidth',2)
hold on
plot(result_stat_mean_nega,'b','LineWidth',2)

xticklabels(result_select_cell_stat(1+shift:size(result_select_cell_stat,1),1) );  % for adding ticklabel
xticks([1:size(result_select_cell_stat,1)-1]); % for adding ticklabel
xtickangle(45) % update 211104
ylabel('Reactivation strength (SD)');
% ylim([0 12])

set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 9, 'FontName','Arial');
hold off; box off
%% result stat maxima
result_stat_maxima = cell2mat(result_select_cell_stat(2:end,18));
result_stat_maxima_posi = mean(result_stat_maxima(:,posi_ID),2);
result_stat_maxima_nega = mean(result_stat_maxima(:,nega_ID),2);

% figure('Position',[400 100 1000 600]); % imagesc(result_stat_maxima);
subplot(222)
plot(result_stat_maxima(:,posi_ID),'m','LineWidth',1)
hold on
plot(result_stat_maxima(:,nega_ID),'c','LineWidth',1)
hold on
plot(result_stat_maxima_posi,'r','LineWidth',2)
hold on
plot(result_stat_maxima_nega,'b','LineWidth',2)

xticklabels(result_select_cell_stat(1+shift:size(result_select_cell_stat,1),1) );  % for adding ticklabel
xticks([1:size(result_select_cell_stat,1)-1]); % for adding ticklabel
xtickangle(45) % update 211104
ylabel('Mean of reactivation maxima (SD)');
% ylim([0 0.5])

set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 8, 'FontName','Arial');
hold off; box off
%% Max
result_stat_max_only = cell2mat(result_select_cell_stat(2:end,11));
result_stat_max_only_posi = mean(result_stat_max_only(:,posi_ID),2);
result_stat_max_only_nega = mean(result_stat_max_only(:,nega_ID),2);

% figure('Position',[400 100 1000 600]); % imagesc(result_stat_max_only);
subplot(223)
plot(result_stat_max_only(:,posi_ID),'m','LineWidth',1)
hold on
plot(result_stat_max_only(:,nega_ID),'c','LineWidth',1)
hold on
plot(result_stat_max_only_posi,'r','LineWidth',2)
hold on
plot(result_stat_max_only_nega,'b','LineWidth',2)

xticklabels(result_select_cell_stat(1+shift:size(result_select_cell_stat,1),1) );  % for adding ticklabel
xticks([1:size(result_select_cell_stat,1)-1]); % for adding ticklabel
xtickangle(45) % update 211104
ylabel('Reactivation max only (SD)');
% ylim([0 0.5])

set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 8, 'FontName','Arial');
hold off; box off
%% Sorting by stat results
posi_ID = find(cell2mat(temp_fold_change_stat_Frq(5,:)) == 1);
nega_ID = 1:num_of_ensemble;
nega_ID(posi_ID) = [];
%% Frequency above threshold
result_stat_frq = cell2mat(result_select_cell_stat(2:end,21));
result_stat_frq_posi = mean(result_stat_frq(:,posi_ID),2);
result_stat_frq_nega = mean(result_stat_frq(:,nega_ID),2);

% figure('Position',[400 100 1000 600]); % imagesc(result_stat_frq);
subplot(224)
plot(result_stat_frq(:,posi_ID),'m','LineWidth',1)
hold on
plot(result_stat_frq(:,nega_ID),'c','LineWidth',1)
hold on
plot(result_stat_frq_posi,'r','LineWidth',2)
hold on
plot(result_stat_frq_nega,'b','LineWidth',2)

xticklabels(result_select_cell_stat(1+shift:size(result_select_cell_stat,1),1) );  % for adding ticklabel
xticks([1:size(result_select_cell_stat,1)-1]); % for adding ticklabel
xtickangle(45) % update 211104
ylabel('Reactivation frequency (Hz)');
% ylim([0 0.5])

set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 8, 'FontName','Arial');
hold off; box off
%% result part
result_MPPCA_stat.fold_change_stat_table       = temp_fold_change_stat_RS;
result_MPPCA_stat.result_stat_mean_nega        = result_stat_mean_nega;
result_MPPCA_stat.result_stat_mean_posi        = result_stat_mean_posi;
result_MPPCA_stat.result_stat_mean_nega_indivi = result_stat_mean(:,nega_ID);
result_MPPCA_stat.result_stat_mean_posi_indivi = result_stat_mean(:,posi_ID);
result_MPPCA_stat.result_stat_maxima           = result_stat_maxima;
result_MPPCA_stat.result_stat_max_only         = result_stat_max_only;
result_MPPCA_stat.result_stat_mean             = result_stat_mean; % reactiovation strength, used generally.
result_MPPCA_stat.result_stat_frq              = result_stat_frq;  % reactiovation frequency, used generally.
%%
end