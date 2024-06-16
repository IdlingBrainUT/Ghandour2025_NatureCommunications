function result = fxn_MPPCA_ensemble_connectivity_double(result_data_cell, result_MPPCA1, result_MPPCA2, prms)
%% comment
% 211118: 1st version
%% decompose prms
session_num_assess  = prms.session_num_assess ;
binarize_thr        = prms.binarize_thr       ;
shuffle_mode        = prms.shuffle_mode       ;
ensemble_num1        = prms.ensemble_num1     ;
ensemble_num2        = prms.ensemble_num2     ;
thr_percentile      = prms.thr_percentile     ;
%%
tic;
result_connect = {};
parfor i = 1:ensemble_num1
    for ii = 1:ensemble_num2
        [result_connect_cell{i,ii}] = fxn_MPPCA_ensemble_connectivity2_double ...
        (result_data_cell, result_MPPCA1, result_MPPCA2, session_num_assess, i, ii, ... 
        shuffle_mode, binarize_thr, thr_percentile);
    end
end
toc;
%% Results section
for i = 1:ensemble_num1
    for ii = 1:ensemble_num2
        index_shuffle{i, ii}    = result_connect_cell{i, ii}.shuffle_index;
        index_connect{i, ii}    = result_connect_cell{i, ii}.connect_index;
        vennIDs_spec_A{i, ii}   = result_connect_cell{i, ii}.vennIDs_spec_A;
        vennIDs_spec_B{i, ii}   = result_connect_cell{i, ii}.vennIDs_spec_B;
        vennIDs_both_AB{i, ii}  = result_connect_cell{i, ii}.vennIDs_both_AB;
        index_thr_bottom{i, ii} = result_connect_cell{i, ii}.thr_bottom;
        index_thr_top{i, ii}    = result_connect_cell{i, ii}.thr_top;
    end
end
%%
figure('Position', [100 100 700 500]);
for i = 1:ensemble_num1
    for ii = 1:ensemble_num2

        % subplot regulation
if ensemble_num1 == ensemble_num2        
subplot(max(ensemble_num1), max(ensemble_num2), (max(ensemble_num2).*(i-1) + ii))
elseif ensemble_num1 > ensemble_num2
subplot(max(ensemble_num1), max(ensemble_num2), (max(ensemble_num2).*(i-1) + ii))
elseif ensemble_nem1 < ensemble_num2
subplot(max(ensemble_num1), max(ensemble_num2), (max(ensemble_num1).*(i-1) + ii))
end

fxn_venn([numel(vennIDs_spec_A{i, ii}), numel(vennIDs_spec_B{i, ii}), numel(vennIDs_both_AB{i, ii})])
xticklabels([]); yticklabels([]); title(['Pattern #' ,num2str(i), ' vs. #', num2str(ii)])
    end
end
%%
figure('Position', [50 50 800 700]) ;
for i = 1:ensemble_num1
    for ii = 1:ensemble_num2
        
% subplot regulation
if ensemble_num1 == ensemble_num2        
subplot(max(ensemble_num1), max(ensemble_num2), (max(ensemble_num2).*(i-1) + ii))
elseif ensemble_num1 > ensemble_num2
subplot(max(ensemble_num1), max(ensemble_num2), (max(ensemble_num2).*(i-1) + ii))
elseif ensemble_nem1 < ensemble_num2
subplot(max(ensemble_num1), max(ensemble_num2), (max(ensemble_num1).*(i-1) + ii))
end
        
        histogram(cell2mat(index_shuffle{i, ii}))
        hold on
        if isnan(index_connect{i,ii})
%             disp('Skip xline drawing')
        else
            xline((index_connect{i,ii})    ,'r','Actual value','LineWidth',2)
            xline((index_thr_bottom{i, ii}),'b--','LineWidth',1)
            xline((index_thr_top{i, ii})   ,'b--','LineWidth',1)
        end
    title(['Pattern #' ,num2str(i), ' vs. #', num2str(ii)])
    end
end
%%
index_real   = cell2mat(index_connect);
index_bottom = cell2mat(index_thr_bottom);
index_top    = cell2mat(index_thr_top);

index_res_bottom = index_real <= index_bottom;
index_res_top    = index_real >= index_top;

index_res_matrix = (double(index_res_bottom).* -1) + double(index_res_top);
figure; imagesc(index_res_matrix); colormap(fxn_redblue);
title('Distribution thresholding (compared to 1000 times shuffling)'); xlabel('Pattern #'); ylabel('Pattern #');
grid on; clim([-1 1]);
%% Output results
result.result_connect_cell = result_connect_cell;
result.index_shuffle =  index_shuffle;
result.index_connect = index_connect;
result.vennIDs_spec_A = vennIDs_spec_A;
result.vennIDs_spec_B = vennIDs_spec_B;
result.vennIDs_both_AB = vennIDs_both_AB;
%%
end