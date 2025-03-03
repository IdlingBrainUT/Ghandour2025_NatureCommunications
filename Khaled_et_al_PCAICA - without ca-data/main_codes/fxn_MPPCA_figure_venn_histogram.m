function fxn_MPPCA_figure_venn_histogram_as_PDF(result, prms, extension)
%% decompose prms
ensemble_num1        = prms.ensemble_num1     ;
ensemble_num2        = prms.ensemble_num2     ;
thr_percentile      = prms.thr_percentile     ;
result_connect_cell = result.result_connect_cell ;

index_shuffle = result.index_shuffle;
index_connect = result.index_connect;
vennIDs_spec_A = result.vennIDs_spec_A;
vennIDs_spec_B = result.vennIDs_spec_B;
vennIDs_both_AB = result.vennIDs_both_AB;
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
ax1 = figure('Position', [100 100 ensemble_num2*120 ensemble_num1*100]);
% ax1 = gca;

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

file_name1 = ['Figure_Venn.', extension];
% ax1 = gca;
exportgraphics(ax1,file_name1, 'ContentType','vector');

%%
ax2 = figure('Position', [100 100 ensemble_num2*150 ensemble_num1*130]);
% ax2 = gca;

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
            xline((index_connect{i,ii})    ,'r','Real','LineWidth',2)
            xline((index_thr_bottom{i, ii}),'b--','LineWidth',1)
            xline((index_thr_top{i, ii})   ,'b--','LineWidth',1)
        end
    title(['Pattern #' ,num2str(i), ' vs. #', num2str(ii)])
    end
end

file_name2 = ['Figure_Histogram.', extension];
% ax2 = gca;
exportgraphics(ax2, file_name2, 'ContentType','vector');

%%

%%

end