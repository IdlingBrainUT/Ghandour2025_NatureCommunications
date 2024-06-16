function [ca_sort_data] = func_sort(speed_digital, ca_filt_data, sort_frame)

%% sorting data selection
% r = xcorr(x,y)

% sort_frame = [55201:60000];
% sort_frame = [57601:60000];
cor_B_threshold = 0.05;

    cor_x = speed_digital(sort_frame,:); % 1CSˆÈ~‚Åfreezing‚Æcorrelation‚ð•]‰¿
    cor_y = ca_filt_data(sort_frame,:);
    cor_x_y = [cor_x, cor_y];
        [cor_matrix, cor_pval] = corrcoef(cor_x_y);

            cor_table(1,:) = cor_matrix(1,2:end);
            [cor_B, cor_I] = sort(cor_table(1,:),'descend');
            cor_table(2,:) = cor_B;
            cor_table(3,:) = cor_I;
            cor_table(4,:) = cor_B > cor_B_threshold;
            
                cor_table(5,:) = cor_pval(1,2:end);
                cor_table(6,:) = cor_pval(1,2:end) < 0.05;
                
            ca_sort_data = ca_filt_data(:,cor_I);
%%
end