function [Ca_data_binned,  Ca_data_binned_mean] = funcHF_temporal_binning(Ca_data, binning_num)
%% Temporal binning (1 sec)

%         Ca_data     = HFcal01kikgr75ca1korepfinalv;
%         binning_num = 20;
    
    [n_of_frame, n_of_cell]  = size(Ca_data);
    Ca_data_binned = zeros(n_of_frame./binning_num, n_of_cell); % m x n initialize
    
    for i = 1:n_of_cell
     Ca_data_temp = reshape(Ca_data(:,i), binning_num, []);
     Ca_data_binned(:,i) = mean(Ca_data_temp,1);
    end
     
    Ca_data_binned_mean   = mean(Ca_data_binned,2);
    
%%
end