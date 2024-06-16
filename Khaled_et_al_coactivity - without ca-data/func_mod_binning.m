%%

function [ca_digital] = func_mod_binning(ca_filt_data, bin_frame_num);

    mod_val = mod(size(ca_filt_data,1), bin_frame_num);
    for i = 1:size(ca_filt_data,2)
    ca_mod(:,i)  = [ca_filt_data(:,i); zeros(bin_frame_num-mod_val,1)];
    end 
    
    [ca_mod_bin,  ca_mod_mean] = funcHF_temporal_binning(ca_mod, bin_frame_num);
    
     ca_logical = ca_mod_bin > 0; 
     ca_digital = double(ca_logical);  
end

%%