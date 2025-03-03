%% 
function [ca_bin_time, ca_cell_num, ca_bin_raster] = fxn_mod_round_binning_time(ca_filt_data, bin_frame_num)
%% debug
fps_sampling = 20;
%%
    mod_val = mod(size(ca_filt_data,1), bin_frame_num);
    
if mod_val == 0
    display('   mod_val = 0, normal binning done!')
    for i = 1:size(ca_filt_data,2)
    ca_mod(:,i)  = ca_filt_data(:,i);
    end
    [ca_mod_bin,  ca_mod_mean] = funcHF_temporal_binning(ca_mod, bin_frame_num);
else    
    display('      mod_val is not = 0, mod_rounded binning done! Last frame was rounded')
    for i = 1:size(ca_filt_data,2)
    ca_mod(:,i)  = [ca_filt_data(:,i); zeros(bin_frame_num-mod_val,1)];
    end 
    [ca_mod_bin,  ca_mod_mean] = funcHF_temporal_binning(ca_mod(1:end-bin_frame_num,:), bin_frame_num);
end
% [ca_mod_bin,  ca_mod_mean] = funcHF_temporal_binning(ca_mod, bin_frame_num);
    
    % for logical
%      ca_logical = ca_mod_bin > 0; 
%      ca_digital = double(ca_logical);  

% for num
ca_bin_raster = ca_mod_bin; 

%      display('finish fxn_mod_binning')
%% add time information
ca_bin_time     = [1:size(ca_bin_raster,1)]/(fps_sampling/bin_frame_num);
ca_cell_num = [1:size(ca_bin_raster,2)];
%%

end