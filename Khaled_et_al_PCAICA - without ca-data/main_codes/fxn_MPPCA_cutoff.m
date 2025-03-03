function [ca_cutoff_data] = fxn_MPPCA_cutoff(ca_raw_data, cutoff_val)
%% for debug,
% cutoff_val   = 0; % 20hz
%% cutoff section
ca_temp = ca_raw_data;
negative = find(ca_temp < cutoff_val); % 2SDˆÈ‰º‚ðƒJƒbƒg
ca_temp(negative) = zeros(size(negative));
ca_cutoff_data = ca_temp;
display('Finish cutoff!')
%% figure for debug
% figure;
% subplot(311); imagesc(ca_raw_data'); caxis([0 3]); title('Before')
% subplot(312); imagesc(ca_cutoff_data'); caxis([0 3]); title('After')
% subplot(313); plot(mean(ca_raw_data', 1),'b'); hold on;
%               plot(mean(ca_cutoff_data', 1), 'r'); xlim([0 size(ca_cutoff_data,1)]); legend('Before','After');
%%    
end