function [ca_filt_data] = fxn_butter_filt(ca_raw_data, sample_fps, highpass, factor)
%% for debug,
% ca_raw_data  = ca_data; % input ca data
% sample_fps   = 20; % 20hz
% highpass     = 0.01; % 30 sec highpass filter
% factor       = 1; % factor
%% filtering code
[b,a] = butter(factor, highpass /(sample_fps/2),'high'); %
% [b,a] = cheby2(factor, 10, highpass /(Ca_fs/2), 'high');% highpass =0.01 good  
ca_filt_data = filtfilt(b, a, ca_raw_data); %ゼロ位相filtfiltで ハイパス　フィルターを通す  
%% figure
figure;
subplot(311); imagesc(ca_raw_data'); caxis([0 3]); title('Before')
subplot(312); imagesc(ca_filt_data'); caxis([0 3]); title('After')
subplot(313); plot(mean(ca_raw_data', 1),'b'); hold on;
              plot(mean(ca_filt_data', 1), 'r'); xlim([0 size(ca_filt_data,1)]); legend('Before','After');
%%    
end