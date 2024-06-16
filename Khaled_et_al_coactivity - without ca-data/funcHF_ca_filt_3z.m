function ca_filt_data = funcHF_ca_filt_3z(ca_raw_data);

%% filtering code

    highpass = 0.01; % 30 sec highpass filter
    Ca_fs    = 20; % 20hz
    
    ca_temp_data = ca_raw_data   ;      %　vlockedデータの代入 ca 縦：時間　ｘ　横：細胞
   
        [b,a] = butter(1, highpass /(Ca_fs/2),'high'); %
        %[b,a] = cheby2(1, 10, highpass /(Ca_fs/2), 'high');% highpass =0.01 good  
        dataOut = filtfilt(b,a,ca_temp_data); %ゼロ位相filtfiltで ハイパス　フィルターを通す
        
%% zscoreも試す
% 	dataIn1 = zscore(dataOut(1:7200,:)); 	    %　縦にフィルターがかかる
%     dataIn2 = zscore(dataOut(7201:19200,:)); 	%　縦にフィルターがかかる
%     dataIn3 = zscore(dataOut(19201:31200,:)); 	%　縦にフィルターがかかる
%     dataIn4 = zscore(dataOut(31201:43200,:)); 	%　縦にフィルターがかかる
%     dataIn5 = zscore(dataOut(43201:50400,:)); 	%　縦にフィルターがかかる
%     dataIn6 = zscore(dataOut(50401:57600,:)); 	%　縦にフィルターがかかる 

%         dataIn = [dataIn1; dataIn2; dataIn3; dataIn4; dataIn5; dataIn6];

           dataIn =  zscore(dataOut(:,:));
        
%%
            ca_filt_data = (dataIn);

                negative = find(ca_filt_data<2); % 2SD以下をカット
    ca_filt_data(negative) = zeros(size(negative));
    
%         % 3SD後に -3 でベースラインに戻すコード
%     ca_filt_data_zero = ca_filt_data-3;                
%     negative = find(ca_filt_data_zero<3);
%     ca_filt_data_zero(negative) = zeros(size(negative));
   
%%
    
end