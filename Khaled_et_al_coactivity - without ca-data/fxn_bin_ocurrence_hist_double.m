%% mostafa binning code

function [result_sum_dur_range, data_for_hist1, data_for_hist2, i, ii] = ...
    func_bin_ocurrence_hist(ca_filt_data, cut_frame_num, bin_frame_num, Cell_type_A_ID, Cell_type_B_ID, cal_range_start, cal_range_end)

% ca_filt_data -> input data
% cut_frame_num = 1; % 0 = no cut, 1 = 50ms-frame-firing cut
% bin_frame_num = 10; % 500ms binning
% cal_range_start = 100; % input range in second scale
% cal_range_end   = 200; % input range in second scale
% Cell_type_A_ID = [1,2,3,4,5];
% Cell_type_B_ID = [6,7,8,9,10,11,12,13,14,15];


%% 
% ca_filt_less_fire = func_ca_filt_less_file(ca_filt_data, cut_frame_num) ; % cut_frame_num = 1, cut 50ms-firing frame 

ca_filt_temp = ca_filt_data;
% cut_frame_num = 1; % 0 = no cut, 1 = 50ms-frame-firing cut, 

        ca_filt_less_fire = func_ca_filt_less_file(ca_filt_temp, cut_frame_num);
        
        % checking 
        compare_firing(:,1) = ca_filt_temp(:,1);
        compare_firing(:,2) = ca_filt_less_fire(:,1);   

%% 500ms binning, 20fps, 1frame = 50ms, 500ms-binnning 10frame binning

% function [ca_digital] = func_mod_binning(ca_filt_data, bin_frame_num);

% bin_frame_num = 10; % 500ms binning
 [ca_mod_bin] = func_mod_binning(ca_filt_less_fire, bin_frame_num);
 
  result_ca_1s_bin_occurance(:,1) = sum(ca_mod_bin,2);
  result_ca_1s_bin_occurance(:,2) = sum(ca_mod_bin,2)/size(ca_mod_bin,2)*100;
  
%% Occurance calculation

    % スタートとエンドを別々に入力する
%     cal_range_start = 100; % input range in second scale
%     cal_range_end   = 200; % input range in second scale
    
    xlim_range = [cal_range_start, cal_range_end];
    cal_range_for_bin_data = [cal_range_start*20./bin_frame_num : cal_range_end*20./bin_frame_num ];
    
    result_cal_data_cut = result_ca_1s_bin_occurance(cal_range_for_bin_data,:);    
    
    [hist_B1, ~] = find(result_cal_data_cut(:,1)>0); % changed
    data_for_hist1 = result_cal_data_cut(hist_B1);
    
    figure; subplot(223);
    hist = histogram(data_for_hist1(:,1));
    
    hist.NumBins = 20;
    hist.BinEdges = [0:20];
%     ylim([0 30])
    
    title_1 = ('cut-frame-num=');
    title_2 = (', cal-range=');
    title_3 = (' in sec');
    title_123 = [title_1, num2str(cut_frame_num), title_2, num2str(xlim_range), title_3];
    
    title(title_123)
    
             xlabel('Neurons','FontSize',12,'FontWeight','bold','Color','k'); 
         ylabel('Occurance','FontSize',12,'FontWeight','bold','Color','k'); 
    
%% figure

clim_range = [0 1];

% figure;
subplot(221);

imagesc( (1:size(ca_mod_bin',2))*bin_frame_num./20 , (1:size(ca_mod_bin',1)) ,ca_mod_bin')

xlim(xlim_range); clim(clim_range);
         title(title_123)

         xlabel('Time (s)','FontSize',12,'FontWeight','bold','Color','k'); 
         ylabel('Neuron#','FontSize',12,'FontWeight','bold','Color','k'); 
%          
         mymap = [1, 1, 1 ;0, 0, 0]; colormap(mymap); grid on; %xticks([xticks_range]); 
%% Input cell ID

% Cell_type_A_ID = [1,2,3,4,5];
% Cell_type_B_ID = [6,7,8,9,10,11,12,13,14,15];
% Cell_type_C_ID = [31,32,33];

%% 

Cells_A = ca_mod_bin(:,Cell_type_A_ID);
Cells_B = ca_mod_bin(:,Cell_type_B_ID);
% Cells_C = ca_1s_bin(:,Cell_type_C_ID);

%%

% Cells_A_B = {};

%%
for i = 1:size(Cells_A,2)
for ii=1:size(Cells_B,2) 
   
    for iii = 1:size(Cells_A,1)
            if Cells_A(iii,i)==1 && Cells_B(iii,ii)==1
%                Cells_A_B{i,ii,iii} =1;
               Cells_A_B(i,ii,iii) =1; 
            else
%                Cells_A_B{i,ii,iii} =0;
                Cells_A_B(i,ii,iii) =0;
            end
    end
end
end


%% cal

Cells_A_B_res_A = sum(Cells_A_B,2); 
Cells_A_B_res_A_sqz = squeeze(Cells_A_B_res_A);
Cells_A_B_res_A_sqz_logical = Cells_A_B_res_A_sqz'> 0;
% Cells_A_B_res_A_sqz_logical_vec = sum(Cells_A_B_res_A_sqz_logical,2) > 0;
% Cells_A_B_res_occurance_vecと一緒だからいらない

Cells_A_B_res_B = sum(Cells_A_B,1); 
Cells_A_B_res_B_sqz = squeeze(Cells_A_B_res_B);
Cells_A_B_res_B_sqz_logical = Cells_A_B_res_B_sqz'> 0;
% Cells_A_B_res_B_sqz_logical_vec = sum(Cells_A_B_res_B_sqz_logical,2) > 0;
% Cells_A_B_res_occurance_vecと一緒だからいらない

%% Calculate entire Occurance

Cells_A_B_res_occurance_raster = [Cells_A_B_res_A_sqz_logical, Cells_A_B_res_B_sqz_logical];
Cells_A_B_res_occurance_sum = sum(Cells_A_B_res_occurance_raster,2);
Cells_A_B_res_occurance_vec = Cells_A_B_res_occurance_sum > 0;

%% figure

    title_0 = ('Sorted ');
    title_0123 = [title_0, title_1, num2str(cut_frame_num), title_2, num2str(xlim_range), title_3];

% figure; 
subplot(222)
imagesc( (1:size(Cells_A_B_res_occurance_raster',2))*bin_frame_num./20 , (1:size(Cells_A_B_res_occurance_raster',1)) ,Cells_A_B_res_occurance_raster')
% imagesc(Cells_A_B_res_occurance_raster')

xlim(xlim_range); clim(clim_range);
         title(title_0123)

         xlabel('Time (s)','FontSize',12,'FontWeight','bold','Color','k'); 
         ylabel('Neuron#','FontSize',12,'FontWeight','bold','Color','k'); 
%          
         mymap = [1, 1, 1 ;0, 0, 0]; colormap(mymap); grid on; %xticks([xticks_range]); 


%% Occurance calculation
% cal_range_in_sec
%     cal_range_in_sec = [100 : 1000]; % input range in second scale
%     xlim_range = [cal_range_in_sec(1) cal_range_in_sec(end)];

data_range_vec = Cells_A_B_res_occurance_vec(cal_range_for_bin_data,:);

Cells_A_B_res_occurance_vec_total_range = sum(data_range_vec);
result_sum_dur_range = num2str(Cells_A_B_res_occurance_vec_total_range);
disp_word_res = ['Occurence value is ', result_sum_dur_range];
disp(disp_word_res);

%% figure Occurance hist

    data_range_sum = Cells_A_B_res_occurance_sum(cal_range_for_bin_data,:);
    [hist_B2, ~] = find(data_range_sum>0);
    data_for_hist2 = data_range_sum(hist_B2);

%     figure; 
    subplot(224)
    hist_sort = histogram(data_for_hist2);%
    hist_sort.NumBins = 15;
    hist_sort.BinEdges = [0:20];
%     ylim([0 10])
    

    title(title_0123)
         xlabel('Neurons','FontSize',12,'FontWeight','bold','Color','k'); xticks([1:20]);
         ylabel('Occurance','FontSize',12,'FontWeight','bold','Color','k'); 
%          
         mymap = [1, 1, 1 ;0, 0, 0]; colormap(mymap); grid on; %xticks([xticks_range]); 

%          disp('show histogram results in sorted data');

         
%%
end
