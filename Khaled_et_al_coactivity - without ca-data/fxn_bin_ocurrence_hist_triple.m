%% mostafa binning code

function [result_sum_dur_range, data_for_hist1, data_for_hist2, i, ii, iii] = ...
    func_bin_ocurrence_hist_triple(ca_filt_data, cut_frame_num, bin_frame_num, Cell_type_A_ID, Cell_type_B_ID, Cell_type_C_ID, cal_range_start, cal_range_end)

% ca_filt_data -> input data
% cut_frame_num = 1; % 0 = no cut, 1 = 50ms-frame-firing cut
% bin_frame_num = 10; % 500ms binning
% cal_range_start = 1000; % input range in second scale
% cal_range_end   = 1300; % input range in second scale
% 
% Cell_type_A_ID = [1,2,3,4,5,6,7];
% Cell_type_B_ID = [10,11,12,13,14,15,16,17,18];
% Cell_type_C_ID = [20,21,22,23,24,25,26,27,28,29];


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
 
  result_ca_1s_bin_occurence(:,1) = sum(ca_mod_bin,2);
  result_ca_1s_bin_occurence(:,2) = sum(ca_mod_bin,2)/size(ca_mod_bin,2)*100;
  
%% occurence calculation

    % スタートとエンドを別々に入力する
%     cal_range_start = 100; % input range in second scale
%     cal_range_end   = 200; % input range in second scale
    
    xlim_range = [cal_range_start, cal_range_end];
    cal_range_for_bin_data = [cal_range_start*20./bin_frame_num : cal_range_end*20./bin_frame_num ];
    
    result_cal_data_cut = result_ca_1s_bin_occurence(cal_range_for_bin_data,:);    
    
    [hist_B1, ~] = find(result_cal_data_cut(:,1)>0);
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
         ylabel('occurence','FontSize',12,'FontWeight','bold','Color','k'); 
    
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

%% 

Cells_A = ca_mod_bin(:,Cell_type_A_ID);
Cells_B = ca_mod_bin(:,Cell_type_B_ID);
Cells_C = ca_mod_bin(:,Cell_type_C_ID);

%%

for i = 1:size(Cells_A,2)
    for ii=1:size(Cells_B,2) 
       for iii=1:size(Cells_C,2)
            for iiii = 1:size(Cells_A,1)

                    if Cells_A(iiii,i)==1 && Cells_B(iiii,ii)==1 && Cells_C(iiii,iii)==1
                        Cells_A_B_C(i,ii,iii,iiii) =1; 
                    else
                        Cells_A_B_C(i,ii,iii,iiii) =0;
                    end
            end
        end
    end
end

%%
% for i = 1:size(Cells_A,2)
% for ii=1:size(Cells_B,2) 
%    
%     for iii = 1:size(Cells_A,1)
%             if Cells_A(iii,i)==1 && Cells_B(iii,ii)==1
% %                Cells_A_B{i,ii,iii} =1;
%                Cells_A_B(i,ii,iii) =1; 
%             else
% %                Cells_A_B{i,ii,iii} =0;
%                 Cells_A_B(i,ii,iii) =0;
%             end
%     end
% end
% end


%% cal

temp_A = sum(Cells_A_B_C,3); 
temp_A_sqz = squeeze(temp_A);
temp_A2 = sum(temp_A_sqz,2); 
temp_A_sqz2 = squeeze(temp_A2);
temp_A_sqz_logical = temp_A_sqz2'> 0;

temp_B = sum(Cells_A_B_C,3); 
temp_B_sqz = squeeze(temp_B);
temp_B2 = sum(temp_B_sqz,1); 
temp_B_sqz2 = squeeze(temp_B2);
temp_B_sqz_logical = temp_B_sqz2'> 0;

temp_C = sum(Cells_A_B_C,1); 
temp_C_sqz = squeeze(temp_C);
temp_C2 = sum(temp_C_sqz,1); 
temp_C_sqz2 = squeeze(temp_C2);
temp_C_sqz_logical = temp_C_sqz2'> 0;

%% Calculate entire occurence

Cells_A_B_C_res_occurence_raster = [temp_A_sqz_logical, temp_B_sqz_logical, temp_C_sqz_logical];
Cells_A_B_C_res_occurence_sum = sum(Cells_A_B_C_res_occurence_raster,2);
Cells_A_B_C_res_occurence_vec = Cells_A_B_C_res_occurence_sum > 0;

%% figure

    title_0 = ('Sorted ');
    title_0123 = [title_0, title_1, num2str(cut_frame_num), title_2, num2str(xlim_range), title_3];

% figure; 
subplot(222)
imagesc( (1:size(Cells_A_B_C_res_occurence_raster',2))*bin_frame_num./20 , (1:size(Cells_A_B_C_res_occurence_raster',1)) ,Cells_A_B_C_res_occurence_raster')
% imagesc(Cells_A_B_res_occurence_raster')

xlim(xlim_range); clim(clim_range);
         title(title_0123)

         xlabel('Time (s)','FontSize',12,'FontWeight','bold','Color','k'); 
         ylabel('Neuron#','FontSize',12,'FontWeight','bold','Color','k'); 
%          
         mymap = [1, 1, 1 ;0, 0, 0]; colormap(mymap); grid on; %xticks([xticks_range]); 


%% occurence calculation
% cal_range_in_sec
%     cal_range_in_sec = [100 : 1000]; % input range in second scale
%     xlim_range = [cal_range_in_sec(1) cal_range_in_sec(end)];

data_range_vec = Cells_A_B_C_res_occurence_vec(cal_range_for_bin_data,:);

Cells_A_B_res_occurence_vec_total_range = sum(data_range_vec);
result_sum_dur_range = num2str(Cells_A_B_res_occurence_vec_total_range);
disp_word_res = ['Occurence value is ', result_sum_dur_range];
disp(disp_word_res);

%% figure occurence hist

    data_range_sum = Cells_A_B_C_res_occurence_sum(cal_range_for_bin_data,:);
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
         ylabel('occurence','FontSize',12,'FontWeight','bold','Color','k'); 
%          
         mymap = [1, 1, 1 ;0, 0, 0]; colormap(mymap); grid on; %xticks([xticks_range]); 

%          disp('show histogram results in sorted data');

         
%%

% figure; 
% subplot(311); plot(temp_A_sqz_logical(:,:))
% subplot(312); plot(temp_B_sqz_logical(:,:))
% subplot(313); plot(temp_C_sqz_logical(:,:))

%%
end
