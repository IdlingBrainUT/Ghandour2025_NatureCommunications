%% input data

    ca_raw_data = finalv_ai12;

    %% funcHF_ca_filt_3z
    ca_filt_data = funcHF_ca_filt_3z(ca_raw_data);
    
%% PCA % [coeff_ca,score_ca,latent_ca] = pca(Ca_temp_data);

     [pca_coeff_ca, pca_score_ca, pca_latent_ca] = pca(ca_filt_data); % normal 
%% time info
ca_time = [1:size(ca_filt_data,1)]/20;
neuron_num = size(ca_filt_data,2);
%% freezing and test Ca
       figure
    figure_range = [0 3535];
    pca_a = 1; pca_b = 2;  pca_c = 3;  pca_d = 4;
    clim_range  = [0 2];
 %   xticks_range = [150:15:600];
    
        subplot(10,1,1); plot(ca_time, pca_score_ca(:,pca_a)); xlim(figure_range); ylim([-1 15]); %title(['PCA#' num2str(pca_i)]);
        hold on
        subplot(10,1,1); plot(ca_time, pca_score_ca(:,pca_b)); xlim(figure_range); ylim([-1 15]); %title(['PCA#' num2str(pca_i)]);
        hold on
        subplot(10,1,1); plot(ca_time, pca_score_ca(:,pca_c)); xlim(figure_range); ylim([-1 15]); %title(['PCA#' num2str(pca_i)]);
        hold on
        subplot(10,1,1); plot(ca_time, pca_score_ca(:,pca_d)); xlim(figure_range); ylim([-1 15]); %title(['PCA#' num2str(pca_i)]);
        legend('PC-a','PC-b','PC-c','PC-d'); grid on; title('\fontsize{12}PCA results'); 
        
     subplot(10,1,3:10); imagesc(ca_time, 1:size(ca_filt_data,2), ca_filt_data');  xlim(figure_range); clim(clim_range);
     title('\fontsize{12}Neuronal activity'); 

         xlabel('Time (s)','FontSize',12,'FontWeight','bold','Color','k'); 
         ylabel('Neuron#','FontSize',12,'FontWeight','bold','Color','k'); 
%          
         mymap = [1, 1, 1 ;0, 0, 0]; colormap(mymap); grid on; %xticks([xticks_range]); 
     
 %% PCA heatmap
  figure
    figure_range = [0 3535];
    clim_range  = [0 2];
    
    figure;imagesc(ca_time, 1:size(ca_filt_data,2),pca_score_ca');xlim(figure_range); clim(clim_range);
    
   %%  
     