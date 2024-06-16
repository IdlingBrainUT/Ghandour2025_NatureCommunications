function funcHF_figure_stamp_3z_raster(ca_data,CS_stamp,US_stamp, figure_range,xticks_range, clim_range)
%% time info
CS_time = [1:57600]/20;
neuron_num = size(ca_data,2);
%% freezing and test Ca
       
%     figure_range = [60 600];
%     pca_a = 1; pca_b = 2;
%     clim_range  = [0 3];
%     xticks_range = [150:15:600];
    
%     figure; 
%     subplot(10,1,1); imagesc(CS_time, 1,CS_digital'); xlim(figure_range); title('Freezing');
%     xticks([60:15:600]); grid on; title('\fontsize{12}CS presentation'); 
    
    figure;
    subplot(10,1,1); plot(CS_time, US_stamp,'r','LineWidth',2); 
    hold on
    subplot(10,1,1); plot(CS_time, CS_stamp,'b','LineWidth',2); xlim(figure_range); ylim([0.1 2])
    title('\fontsize{12}CS and US time stamp'); 
    xticks(xticks_range)
    mymap = [1, 1, 1 ;0, 0, 0]; colormap(mymap); grid on; yticks([]); ;grid on
    
%         subplot(10,1,2); plot(CS_time, pca_score_ca(:,pca_a)); xlim(figure_range); %ylim([-1 15]); title(['PCA#' num2str(pca_i)]);
%         hold on
%         subplot(10,1,2); plot(CS_time, pca_score_ca(:,pca_b)); xlim(figure_range); %ylim([-1 15]); title(['PCA#' num2str(pca_i)]);
%         legend('PC-1','PC-2');xticks([60:15:600]); grid on; title('\fontsize{12}CS presentation'); 
        
     subplot(10,1,2:10); imagesc(CS_time, 1:size(ca_data,2), ca_data'); clim([clim_range]); xlim(figure_range);
         title('\fontsize{12}Neuronal activity'); 
         xlabel('Time (s)','FontSize',12,'FontWeight','bold','Color','k'); 
         ylabel('Neuron#','FontSize',12,'FontWeight','bold','Color','k'); 
         mymap = [1, 1, 1 ;0, 0, 0]; colormap(mymap); grid on; xticks([xticks_range]); 