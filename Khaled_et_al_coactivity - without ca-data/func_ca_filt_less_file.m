%%
function ca_filt_less_fire = func_ca_filt_less_file(ca_filt_data, cut_frame_num) ; % cut_frame_num = 1, cut 50ms-firing frame 

ca_filt_temp = ca_filt_data;
% interval_range = [1:size(ca_filt_temp,1)];
thr_val = 0;

        A1 = ca_filt_temp > thr_val;

        A1_vector = reshape(A1,1,[]);
        A1_temp     = ToIntervals(A1_vector');

        A1_temp(:,3) = A1_temp(:,2)-A1_temp(:,1)+1;
                A1_temp(:,4) = A1_temp(:,3) > cut_frame_num; % 4 frame difference means 100ms.
        
        A1_vector_cut = A1_vector;        
                
        for i = 1:size(A1_temp,1)
                if A1_temp(i,4) ==  0
                    A1_vector_cut(A1_temp(i,1):A1_temp(i,2)) = 0  ;
                end
        end
    
        compare_vectors = [A1_vector', A1_vector_cut'];
        temp_reshape = double(A1_vector_cut);
        
        ca_filt_less_fire = reshape(temp_reshape,size(A1,1),[]);
        
        % checking 
        pick_fire(:,1) = ca_filt_temp(:,1);
        pick_fire(:,2) = ca_filt_less_fire(:,1); 
end

%%