function [FR_frame_digi] = func_FR_input(session_frame_size, FR_posi_second)

%% Murayama freeze

% session_frame_size = 8000;
% FR_posi_second = [7,11  13,16  17,19  28,31  34,35  36,44  48,49  50,52  53,54  59,60  62,64  65,68  69,70  71,72  73,74  76,77  78,79  81,82  84,85  89,90  92,93  95,96  98,107  108,111  113,115  116,128  129,130  131,136  141,144  145,148  163,164  169,170  173,174  175,178  180,181  185,187  192,194  195,198  200,201  202,205  209,211  212,216  217,218  220,221  222,229  230,237  ];


        FR_posi_second = FR_posi_second*20;
        FR_frame1 = ones([session_frame_size,1])*10; mura_size = size(FR_posi_second,2);
                    for ii = 1:2:mura_size
                        frame = [FR_posi_second(ii):FR_posi_second(ii+1)];
                        FR_frame1(frame) = 0; 
                    end
                    
        FR_frame2 = FR_frame1 == 0;
        FR_frame_digi = double(FR_frame2);
                    %%
end