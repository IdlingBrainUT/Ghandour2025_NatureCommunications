%% Circshift shuffling
function shuffled_data = fxn_data_shuffling(data_input,iteration,shuffle_mode)
%% comment
% 211111: make 1st version
%% debug
% clc; clear; close all;
% % data_input = rand(100,2);
% data_input = [0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0]';
% iteration = 100;
%% 
[shift_max, shifted_column_num] = size(data_input); % time

%% circshift
%Y = circshift(A,K,dim)
rng = 1; 
data_shuffled = {};
%% data shuffling

% randperm shuffling
if shuffle_mode == 1
    disp('randperm shuffling')
    for i = 1:iteration   
        for ii = 1:shifted_column_num
            shift_each = randperm(shift_max);
            shuffled_data{i,1}(:,ii) = data_input(shift_each,ii);
        end
    end
    
% circshift shuffling
elseif shuffle_mode == 2
    disp('circshift shuffling')
    for i = 1:iteration
        shift_each = randi(shift_max, shifted_column_num, 1);
        for ii = 1:shifted_column_num
            shuffled_data{i,1}(:,ii) = circshift(data_input(:,ii), shift_each(ii), 1);
        end
    end
end
%%

%%
end
