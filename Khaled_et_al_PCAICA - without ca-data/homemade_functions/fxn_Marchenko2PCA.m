%% Maharanobis Population vector distance calculation code
% 1bin_ca_data: time x neuron matrix
function [PCA_X_MPsig_num, res_MPPCA] = fxn_Marchenko2PCA(ca_data, thrcov_PC_percent)
%% for debug,
% PC_threshold = 90; % 90 means PCs compose 90 percentile 
X = ca_data;
PC_threshold = thrcov_PC_percent;
%% MP-PCA dimension reduction
% (m)252 trading days of (n) 10 stock returns, 
% 10回の株式リターンの1年にわたるシリーズ（つまり、252取引日) square(1+sqrt(10/252)
% n × B matrix, (1+sqrt(n/B)).^2
% PCA calculation, [coeff,score,latent] = pca(X);

% Marchenko–Pastur cal
[size_X(1,1), size_X(1,2)] = size(X); % dim1 -> time, dim2 -> neuron
% [size_Y(1,1), size_Y(1,2)] = size(Y); % dim1 -> time, dim2 -> neuron

% (1+sqrt(n/t)).^2; % n -> neuron, t -> time;
[PCA_X_coeff, PCA_X_score, PCA_X_latent] = pca(X);
n1 = size_X(1,2); t1 = size_X(1,1);
MP_val_X = (1+sqrt(n1/t1)).^2; % n -> neuron, t -> time;
PCA_X_sig = find(PCA_X_latent > MP_val_X);
PCA_X_MPsig_num = numel(PCA_X_sig);

PCA_X_latent_cumsum(:,1) = PCA_X_latent;
PCA_X_latent_cumsum(:,2) = cumsum(PCA_X_latent);
PCA_X_latent_cumsum(:,3) = PCA_X_latent_cumsum(:,2)/sum(PCA_X_latent)*100;
PCA_X_latent_cumsum(:,4) = PCA_X_latent_cumsum(:,3) < PC_threshold;
PCA_X_PC_thr_num = sum(PCA_X_latent_cumsum(:,4),1); % original
% PCA_X_PC_thr_num = 3; % #### forced top3 PCs ON ###

MPPCA_X_score     = PCA_X_score(:,1:PCA_X_MPsig_num);
thr_cover_X_score = PCA_X_score(:,1:PCA_X_PC_thr_num);

figure('Position',[600,50,400,200]); 
% subplot(211); 
plot(PCA_X_latent,'k'); hold on
area(PCA_X_MPsig_num, max(PCA_X_latent),'EdgeColor', 'red','FaceColor', 'red'); 
title('Marchenko–Pastur Max Lambda for reference matrix'); ylabel('Eigenvalue'); xlabel('Latent')
area(PCA_X_PC_thr_num, max(PCA_X_latent),'EdgeColor','blue','FaceColor', 'blue'); 
title('Dimension reduction thresholding'); ylabel('Eigenvalue'); xlabel('Latent')
legend PCA-latent-eigenval  Marchenko–Pastur-thr PC-coverage-thr

ax = gca;
set(gca, 'FontSize', 10, 'FontName','Arial'); colormap('parula'); grid on; ax.TickDir = 'both';

% hold off
%%

%% 
res_MPPCA.MPPCA_sig_num          = PCA_X_MPsig_num;
res_MPPCA.MPPCA_sig_score           = MPPCA_X_score;
res_MPPCA.thrcov_PCA_latent_cumsum     = PCA_X_latent_cumsum;
res_MPPCA.MP_Lambda_max_val     = MP_val_X;

res_thrcov_PCA.thrcov_PCA_thr_num        = PCA_X_PC_thr_num;
res_thrcov_PCA.thrcov_PCA_score          = thr_cover_X_score;
res_thrcov_PCA.thrcov_PCA_latent_cumsum  = PCA_X_latent_cumsum;
res_thrcov_PCA.thrcov_PCA_threshold      =  PC_threshold;
%%
end
%%
