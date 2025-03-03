SigPtrnNo = [1:5,8:20,22:27,29,30];
SigPtrn = result_MPPCA.neuron_weight_sigID(:,SigPtrnNo);
[rows_SigPtrn,cols_SigPtrn] = find(SigPtrn == 1);
SeqSigPtrn = sort(rows_SigPtrn);
%% Excel
% % Export as an excel file
% xlswrite('SeqSigPtrn.xlsx',SeqSigPtrn)