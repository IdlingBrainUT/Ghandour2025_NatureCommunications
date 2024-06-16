FolderName = 'G:\Khaled\For_lab_men_MPPCA_ICA_V05_v11\PCA ICA results\KD109';   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
FigHandle = FigList(iFig);
FigName   = num2str(get(FigHandle, 'Number'));
set(0, 'CurrentFigure', FigHandle);
savefig(fullfile(FolderName, [FigName '.fig']));
end