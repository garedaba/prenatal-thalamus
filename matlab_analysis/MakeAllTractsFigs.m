outputDir = fullfile(pwd,'Thal_Streamline_Figures');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

figure('Position',[125.4444  136.1111  743.5556  538.6667])

axes('Position',[0.1   0.05  .8  .8])

ThalInd2Use = 1:40;
ThalInd2Use(ismember(ThalInd2Use,Thal2Ignore))=[];
Thal2Ignore = [10 11 12 29 30];
ThalLUT = readtable('thalamic_all-regions_lut.csv');

for i = 1:length(ThalNucName)

nucID = ThalInd2Use(i);
NucCOlour = [ThalLUT.R(nucID+1) ThalLUT.G(nucID+1)  ThalLUT.B(nucID+1)]./255;

AlphaData = norm_TCprop(i,:);

plot_streamlines_from_nucleus(nucID,NucCOlour,AlphaData.^3)

title(ThalNucName{i},'FontSize',36)

filename = fullfile(outputDir, sprintf('Nuc_%02d.png', i));
    %exportgraphics(gcf, filename, 'Resolution', 300);

print(filename,'-dpng','-r300')

cla
end

ImgMagLoc = '"C:\Program Files\ImageMagick-7.1.1-Q16-HDRI\magick.exe"';

inputPattern = fullfile(outputDir, 'Nuc_*.png');
outputFile   = fullfile(outputDir, 'All_Nuclei_Montage.png');

cmd = sprintf('%s montage "%s" -tile 5x -geometry +10+10 "%s"', ...
    ImgMagLoc, inputPattern, outputFile);

system(cmd);