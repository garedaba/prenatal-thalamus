outputDir = fullfile(pwd,'ThalConnCort_Figures');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

gR = gifti('./dHCP40wk_surface/week-40_hemi-right_space-dhcpSym_dens-32k_vinflated.surf.gii');

for i = 1:35

data = norm_TCprop(i,:);
ExampleSurfacePlotFunction({gR.vertices,gR.faces},parc,data,summer,[ThalNucName{i},' connectivity']);

filename = fullfile(outputDir, sprintf('Nuc_%02d.png', i));
print(filename,'-dpng', '-r300');
close all
end


ImgMagLoc = '"C:\Program Files\ImageMagick-7.1.1-Q16-HDRI\magick.exe"';

inputPattern = fullfile(outputDir, 'Nuc_*.png');
outputFile   = fullfile(outputDir, 'All_TCconn_Montage.png');

cmd = sprintf('%s montage "%s" -tile 5x -geometry +10+10 "%s"', ...
    ImgMagLoc, inputPattern, outputFile);

system(cmd);