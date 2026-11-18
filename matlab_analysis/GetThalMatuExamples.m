fig = figure;
fig.Units = 'pixels';
fig.Position = [100 100 1920 1080]; 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
dpi = 100;  % You can choose this
fig.PaperPosition = [0 0 1920/dpi 1080/dpi];
fig.PaperSize     = [1920/dpi 1080/dpi];

cmap = flipud(coolwarm(256));

ThalColours = map_values_to_colormap(ThalMatEst,cmap);
ThalColours(isnan(ThalMatEst), :) = repmat([.5 .5 .5], sum(isnan(ThalMatEst)), 1);

[~,I] = nanmax(ThalMatEst);

PredAgeMat = repmat(ThalMatEst,[1 100]);

bbb = zeros(size(PredAgeMat));
bbb(I,:) = PredAgeMat(I,:);

[pp,sA,sB] = plot_thalamocortical_connectogram(bbb, ...
'orderA', ThalOrd, ...
'orderB', CortOrd, ...
'Asize', 100, ...
'Bsize', 140, ...
'thal_frac',.2,...
'cortex_frac',.75,...
'EdgeAlpha', norm_TCprop, ...
'Acolours', ThalColours, ...
'Bcolours', [0.3, 0.3, 0.3], ...
'Shuffle', 1);
colormap(cmap)
axis off
clim([nanmin(ThalMatEst) nanmax(ThalMatEst)])
xlimit=xlim;
ylimit=ylim;
% delete(sA)
% delete(sB)
axis equal
xlim(xlimit)
ylim(ylimit)
print(fig,['./figures/HIGH_MAT_EXAMPLE.png'],'-dpng',['-r' num2str(dpi)])


clf
[iI,II] = sort(ThalMatEst,'ascend');
I=II(5);
bbb = zeros(size(PredAgeMat));
bbb(I,:) = PredAgeMat(I,:);

[pp,sA,sB] = plot_thalamocortical_connectogram(bbb, ...
'orderA', ThalOrd, ...
'orderB', CortOrd, ...
'Asize', 100, ...
'Bsize', 140, ...
'thal_frac',.2,...
'cortex_frac',.75,...
'EdgeAlpha', norm_TCprop, ...
'Acolours', ThalColours, ...
'Bcolours', [0.3, 0.3, 0.3], ...
'Shuffle', 1);
colormap(cmap)
axis off
clim([nanmin(ThalMatEst) nanmax(ThalMatEst)])
xlimit=xlim;
ylimit=ylim;
% delete(sA)
% delete(sB)
axis equal
xlim(xlimit)
ylim(ylimit)
print(fig,['./figures/LOW_MAT_EXAMPLE.png'],'-dpng',['-r' num2str(dpi)])