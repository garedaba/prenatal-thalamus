
%% Load connectivity data

ImgMagLoc = '"C:\Program Files\ImageMagick-7.1.1-Q16-HDRI/magick.exe"';
%system([ImgMagLoc,' montage ./ResultFigs/CortMatu_thr* -tile 4x2 -geometry +0+0 ./ResultFigs/CortMatu.png']);

load('Data4Analysis.mat','Cconn','Deg','ThalNucName','ThalMatEst','norm_TCprop','RDeg','Thal2Ignore')
close all

IsValidMatuEst = ~isnan(ThalMatEst);
NValid = sum(IsValidMatuEst);
ThalMatEstValid = ThalMatEst(IsValidMatuEst);

if exist("./NucRandPerm.mat","file")==2
    load("./NucRandPerm.mat","NucPerm")
else
    NucPerm = zeros(NValid,1000);
    for i = 1:1000
    NucPerm(:,i) = randperm(NValid,NValid);
    end
end

%%

gL = gifti('./dHCP40wk_surface/week-40_hemi-right_space-dhcpSym_dens-32k_vinflated.surf.gii');
verts = double(gL.vertices);
verts4centroid = verts;
faces = double(gL.faces);

parc = dlmread('dHCP_wk40_100random_parc.txt');

nparc = max(parc);
lh_centroid_vertIND = zeros(nparc, 1);
lh_centroid_coord = zeros(nparc, 3);

G = GfromSurf(verts,faces,parc);
for i = 1:nparc

    InParc = find(parc == i);
    g = subgraph(G,InParc);    
    d = distances(g);
    meanD = mean(d);    
    [~, I] = min(meanD);
    lh_centroid_vertIND(i) = InParc(I);
    lh_centroid_coord(i,:) = verts4centroid(lh_centroid_vertIND(i), :);
end

%subplot(2,6,[1:3 7:9])

[~,ThalOrd] = sort(ThalMatEst,'ascend');

[~,CortOrd] = sort(lh_centroid_coord(:,2),'descend');

PredAgeMat = repmat(ThalMatEst,[1 100]);

W = PredAgeMat.*double(norm_TCprop>.5);

ThalColours = map_values_to_colormap(ThalMatEst,flipud(coolwarm));
ThalColours(isnan(ThalMatEst), :) = repmat([.5 .5 .5], sum(isnan(ThalMatEst)), 1);
%CortColors = map_values_to_colormap(lh_centroid_coord(:,2),viridis);

figure('Position',[0 0 1920 1767])
clf
ConGramPos = [0.0134259259259259 0.44-.28 0.466203703703704 0.551025991792065];

ConGramPosOld = [0.0134259259259259 0.44 0.466203703703704 0.551025991792065];
AxisForaPanel = axes('Position',ConGramPosOld);
axis square
axis off
AxisForaPanel.TitleHorizontalAlignment = 'left';
Titleaa = title("\fontsize{28}a | thalamocortical connectivity examples");

ConnGramAxBack = axes('Position',ConGramPos);
axis square
axis off
MissingNucMat = double(isnan(PredAgeMat));
plot_thalamocortical_connectogram(MissingNucMat, ...
'orderA', ThalOrd, ...
'orderB', CortOrd, ...
'Asize', 10, ...
'Bsize', 10, ...
'thal_frac',.2,...
'cortex_frac',.75,...
'EdgeAlpha', norm_TCprop.^3, ...
'Acolours', ThalColours, ...
'Bcolours', [0.3, 0.3, 0.3], ...
'Shuffle', 1);
colormap(ConnGramAxBack,[.5 .5 .5])
axis off

ConnGramAx = axes('Position',ConGramPos);
ConnGramAx.Color='None';
axis square
axis off
pp = plot_thalamocortical_connectogram(PredAgeMat, ...
'orderA', ThalOrd, ...
'orderB', CortOrd, ...
'Asize', 100, ...
'Bsize', 140, ...
'thal_frac',.2,...
'cortex_frac',.75,...
'EdgeAlpha', norm_TCprop.^3, ...
'Acolours', ThalColours, ...
'Bcolours', [0.3, 0.3, 0.3], ...
'Shuffle', 1);
colormap(ConnGramAx,flipud(coolwarm))
axis off

connGramCbar = colorbar(ConnGramAx,'southoutside');
connGramCbar.Label.String = 'Maturation index';
connGramCbar.Label.FontSize=16;
connGramCbar.Position = [0.0141203703703704 0.440944293604344-.28 0.464814814814815 0.0120754716981132];
set(gca,'FontSize',16)

ConnGramAx.TitleHorizontalAlignment = 'left';
Titlea = title("\fontsize{28}b | thalamocortical connectogram");


%% Correlation between total sum of nuclei's connectivity and maturation

NucTotalStr = sum(norm_TCprop,2);

[StrMatu_c,StrMatu_p]=corr(ThalMatEst,NucTotalStr,'rows','pairwise','type','Pearson');

n = sum(IsValidMatuEst);
z = atanh(StrMatu_c);
SE = 1 / sqrt(n - 3);
% 95% confidence intervals
StrMatu_CI = tanh(z + [-1 1] * 1.96 * SE);

StrMatu_perm = zeros(1,1000);

for k = 1:1000
PermThalage = ThalMatEstValid(NucPerm(:,k));
StrMatu_perm(1,k)  = corr(NucTotalStr(IsValidMatuEst),PermThalage,'rows','pairwise','type','Pearson');
end

r_emp = StrMatu_c;
r_null = StrMatu_perm;
%two sided permutation test
StrMatu_perm_p = (sum(abs(r_null) >= abs(r_emp)) + 1) / (numel(r_null) + 1);

fprintf('Total strength-maturation result:\n');
fprintf('Correlation: %.2f\n', StrMatu_c);
fprintf('p-value (perm): %.3g\n', StrMatu_perm_p);
fprintf('Confidence interval: [%.2f, %.2f]\n', StrMatu_CI);

%% Correlation between the number of connections a nuclei has above a given threshold, and matuation
TCThrNConn_ax = axes('Position',[.525    0.511    0.2260    0.18]);

thrs = 0:.01:1;
for i = 1:length(thrs)
    T_conn_thr(:,i) = sum(norm_TCprop>=thrs(i),2);
end

NucMatColours = map_values_to_colormap(ThalMatEst,flipud(coolwarm));
NucMatColours(isnan(ThalMatEst),:) = repmat([.5 .5 .5],[sum(isnan(ThalMatEst)) 1]);
hold on

[strThr_corr,strThr_p]=corr(ThalMatEst,T_conn_thr,'rows','pairwise','type','Pearson');

CI_strThr = zeros(length(thrs),2);

n = sum(IsValidMatuEst);
for i = 1:length(thrs)
% Fisher z-transformation
z = atanh(strThr_corr(i));
SE = 1 / sqrt(n - 3);
% 95% confidence intervals
CI_strThr(i,:) = tanh(z + [-1 1] * 1.96 * SE);
end

for i = 1:35
plot(thrs,T_conn_thr(i,:),'LineWidth',2,'Color',NucMatColours(i,:))
end
ylabel('Number of connections')
xlabel('Strength threshold')
set(gca,'FontSize',16)
[strThr_p_pFDR_isSig,~,strThr_p_pFDR] = BF_FDR(strThr_p,.05,1);

alpha = 0.05;

yl = ylim;

% Logical mask of significance
sigMask = strThr_p_pFDR <= alpha;

% Find contiguous significant segments
d = diff([0 sigMask' 0]);
startIdx = find(d == 1);
endIdx   = find(d == -1) - 1;

for i = 1:length(startIdx)
    xPatch = [thrs(startIdx(i))-.005 thrs(endIdx(i))+.005 ...
              thrs(endIdx(i))+.005 thrs(startIdx(i))-.005];
    yPatch = [yl(1) yl(1) yl(2) yl(2)];
    
    patch(xPatch, yPatch, [255, 0, 0]./255, ...
          'EdgeColor','none', ...
          'FaceAlpha',0.4);
end

% Plot the lines of each nuclei coloured according to their maturation
for i = 1:35
plot(thrs,T_conn_thr(i,:),'LineWidth',2,'Color',NucMatColours(i,:))
end
ylabel('Number of connections')
xlabel('Strength threshold')

box on

strThr_corr_dup = strThr_corr;
strThr_corr_dup(~sigMask) = NaN;
[~,BestCorr] = nanmin(strThr_corr_dup);

fprintf('Best strength threshold result:\n');
fprintf('Threshold: %.2f\n', thrs(BestCorr));
fprintf('Correlation: %.2f\n', strThr_corr_dup(BestCorr));
fprintf('p-value: %.3g\n', strThr_p_pFDR(BestCorr));
fprintf('Confidence interval: [%.2f, %.2f]\n', CI_strThr(BestCorr,:));

TCThrNConn_ax.TitleHorizontalAlignment = 'left';
Titleb = title("\fontsize{28}c | less mature nuclei have more diffuse connectivity");

%subplot(2,6,6)
TCThrNConnBest_ax = axes('Position',[.8    0.511    0.19    0.18]);

scat = scatterfit(ThalMatEst,T_conn_thr(:,BestCorr),160,ThalMatEst,'',1);
xlabel('Maturation index')
ylabel('Number of connections')
%title(['Strength threshold > ',num2str(thrs(BestCorr))])
scat.MarkerEdgeColor='k';
set(gca,'FontSize',16)
colormap(flipud(coolwarm))

%%

DegThr = prctile(Deg,90);
Rinds = (length(Cconn)/2)+1:length(Cconn);

ValName = {'Hub conn. strength',...
           'Non-hub conn. strength'...
           'Participation coefficient'};

IsHub = double(RDeg>DegThr)+1;

P = bipartite_PC(norm_TCprop,IsHub);
Conn2Hub = sum(norm_TCprop(:,IsHub==2),2);
Conn2NonHub = sum(norm_TCprop(:,IsHub==1),2);

%% Real correlations (rho + parametric p-value)
[HubMatu_corr(3), HubMatu_corr_p(3)] = ...
    corr(P,ThalMatEst,'rows','pairwise','type','Pearson');

[HubMatu_corr(1), HubMatu_corr_p(1)] = ...
    corr(Conn2Hub,ThalMatEst,'rows','pairwise','type','Pearson');

[HubMatu_corr(2), HubMatu_corr_p(2)] = ...
    corr(Conn2NonHub,ThalMatEst,'rows','pairwise','type','Pearson');

HubMatu_corr_perm = zeros(3,1000);
HubMatu_corr_perm_p = zeros(3,1);
HubMatu_CI = zeros(3,2);

n = sum(IsValidMatuEst);
for i = 1:3
% Fisher z-transformation
z = atanh(HubMatu_corr(i));
SE = 1 / sqrt(n - 3);
% 95% confidence intervals
HubMatu_CI(i,:) = tanh(z + [-1 1] * 1.96 * SE);
end

for k = 1:1000

PermThalage = ThalMatEstValid(NucPerm(:,k));

HubMatu_corr_perm(3,k) = corr(P(IsValidMatuEst),PermThalage,'rows','pairwise','type','Pearson');
HubMatu_corr_perm(1,k)  = corr(Conn2Hub(IsValidMatuEst),PermThalage,'rows','pairwise','type','Pearson');
HubMatu_corr_perm(2,k) = corr(Conn2NonHub(IsValidMatuEst),PermThalage,'rows','pairwise','type','Pearson');
end

for i = 1:3
    r_emp = HubMatu_corr(i);
    r_null = HubMatu_corr_perm(i,:);
%two sided permutation test
HubMatu_corr_perm_p(i) = (sum(abs(r_null) >= abs(r_emp)) + 1) / (numel(r_null) + 1);
end

HubMatuType = {'Hub connectivity-maturation','Non-hub connectivity-maturation','Hub participation-maturation'};

for i = 1:3
fprintf([HubMatuType{i},' result:\n']);
fprintf('Correlation: %.2f\n', HubMatu_corr(i));
fprintf('p-value: %.3g\n', HubMatu_corr_perm_p(i));
fprintf('Confidence interval: [%.2f, %.2f]\n', HubMatu_CI(i,:));
end

% pos = get(gca,'Position');
% ax = gca;
% ax.Position(1) = pos(1)+.05;
% ax.Position(3) = pos(3)-.02;
% ax.Position(2) = pos(2)+.01;

%ax.XRuler.TickLabelGapOffset = -3;

Scat1=axes('Position',[.525    0.324    0.2    0.117823529411761]);
%scatter(Conn2Hub,ThalMatEst,40,ThalMatEst,'filled','MarkerEdgeColor','k')
scat=scatterfit(Conn2Hub,ThalMatEst,160,ThalMatEst,'',1);
scat.MarkerEdgeColor='k';
colormap(Scat1,flipud(coolwarm))
xlabel('Hub connectivity')
ylabel('Maturation index')
set(gca,'FontSize',16)
Titlec = title("\fontsize{28}d | maturation is unrelated to thalamic-hub connectivity");

Scat15=axes('Position',[.525    0.166    0.2    0.117823529411761]);
%scatter(Conn2Hub,ThalMatEst,40,ThalMatEst,'filled','MarkerEdgeColor','k')
scat=scatterfit(Conn2NonHub,ThalMatEst,160,ThalMatEst,'',1);
scat.MarkerEdgeColor='k';
colormap(Scat15,flipud(coolwarm))
xlabel('Non-hub connectivity')
ylabel('Maturation index')
set(gca,'FontSize',16)

Scat2=axes('Position',[.78    0.166    0.2    0.273]);

%scatter(P,ThalMatEst,40,ThalMatEst,'filled','MarkerEdgeColor','k')
scat=scatterfit(P,ThalMatEst,160,ThalMatEst,'',1);
scat.MarkerEdgeColor='k';
colormap(Scat2,flipud(coolwarm))
xlabel('Participation coefficient')
ylabel('Maturation index')
set(gca,'FontSize',16)

Titlea.Position(1) = -1;
Titleb.Position = [-.18 257.775688880164];
Titlec.Position = [7.640 0.861207625836436];
%Titleaa.Position(1) = -1;

%% Uncomment below to save the figure
%exportgraphics(gcf,'./figures/Fig3Partial.png','Resolution','300')