function [pp,sA,sB] = plot_thalamocortical_connectogram(W,varargin)

p = inputParser;
addParameter(p,'orderA',[])
addParameter(p,'orderB',[])
addParameter(p,'EdgeAlpha',1)
addParameter(p,'Shuffle',0)
addParameter(p,'thal_frac',.3)
addParameter(p,'cortex_frac',.6)
addParameter(p,'Acolours',[0.2 0.4 0.8])
addParameter(p,'Bcolours',[0.8 0.3 0.3])
addParameter(p,'Asize',100)
addParameter(p,'Bsize',100)
parse(p,varargin{:})

orderA = p.Results.orderA;
orderB = p.Results.orderB;
EdgeAlpha = p.Results.EdgeAlpha;
Shuffle = p.Results.Shuffle;
thal_frac     = p.Results.thal_frac;
cortex_frac     = p.Results.cortex_frac;
Acolours = p.Results.Acolours;
Bcolours = p.Results.Bcolours;
Asize = p.Results.Asize;
Bsize = p.Results.Bsize;

% thal_frac = 0.3; % proportion of circle for thalamus
% cortex_frac = 0.6; % proportion of circle for cortex

[Nthal,Ncortex] = size(W);

if isempty(orderA), orderA = 1:Nthal; end
if isempty(orderB), orderB = 1:Ncortex; end

if size(Acolours,1)>1
    Acolours = Acolours(orderA,:);
end

if size(Bcolours,1)>1
    Bcolours = Bcolours(orderB,:);
end

W = W(orderA,orderB);


Circle = 2*pi;

CortProp = cortex_frac*Circle;
CortPropHalf = (CortProp)/2;
Circle_CortCent = pi/2;

cortex_theta = linspace(Circle_CortCent+CortPropHalf,Circle_CortCent-CortPropHalf,Ncortex);

ThalProp = thal_frac*Circle;
ThalPropHalf = ThalProp/2;
Circle_ThalCent = (pi/2)+pi;

thal_theta = linspace(Circle_ThalCent-ThalPropHalf,Circle_ThalCent+ThalPropHalf,Nthal);

x_thal = cos(thal_theta);
y_thal = sin(thal_theta);

x_ctx = cos(cortex_theta);
y_ctx = sin(cortex_theta);

%Wplot = W.*double(W>.7);
Wplot=W;

% ---- Plot edges ----
[x,y,w] = find(Wplot);

if Shuffle
    Nord = length(w);
    randEdgeOrd = randperm(Nord);
    x = x(randEdgeOrd);
    y = y(randEdgeOrd);
    w = w(randEdgeOrd);
end

LinePts = 25;

totalPoints = (length(w)*LinePts)+length(w);
allLines = nan(totalPoints,2);
LineVal = nan(totalPoints,1);
LineAlphaVal = nan(totalPoints,1);

if ~isscalar(EdgeAlpha)
    EdgeAlphaMat = EdgeAlpha(orderA,orderB);
else
    EdgeAlphaMat = W;
end

for k = 1:length(w)
t = linspace(0,1,LinePts);
% Bezier curve via centre
cx = (1-t).^2 * x_thal(x(k)) + 2*(1-t).*t*0 + t.^2 * x_ctx(y(k));
cy = (1-t).^2 * y_thal(x(k)) + 2*(1-t).*t*0 + t.^2 * y_ctx(y(k));
%lw = 0.5 + 2 * w(k) / max(w);

SegmentsPerLine = (ones(length(w),1)*LinePts)+1;

LineIndStart = [0;cumsum(SegmentsPerLine(1:end))];

CurrentPoints = LineIndStart(k)+1:LineIndStart(k+1)-1;

allLines(CurrentPoints,:) = [cx' cy'];

LineVal(CurrentPoints) = w(k);

LineAlphaVal(CurrentPoints) = EdgeAlphaMat(x(k),y(k));

end

%EdgeAlpha = .5;

if isscalar(EdgeAlpha)
    EdgeAlphaType = EdgeAlpha;
    FaceVertexAlphaData = LineVal;
else
    FaceVertexAlphaData = LineAlphaVal;
    EdgeAlphaType = 'flat';
end

pp = patch(allLines(:,1),allLines(:,2),LineVal,'FaceAlpha','flat','EdgeColor','flat','LineWidth',2,'FaceVertexAlphaData',FaceVertexAlphaData,'EdgeAlpha',EdgeAlphaType,'AlphaDataMapping','none');
alim([0 1])
hold on
sA = scatter(x_thal, y_thal, Asize,Acolours,'filled','MarkerEdgeColor','k');
sB = scatter(x_ctx, y_ctx, Bsize,Bcolours,'filled','MarkerEdgeColor','k');

% axis equal
% axis off