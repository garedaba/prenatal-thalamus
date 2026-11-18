%% Creates .csv files of the data to use 
% Run after CompileData is run

mkdir ./Data4GenMdl/

%% read in surface and parcellaiton data

parc = dlmread('dHCP_wk40_100random_parc.txt');
gL = gifti('./dHCP40wk_surface/week-40_hemi-left_space-dhcpSym_dens-32k_wm.surf.gii');
lh_verts = double(gL.vertices);
lh_faces = double(gL.faces);

gR = gifti('./dHCP40wk_surface/week-40_hemi-right_space-dhcpSym_dens-32k_wm.surf.gii');
rh_verts = double(gR.vertices);
rh_faces = double(gR.faces);
%% Left hemisphere
verts = lh_verts;
verts4centroid = verts;

nparc = max(parc);
lh_centroid_vertIND = zeros(nparc, 1);
lh_centroid_coord = zeros(nparc, 3);

G = GfromSurf(verts,lh_faces,parc);
for i = 1:nparc

    InParc = find(parc == i);
    %G_ = VFtoG(verts,lh_faces_noMedwall,parc,i);
    g = subgraph(G,InParc);    
    d = distances(g);
    meanD = mean(d);
    %[~, I] = min(meanD);

    %d = distances(G_,InParc,InParc);
    
    [~, I] = min(meanD);
    lh_centroid_vertIND(i) = InParc(I);
    lh_centroid_coord(i,:) = verts4centroid(lh_centroid_vertIND(i), :);
end

lh_GeoD = distances(G,lh_centroid_vertIND,lh_centroid_vertIND);

%% Right hemisphere
verts = rh_verts;
verts4centroid = verts;

nparc = max(parc);
rh_centroid_vertIND = zeros(nparc, 1);
rh_centroid_coord = zeros(nparc, 3);

G = GfromSurf(verts,rh_faces,parc);
for i = 1:nparc

    InParc = find(parc == i);
    %G_ = VFtoG(verts,lh_faces_noMedwall,parc,i);
    g = subgraph(G,InParc);    
    d = distances(g);
    meanD = mean(d);
    %[~, I] = min(meanD);

    %d = distances(G_,InParc,InParc);
    
    [~, I] = min(meanD);
    rh_centroid_vertIND(i) = InParc(I);
    rh_centroid_coord(i,:) = verts4centroid(rh_centroid_vertIND(i), :);
end

rh_GeoD = distances(G,rh_centroid_vertIND,rh_centroid_vertIND);

D = squareform(pdist([lh_centroid_coord;rh_centroid_coord]));

writematrix(Cconn,'./Data4GenMdl/Connectivity.csv')
writematrix(parc,'./Data4GenMdl/parc.csv')
writematrix(D,'./Data4GenMdl/Distance_euclidean.csv')
writematrix(lh_GeoD,'./Data4GenMdl/Left_geo_distance.csv')
writematrix(rh_GeoD,'./Data4GenMdl/Right_geo_distance.csv')

%% Get estimate of innervation timing

IsValidMatuEst = ~isnan(ThalMatEst);

ThalPredAge_mat = repmat(ThalPredAge,[1 100]);
WeiAvg_PredAge = nansum(norm_TCprop.*ThalPredAge_mat)./nansum(norm_TCprop(IsValidMatuEst,:));
writematrix(WeiAvg_PredAge,'./Data4GenMdl/cortical_innervation_PredAge.csv')

%% Write estimates to directory expected by the python code

writematrix(Cconn,'../data/connectivity/Connectivity.csv')
writematrix(parc,'../data/connectivity/parc.csv')
writematrix(D,'../data/connectivity/Distance_euclidean.csv')
writematrix(lh_GeoD,'../data/connectivity/Left_geo_distance.csv')
writematrix(rh_GeoD,'../data/connectivity/Right_geo_distance.csv')
writematrix(WeiAvg_PredAge,'./Data4GenMdl/cortical_innervation_timing.csv')

% ThalMatEst_mat = repmat(ThalMatEst,[1 100]);
% WeiAvg_MatEst = nansum(norm_TCprop.*ThalMatEst_mat)./nansum(norm_TCprop(IsValidMatuEst,:));
% writematrix(WeiAvg_MatEst,'./Data4GenMdl/cortical_innervation_MatuEst.csv')
