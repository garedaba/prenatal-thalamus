%% Load connectivity data

% Exclude some nuclei because they are not in the dorsal thalamus. This will be done for thalamic IDs 
Thal2Ignore = [10 11 12 29 30];

Get_dHCP_Indv_TC
ADJs = GetIndvRand100();

RInds = 101:200;

Cconn = MakeGroupAverage(ADJs,.2,.3);
Deg = sum(Cconn>0,2);

RDeg = Deg(RInds);
% 
DegThr = prctile(Deg,90);

ImgMagLoc = '"C:\Program Files\ImageMagick-7.1.1-Q16-HDRI/magick.exe"';

%% Load thalamic maturation

Thal2Ignore = [10 11 12 29 30];

ThalLUT = readtable('thalamic_all-regions_lut.csv');
ThalMatur = readtable('thalamic_nucleus_maturation_1.3thalamic_genes_0.3DS.csv');

ThalMatur.Properties.VariableNames(1) = "labels";
T_joined = join(ThalMatur, ThalLUT, 'Keys', 'labels');

ThalPredAge = zeros(40,1);
ThalMatEst = zeros(40,1);
for i = 1:40
    ind = T_joined.label_id==i;
    ThalPredAge(i) = mean(T_joined.predicted_age(ind));
    ThalMatEst(i) = mean(T_joined.maturation_index(ind));
end

ThalPredAge(Thal2Ignore)=[];
ThalMatEst(Thal2Ignore)=[];

ThalNucName = ThalLUT.labels(2:end);
ThalNucName(Thal2Ignore)=[];

IsValidMatuEst = ~isnan(ThalMatEst);

save('./data/Data4Analysis.mat','Cconn','Deg','ThalNucName','ThalMatEst','ThalPredAge','norm_TCprop','RDeg','Thal2Ignore','IsValidMatuEst')

mkdir ./figures