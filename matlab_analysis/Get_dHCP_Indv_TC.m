
% Thal2Ignore = [10 11 12 29 30];

files = dir('./SubTCedits_mats/*.csv');
names = {files.name};

ids = regexp(names, '^(.*)_rand\d+\.csv$', 'tokens', 'once');
ids = [ids{:}];   % convert from cell-of-cells to cell array of strings

TCs = zeros(40,100,length(names));
TCprop = zeros(40,100,length(names));

for i = 1:length(names)
a = dlmread(['./SubTCedits_mats/',names{i}]);
A = (a+a').*~eye(length(a));
tc = A(201:end,101:200);
if size(tc,1) == 40
TCs(:,:,i) = tc;
else
TCs(1:size(tc,1),:,i) = tc;
end
tc = TCs(:,:,i);
% Normalises each row (i.e., each nuclei) such that its connectivity equals
% 1. Indivcates the proportion of connectivity to each cortical ROI
TCprop(:,:,i) = tc./repmat(sum(tc,2),[1 size(tc,2)]);

end

ThalLUT = readtable('thalamic_all-regions_lut.csv');

TCprob_orig = TCprop;

TCprob_orig(isnan(TCprob_orig)) = 0; 

TCprop(Thal2Ignore,:,:) = [];
TCs_orig = TCs;
TCs(Thal2Ignore,:,:) = [];

TCprop(isnan(TCprop)) = 0; 

meanTCprop = mean(TCprop,3);

nconn = sum(TCs>0,3);

norm_TCprop = BF_NormalizeMatrix(meanTCprop,'scaledSigmoid');