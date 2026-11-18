function ADJs = GetIndvRand100()

files = dir('./rand100_conns/*.csv');
names = {files.name};

ids = regexp(names, '^(.*)_rand\d+\.csv$', 'tokens', 'once');
ids = [ids{:}];   % convert from cell-of-cells to cell array of strings

ADJs = zeros(200,200,length(ids));

for i = 1:length(ids)
a = dlmread(['./rand100_conns/',ids{i},'_rand100.csv']);
A = (a+a').*~eye(length(a));
mu = dlmread(['./rand100_conns/',ids{i},'_mu.txt']);
ADJs(:,:,i) = A.*mu;
end
% 
% load('DemographicData.mat')
% 
% age = zeros(numel(ids),1);
% for i = 1:numel(ids)
%     parts = split(ids{i}, "_");
%     idx = strcmp(DemographicData.participant_id,parts{1}) & DemographicData.session_id == str2num(parts{2});
%     age(i) = DemographicData.ScanAge(idx);
% end
% 
% run = 0;
% if run
% 
% A = mean(ADJs,3);
% 
% % for i = 1:182
% % a = ADJs(:,:,i);
% % c(i) = corr(triu2vec(a),triu2vec(A));
% % end
% 
% Athr = strength_threshold(A,.1);
% for i = 1:182
% athr = strength_threshold(ADJs(:,:,i),.1);
% DegCorr(i) = corr(sum(Athr>0,2),sum(athr>0,2),'type','Spearman');
% end
% 
% gL = gifti('./ancillary_data/week-40_hemi-left_space-dhcpSym_dens-32k_vinflated.surf.gii');
% 
% parc = dlmread('dHCP_wk40_100random_parc.txt');
% 
% for i = 1:10
% athr = strength_threshold(ADJs(:,:,i),.1);
% deg = sum(athr>0,2);
% 
% ExampleSurfacePlotFunction({gL.vertices,gL.faces},parc,deg(101:200),turbo(264),'Degree');
% end
% 
% for i = 1:182
% athr = strength_threshold(ADJs(:,:,i),.1);
% deg = sum(athr>0,2);
% DEGS(:,i) = deg;
% end
% 
% IndvDegCorrs = corr(DEGS(101:200,:),'type','Spearman');
% 
% meanC = mean(IndvDegCorrs);
% 
% AgeDiff = abs(age-age');
% 
% scatter(AgeDiff(:),IndvDegCorrs(:))
% 
% ranked = tiedrank(DEGS(101:200,:));
% mean_rank = mean(ranked,2);
% 
% ExampleSurfacePlotFunction({gL.vertices,gL.faces},parc,mean_rank,turbo(264),'Mean degree rank');
% 
% Athrdeg = sum(Athr,2);
% figure
% scatter(Athrdeg(101:200),mean_rank);
% 
% end