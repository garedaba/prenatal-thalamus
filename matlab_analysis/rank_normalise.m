function y = rank_normalise(x)
% Rank-transform to approximately [0, 1).
% Tied values receive separate, consecutive ranks.

[~, order] = sort(x);
[~, ranks] = sort(order);

y = (ranks - 1) / numel(x) + 1e-8;
end