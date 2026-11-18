function P = bipartite_PC(B,M)

SumB = sum(B,2);

P = zeros(length(SumB),1);
for m = 1:max(M)
P = P+(sum(B(:,M==m),2)./SumB).^2;
end
P = 1-P;