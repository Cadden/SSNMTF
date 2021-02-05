%create a matrix for a undirected network
% n is the number of veticles in network, netFile is the network file
function [G,M,n,pc] = createMatrix(netFile,mustNetFile)
data = load(netFile);
data = data.protein;
% the number of proteins in PPI network
n = max(max(data)) + 1;

x = zeros(n);
len = length(data);
for i = 1:len
    x(data(i,1) + 1, data(i,2) + 1) = 1;
    x(data(i,2) + 1, data(i,1) + 1) = 1;
end

G = x;

data = load(mustNetFile);
data = data.W;

y = zeros(n);
len = length(data);
for i = 1:len
    y(data(i,1) + 1, data(i,2) + 1) = 1;
    y(data(i,2) + 1, data(i,1) + 1) = 1;
end
M = y;
% the number of proteins in complex database
pc = sum(sum(G,1)~=0,2);



