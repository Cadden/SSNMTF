%create a matrix for a undirected network
% n is the number of veticles in network, netFile is the network file
function ml = createMatrix_ml(x, y)

ml = x;
len = length(y);
for i = 1:len
    
    ml(y(i,1) + 1, y(i,2) + 1) = 2;
    ml(y(i,2) + 1, y(i,1) + 1) = 2;    
    
end