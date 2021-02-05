function []=community_detect(S)
[n,m] = size(S);
S=convert3(S);
% S=convert2(S);

nodecom = diag(S' * S);

com = zeros(n,m);
for i = 1 : m
    
    % outlier
    if(nodecom(i) == 0)
        continue;
    end
    
    comsize = nodecom(i);
    [strength, index]=sort(S(:,i),'descend');
    sum  = 0;
    jth = 0; % jth protein
    while sum < comsize
        jth = jth +1;
        sum = sum + strength(jth);
    end
    
    for proth =1 : jth
        [~,location] = find(com(index(proth),:) == 0);
        com(index(proth),location(1)) = i;
    end
    
end

labels = com;
l_file = strcat('data/overlapping.txt');
write_infile( l_file,labels );

end

function S=convert3(X)
b= sqrt(sum(X.^2,2));
% b=power(sum(X,2),2);
% b=sum(X,2);
c=size(X,2);
D=repmat(b, [1,c]);
S=X./(D+eps);
end
