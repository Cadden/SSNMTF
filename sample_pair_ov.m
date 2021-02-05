function [ W ] = sample_pair_ov(  nNode, rlabels, labels, percent, postive_para )
disp(['sampling pair at ratio:' num2str(percent)])
nClass = length(rlabels);
W = zeros(nNode);
nPair_total = 0;
for i=1:nClass
    class_member_single = rlabels{1,i};
    
    overlapped = [];
    for node=1:length(class_member_single)
        if(length(find(labels(class_member_single(node),:)~=0))>1)
            overlapped = [overlapped class_member_single(node)];
        end
    end    
    
    [c, ia, ib] = intersect(class_member_single,overlapped);
    class_member_single(ia) = [];
    
    class_size = length(class_member_single);
    nPairs = class_size * (class_size - 1)/2;
    nPair_total = nPair_total + nPairs*2;
    nSelected = ceil(nPairs*percent);
    u = randperm(nPairs);
    u=u(1:nSelected);
    pairs = zeros(nSelected,2);
    index = 1;
    for p=1:class_size;
        for q=1:p-1
            pairs(index,:) = [p q];
            index = index + 1;
        end
    end
    paris = pairs(u,:);
    
    for j =1:nSelected
        pair = paris(j,:);
        a = class_member_single(pair(1));
        b = class_member_single(pair(2));
        W(a,b) = postive_para;
        W(b,a) = postive_para;
    end
    
end
sum(sum(W))/nPair_total;

end


