function [n, M, nClass, labels, rlabels] = load_real_ov(dataset)
cd =['data' filesep ];

str1=[cd,dataset, '.dat']; %network filename
str2=[cd,dataset, '_rlabels.mat']; %network filename

load(str2);
n = length(labels);

A=load(str1);
m=size(A,1);
M=sparse(A(:,1),A(:,2),ones(m,1),n,n);
M=full(M);

isolated_set = find(sum(M)==0);
selected_set = setdiff(1:n, isolated_set);

M = M(selected_set, selected_set); %remove isolated point

M=M+M';
[x,y,z]=find(M);
M=[x,y];

% labels = labels(selected_set); %remove isolated point
% n = length(labels);

% for overlapp modules
nClass = max(max(labels));
rlabels = cell(1,nClass);
[row cow] = size(labels);


for i = 1:row
   mou = find(labels(i,:)~=0);
   for j =1:length(mou)
       rlabels{1,labels(i,j)} = [rlabels{1,labels(i,j)} i];
   end
end
filename=[cd,dataset,'networkBenchmar.txt'];
write_cell_infile(rlabels,filename);
%for non-overlapp modules
% nClass = max(labels);
% rlabels = cell(1,nClass);
% [Y,I] = sort(labels);
% for i=1:nClass
%     rlabels{1,i} = I(Y==i);
% end

end

function write_cell_infile(rlabels,filename)
%  write matrix to file
fid=fopen(filename,'w');

[x,y]=size(rlabels);

for i=1:x
    for j=1:y
        member = rlabels{i,j};
        for m=1:length(member)
            fprintf(fid,'%d\t',member(m));
        end
        fprintf(fid,'%d\n', '');%每一行回车\n
    end
    
end

fclose(fid);

end