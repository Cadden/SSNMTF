% clear;

methods = {'trinmf'};
k = [10];
for m = 1 : length(methods)
    for j =1 : length(k)
%         file = strcat(methods{m},num2str(k(j)));
%         read = load(file);
%         labels = read.labels;
        
%         l_file = strcat('E:/project/GOAnalysis/material/methods/trinmf10000/', methods{m},num2str(k(j)),'.txt');%write lables into a new txt file
        l_file = strcat('E:/project/GOAnalysis/material/methods/trinmf10000/trinmfol500.txt');
        write_infile( l_file,labels );
    end

end

clear;
