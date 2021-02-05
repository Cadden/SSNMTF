function [n, G, nClass, labels, rlabels] = load_file(dataset,type)
load(['data' filesep dataset '.mat']);
labels = labels_real;
end

