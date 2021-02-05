function [G, n, pc] = semi_nmf(dataset, mlset, ratios)
    n = 0;
    pc = 0;

    cd = ['data' filesep];
    str1 = [cd, dataset, '.mat'];
    str2 = [cd, mlset, '.mat'];
    % for un-weighted network
    [G, W, n, pc] = createMatrix(str1, str2);
    % for weighted network ()
    %     [G,W,n,pc] = createMatrix_w(str1,str2);

    for i = 1:length(ratios)
        k = ratios(i);

        WW(G, k, W);
    end

end

function WW(G, k, W)
    nClass = k;

    n = size(G, 1);

    %initilize V
    [V, new_matrix] = K_rank_D_new(G, nClass);
    vec = V';
    vec = matrix2norm(vec, 2);
    V = vec';

    %initilize U
    [n, m] = size(V);
    RandStream.setGlobalStream(RandStream('mt19937ar', 'seed', sum(100 * clock)));
    U = rand(m, m);
    U(:, :) = 1 + U(:, :) / 10; %high entropy

    lambda = 10;

    [X, U] = SSNMTF(G, nClass, W, U, V, lambda);

    S = X;

    a = strcat(num2str(lambda), '.mat');
    % save(a,'S');

    %non-overlapping: choose the maximum as the community ID
    [tmp, labels] = max(S, [], 2);

    l_file = strcat('data/non-overlapping.txt');
    write_infile(l_file, labels);

    %overlapping
    community_detect(S);
    l_file = strcat('data/overlapping.txt');
    write_infile(l_file, labels);
end

function [mX, U] = SSNMTF(network, nClass, W, U, V, lambda)
    options = [];
    options.maxIter = 1000;

    options.alpha = 1;
    options.nRepeat = 1;
    options.V = V;
    options.U = U;
    options.lambda = lambda;

    [U, V] = GNMF_SGD(network, nClass, W, options);
    mX = V;
end
