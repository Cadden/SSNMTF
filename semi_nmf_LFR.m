function [] = semi_nmf_LFR(dataset, type, alpha, ratios, nRpt)

    if type == 0% non-overlapping LFR network
        [n, M, nClass, labels, rlabels] = load_real(dataset);
    elseif type == 1% overlapping LFR network
        [n, M, nClass, labels, rlabels] = load_real_ov(dataset);
    end

    % initialize F and G
    m2 = size(M, 1);
    G = sparse(M(:, 1), M(:, 2), ones(m2, 1));
    G = full(G);

    %initilize V (F in SSNMTF)
    [V, new_matrix] = K_rank_D_new(G, nClass);
    vec = V';
    vec = matrix2norm(vec, 2);
    V = vec';

    %initilize U (G  in SSNMTF)
    [n, m] = size(V);
    RandStream.setGlobalStream(RandStream('mt19937ar', 'seed', sum(100 * clock)));
    U = rand(m, m);
    U(:, :) = 1 + U(:, :) / 10; %high entropy

    for i = 1:length(ratios)
        ratio = ratios(i);
        tt = datestr(now, 0);
        disp([tt, '    semi-supervised nmf on ' dataset ' with ratio ' num2str(ratio)]);

        if type == 0% non-overlapping LFR network
            [W] = sample_pair(n, rlabels, ratio, 1);
        elseif type == 1% overlapping LFR network
            [W] = sample_pair_ov(n, rlabels, labels, ratio, 1);
        end

        WW(G, nClass, W, V, U);

    end

end

function WW(G, k, W, V, U)
    nClass = k;

    n = size(G, 1);

    lambda = 10;

    [X, U] = SSNMTF(G, nClass, W, U, V, lambda);

    S = X;

    a = strcat(num2str(lambda), '.mat');
    save(a, 'S');

    % choose the maximum as the community ID
    [tmp, labels] = max(S, [], 2);

    a = strcat(num2str(lambda), '.mat');
    save(a, 'labels');

    a = strcat(num2str(lambda), '.mat');
    save(a, 'U')

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
