function [U_final, V_final, nIter_final, objhistory_final] = GNMF_Multi_SGD_ZU(X, k, W, options, U, V)

differror = options.error;
maxIter = options.maxIter;
nRepeat = options.nRepeat;
minIter = options.minIter - 1;
if ~isempty(maxIter) && maxIter < minIter
    minIter = maxIter;
end
meanFitRatio = options.meanFitRatio;

[mFea,nSmp]=size(X);

lambda = options.lambda;

W = lambda * W;
DCol = full(sum(W,2));
D = spdiags(DCol,0,nSmp,nSmp);
L = D - W;

selectInit = 1;
U = options.U;
V = options.V;


if nRepeat == 1
    selectInit = 0;
    minIter = 0;
    if isempty(maxIter)
        objhistory = CalculateObj(X,V,U,L);
        meanFit = objhistory*10;
    else
        if isfield(options,'Converge') && options.Converge
            objhistory = CalculateObj(X,V,U,L);
        end
    end
else
    if isfield(options,'Converge') && options.Converge
        error('Not implemented!');
    end
end



tryNo = 0;
nIter = 0; 

% learning_rate1 = 0.5;
% learning_rate2 = 0.7;
ov =[];

while tryNo < nRepeat
    tryNo = tryNo+1;
    maxErr = 1;
    while(maxErr > differror)
        % ===================== update V (Z)========================
        vu = V * U;
        wv = W * V;
        dv = D * V;
        uu = U * U;
        
        
        VUU = 2*(V * (vu' * vu)) + dv * uu+ dv;
        
        V = V.*((2*(X * vu) + wv * uu + wv)./max(VUU,1e-10));

        
        % ===================== update U ========================
        XV = (V' * X) * V;
        uv = U * V';
        vv = V' * V;
        wv = W * V;
        dv = D * V;
        UVV = vv * (U * vv) + uv * dv;
        U = U.*((XV + uv * wv)./max(UVV,1e-10)); % 
        
%          HeatMap(U);
%         ov = [ov CalculateObj(X,V,U,L)];
        nIter = nIter + 1;
        if nIter > minIter
            if selectInit
                objhistory = CalculateObj(X,V,U,L);
                %                 objhistory = obj_new;
                maxErr = 0;
            else
                if isempty(maxIter)
                    newobj = CalculateObj(X,V,U,L);
                    objhistory = [objhistory newobj]; %#ok<AGROW>
                    meanFit = meanFitRatio*meanFit + (1-meanFitRatio)*newobj;
                    maxErr = (meanFit-newobj)/meanFit;
                else
                    if isfield(options,'Converge') && options.Converge
                        newobj = CalculateObj(X,V,U,L);
                        objhistory = [objhistory newobj]; %#ok<AGROW>
                    end
                    maxErr = 1;
                    if nIter >= maxIter
                        maxErr = 0;
                        if isfield(options,'Converge') && options.Converge
                        else
                            objhistory = 0;
                        end
                    end
                end
            end
        end
    end
    
    if tryNo == 1
        U_final = U;
        V_final = V;
        nIter_final = nIter;
        objhistory_final = objhistory;
    else
        if objhistory(end) < objhistory_final(end)
            U_final = U;
            V_final = V;
            nIter_final = nIter;
            objhistory_final = objhistory;
        end
    end
    
    if selectInit
        if tryNo < nRepeat
            %re-start
            U = options.U;
            V = options.V;

            nIter = 0;
        else
            tryNo = tryNo - 1;
            nIter = minIter+1;
            selectInit = 0;
            U = U_final;
            V = V_final;
            objhistory = objhistory_final;
            meanFit = objhistory*10;
        end
    end
%     ov;
end



%==========================================================================

function [obj] = CalculateObj(A,Z,U,L)

zu = Z * U;
zuz = zu * Z';
zlz = zu' * L * zu;

obj = sum(sum((A - zuz).^2)) + trace(zlz) + trace(Z'*L*Z);


