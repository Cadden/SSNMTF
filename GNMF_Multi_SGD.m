function [U_final, V_final, nIter_final, objhistory_final] = GNMF_Multi(X, k, W, options, U, V)

differror = options.error;
maxIter = options.maxIter;
nRepeat = options.nRepeat;
minIter = options.minIter - 1;
if ~isempty(maxIter) && maxIter < minIter
    minIter = maxIter;
end
meanFitRatio = options.meanFitRatio;

alpha = options.alpha;

Norm = 2;
NormV = 0;

[mFea,nSmp]=size(X);

if alpha > 0
    W = W;
    DCol = full(sum(W,2));
    D = spdiags(DCol,0,nSmp,nSmp);
    L = D - W;
    if isfield(options,'NormW') && options.NormW
        D_mhalf = spdiags(DCol.^-.5,0,nSmp,nSmp) ;
        L = D_mhalf*L*D_mhalf;
    end
elseif alpha < 0
    L = (-alpha)*W;
else
    L = [];
end

selectInit = 1;
if isempty(U)
    U = abs(rand(k,k));

%     U = 0.5 * (U + U');
%     U = tril(U,-1)+triu(U',0);
    
%     U(logical(eye(size(U))))=1;
    
%     V = abs(rand(nSmp,k));
      V = options.V;
else
    nRepeat = 1;
end
% V = options.V;
% if isempty(V)
%     V = abs(rand(nSmp,k));
% else
%     V = options.V;
% end
%归一化U & V
[U,V] = NormalizeUV(U, V, NormV, Norm);

eta = options.eta;%1e-8
gama = options.gama;
lambda = options.lambda;

if nRepeat == 1
    selectInit = 0;
    minIter = 0;
    if isempty(maxIter)
        objhistory = CalculateObj(X,V,U,L,eta,lambda, gama);
        meanFit = objhistory*10;
    else
        if isfield(options,'Converge') && options.Converge
            objhistory = CalculateObj(X,V,U,L,eta,lambda, gama);
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

while tryNo < nRepeat   
    tryNo = tryNo+1;
    maxErr = 1;
    while(maxErr > differror)
        % ===================== update V (Z)========================              
        vu = V * U;
        wv = W * V;
        dv = D * V;
       
        VUU = 2 * V * (vu' * vu) + lambda * dv + gama * V;
        
        V = V.*((2*X * vu + lambda * wv)./max(VUU,1e-10));
        
        [U,V] = NormalizeUV(U, V, NormV, Norm);
        
%        [V,obj_new,learning_rate] = update_z_armijo(X,V,U,L,eta,lambda,learning_rate1);
        
        % ===================== update U ========================
        XV = V' * X * V;
        vv = V' * V;
        UVV = vv * U * vv + eta * U;
        U = U.*(XV./max(UVV,1e-10)); % 3mk
    
%         [U,obj_new,learning_rate2] = update_u_armijo(X,V,U,L,eta,lambda,learning_rate2);
        
        nIter = nIter + 1;
        if nIter > minIter
            if selectInit
                objhistory = CalculateObj(X,V,U,L,eta,lambda, gama);
%                 objhistory = obj_new;
                maxErr = 0;
            else
                if isempty(maxIter)
                    newobj = CalculateObj(X,V,U,L,eta,lambda, gama);
                    objhistory = [objhistory newobj]; %#ok<AGROW>
                    meanFit = meanFitRatio*meanFit + (1-meanFitRatio)*newobj;
                    maxErr = (meanFit-newobj)/meanFit;
                else
                    if isfield(options,'Converge') && options.Converge
                        newobj = CalculateObj(X,V,U,L,eta,lambda, gama);
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
            U = abs(rand(k,k));      
            V = abs(rand(nSmp,k));
%              V = options.V;
       
            [U,V] = NormalizeUV(U, V, NormV, Norm);
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
end

[U_final,V_final] = NormalizeUV(U_final, V_final, NormV, Norm);


%==========================================================================

function [obj] = CalculateObj(A,Z,U,L,eta,lambda,gama)

zu = Z * U;
zuz = zu * Z';
zlz = Z' * L * Z;

obj = sum(sum((A - zuz).^2)) + lambda * trace(zlz) + eta * sum(sum(U.^2)) + gama * sum(sum(Z.^2));

function [U, V] = NormalizeUV(U, V, NormV, Norm)

   %对V(列)进行归一化
    K = size(V,2);  
    N = size(V,1);
%     
% %     column
    norms = max(1e-15,sqrt(sum(V.^2,1)));
    normV = repmat(norms.^-1,N,1);
    V = V .* normV;
   
    normsU = repmat(norms',1,K);
    U = U .* normsU;
    
    U = U .* repmat(norms,K,1);