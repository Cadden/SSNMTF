% A is the adjcancy matrix, K is the real number of communities.

function [result,new_matrix]=K_rank_D_new(A,K)
if nargin < 2 % 进行初始化
    flag_need=0;
else
    flag_need=1;
end
ND=size(A,1);%节点数
new_matrix=adj2vec(A);
dis=batch_distance(new_matrix,new_matrix);
dis=dis-diag(diag(dis));
for i=1:ND %kernel
    N_inf(i)=0.;
end
N_inf2=pagerank_alg(A,0.85);%pagerank_new

if ND<1000
    flag=1;
else
    flag=2;
end
% Two cases: 1)
% with few nodes and pagerank value demonstrate litte divergence
if flag==1
    for i=1:ND-1
        for j=i+1:ND
            N_inf(i)=N_inf(i)+exp(-((dis(i,j)*dis(i,j))/(N_inf2(i))));
            N_inf(j)=N_inf(j)+exp(-((dis(i,j)*dis(i,j))/(N_inf2(j))));
        end
    end
    N_inf=N_inf+N_inf2';
else
    %otherwise 2)
    N_inf=N_inf2;
end
nneigh=zeros(1,ND);
delta=zeros(1,ND);
maxd=max(max(dis));
[N_sorted,ordN]=sort(N_inf,'descend');%降序排列，[value,index]
delta(ordN(1))=-1.;
nneigh(ordN(1))=0;
for ii=2:ND
    list=ordN(1:ii-1);
    distance=dis(ordN(ii),list);
    [value,index]=min(distance);
    delta(ordN(ii))=value;
    nneigh(ordN(ii))=list(index);
end
delta(ordN(1))=max(delta(:)); %author provided
if flag_need==0
    
    disp('Generated file:DECISION GRAPH')
    disp('column 1:Density')
    disp('column 2:Delta')
    
    disp('Select a rectangle enclosing cluster centers')
    
    tt=plot(N_inf(:),delta(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
    title ('Decision Graph of K-rank-D','FontSize',15.0)
    rect = getrect(1); %rect=[x_min,y_min,weight,height]
    rhomin=rect(1);
    deltamin=rect(2);
    NCLUST=0;
    for i=1:ND
        cl(i)=-1;
    end
    icl=[];
    for i=1:ND
        if ( (N_inf(i)>rhomin) && (delta(i)>deltamin))
            NCLUST=NCLUST+1;
            cl(i)=NCLUST; % 对类中心进行类别指派
            icl(NCLUST)=i;% store centers
        end
    end
    
     V = full(new_matrix(:,icl));
elseif flag_need==1
    %     N_inf_t=(N_inf-min(N_inf))/(max(N_inf)-min(N_inf));
%     delta_t=(delta-min(delta))/(max(delta)-min(delta));

    N_inf_t=N_inf / max(N_inf);
    delta_t=delta / max(delta);
    for i=1:ND
        ind(i)=i;
        gamma(i)=N_inf_t(i)* delta_t(i);
%           gamma(i)=(N_inf_t(i)* delta_t(i))/(N_inf_t(i)+ delta_t(i));
    end
    [~,index]=sort(gamma,'descend');
    NCLUST=K;% the number of clusters
    icl=index(1:NCLUST);
    V = full(new_matrix(:,icl));
%     save('V','V');
end
result=V;

