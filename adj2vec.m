%ͨ��signal_processing�����ڽӾ���ת��Ϊ�����ռ�
%adj ����G���ڽӾ���
function vec=adj2vec(adj)
n=size(adj,1);
v=adj+speye(n,n);
v=sparse(v);
vec=v*v*v;
%vec=vec-diag(diag(vec))+eye(n);

%�����б�׼���������ɵľ�����ÿ����Ϊһ�����������е�Ԫ��ֵ��ʾ��Ը����ڵ��Ӱ��̶�
% vec=matrix2norm(vec,2);
end