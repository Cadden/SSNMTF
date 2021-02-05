%ÇóÁÚ½Ó¾ØÕóG
function A=Adj(G)
n=max(unique(G));
A=sparse(n,n);
for i=1:size(G,1)
m=G(i,1);
n=G(i,2);
A(m,n)=1;
A(n,m)=1;
end
end