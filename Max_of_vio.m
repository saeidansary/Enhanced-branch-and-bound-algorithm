% Computes a constraint (ind>i) with maximum violating candidates number
function ind=Max_of_vio(A,b,i,cands)
%n=size(A,2);
m=size(A,1);
%num=zeros(m-i,1);
K=size(cands,2);
cands_ac=cands{1};
if i<m
    for k=1:K
        cands_ac=union(cands{k}',cands_ac','rows')';
    end
    c=size(cands_ac,2);
    num=sum((A(i+1:m,:)*cands_ac>b(i+1:m)*ones(1,c))-1e-2,2);
    [~,ind]=max(num);
    ind=ind+i;
end