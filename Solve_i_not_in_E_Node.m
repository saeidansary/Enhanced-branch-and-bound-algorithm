function [cands_new,vcands_new,scands_new,U,cur_sol]=Solve_i_not_in_E_Node(A,b,Q0,b0,b2,c2,cands,vcands,scands,cands_new,vcands_new,scands_new,U,cur_sol,node,nodes_new,j)
%Solving a general node of the type [i 2k+1]
% In those nodes one needs to check for candidates from previous layer i-1
% whether they satisfy the ith constraint. In addition, it also takes the
% candidates from the node [i 2k], which is the node where the ith
% constraint is in E (where in the [i 2k+1] node it is not).
i=node(1);
k=(node(2)-1)/2;
%node=[i 2*k+1];%node(1)+1==i+1, node(2)+1==2k+2
KN=(nodes_new==2*k+1);
if nnz(j==k)>0&&size(cands{j==k},2)>0
    feas=(A(i,:)*cands{j==k}<b(i)*ones(1,size(cands{j==k},2))+1e-4);
    cands_new{KN}=[cands_new{KN},cands{j==k}(:,feas)];
    vcands_new{KN}=[vcands_new{KN},vcands{j==k}(feas)];
    l=size(cands{j==k}(:,feas),2);
else
    l=0;
end
s=scands_new(nodes_new==2*k);%node(1)+1==i+1, node(2)-1+1==2k+1
cands_new{KN}=[cands_new{KN},cands_new{nodes_new==2*k}];
vcands_new{KN}=[vcands_new{KN},vcands_new{nodes_new==2*k}];
if min(size(vcands_new{KN}))>0
    [val,ind]=min(vcands_new{KN});
    x=cands_new{KN}(:,ind);
    if A*x<b+1e-4
        if val<U
            cur_sol=x;
            U=val;
        end
    end
end
scands_new(KN)=l+s;
% 'cands_new' contains the local and global minimizers of the subproblems
% 'vcands_new' for their values and 'scands_new' for their numbers.
% 'cur_sol' and 'U' are updated when a better solution which is feasible
% for the original problem is detected.
%--------------------------------------------------------------
