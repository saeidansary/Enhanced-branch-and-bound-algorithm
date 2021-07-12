function [cands,vcands,scands,U,cur_sol,x,node,val]=Solve_i_in_E_Node(A,b,Q0,b0,b2,c2,cands,vcands,scands,U,cur_sol,node,nodes_new)
%Solving a genral node of the type [i 2k]
% In such nodes one needs to compute candidates of the problem 
% containing all the equality constraints which are feasible for the
% original, and then add al candidates from nodes in the set S_{i,E}^1.
i=node(1);
bin=de2bi(node(2));
KN=(nodes_new==node(2));
if length(bin)<i
    bin=[bin,zeros(1,i-length(bin))];
end
bin=fliplr(bin);
eq=find(1-bin);
in=find(bin);
%---------------
%%Computing the candidates in R_{i,E}:
[x,val,ca,valc,flag_GeSy,erg]=Lin_Reduced_dim_TRS(Q0,b0,b2,c2,A(eq,:),b(eq));
if flag_GeSy==1
    fprintf('Possibly inaccurate solution in node [%d,%d]\n',node(1),node(2))
    fprintf('Max violation of excluded candidates: %f\n',erg)
end
mmm=size(A,1);
if size(x,2)>0
    if sum(A*x(:,1)<b+1e-4)==mmm
        if val<U+1e-6
            cur_sol=x(:,1);
            U=val;
        end
    elseif size(x,2)>1& sum(A*x(:,2)<b+1e-4)==mmm
        if val<U+1e-6
            cur_sol=x(:,2);
            U=val;
        end
    end
end
s=size(x,2);
k=0;
if s>0
    if min(size(in))>0
        for j=1:s 
            if A(in,:)*x(:,j)<b(in)+1e-4
                cands{KN}=[cands{KN},x(:,j)];
                vcands{KN}=[vcands{KN},val];
                k=k+1;
            end
        end
    else
        for j=1:s 
            cands{KN}=[cands{KN},x(:,j)];
            vcands{KN}=[vcands{KN},val];
            k=k+1;
        end
    end
end
j=k;
if size(ca,2)>0
    for k=1:size(ca,2)
        if A*ca(:,k)<b+1e-4
            [U,ind]=min([U,valc(k)]);
            if ind==2
                cur_sol=ca(:,k);
            end
        end
        if min(size(in))>0
            if A(in,:)*ca(:,k)<b(in)+1e-4
                j=j+1;
                cands{KN}=[cands{KN},ca(:,k)];
                vcands{KN}=[vcands{KN},valc(k)];
            end
        else
            j=j+1;
            cands{KN}=[cands{KN},ca(:,k)];
            vcands{KN}=[vcands{KN},valc(k)];
        end
    end
end
%-----------
%Computing the set S_{i,E}^1 and adding the relevant candidates
[is,js]=S_iE_1(node,nodes_new); 
sp=zeros(size(js));
for t=1:length(js)
    if nnz(nodes_new==js(t))>0
        sp(t)=scands(nodes_new==js(t));
    end
    if sp(t)>0
        if size(cands{KN},2)>0
            CKN=[cands{KN};vcands{KN}];
            CKN=union(CKN',[cands{nodes_new==js(t)}(:,1:sp(t));vcands{nodes_new==js(t)}(1:sp(t))]','rows')';
            cands{KN}=CKN(1:end-1,:);
            vcands{KN}=CKN(end,:);
        else
            cands{KN}=[cands{KN},cands{nodes_new==js(t)}(:,1:sp(t))];
            vcands{KN}=[vcands{KN},vcands{nodes_new==js(t)}(1:sp(t))];
        end
        %j=j+sp(t)-repetitions
    end
end
j=size(vcands{KN},2);
scands(KN)=j;
[val,ind]=min(vcands{KN});
x=cands{KN}(:,ind);
% 'cands' contains the local and global minimizers of the subproblems
% 'vcands' for their values and 'scands' for their numbers.
% 'cur_sol' and 'U' are updated when a better solution which is feasible
% for the original problem is detected.
%--------------------------------------------------------------
