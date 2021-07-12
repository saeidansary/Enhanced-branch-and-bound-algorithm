function [cands,vcands,scands,U,cur_sol,x,node,val]=Solve_Root_Node(A,b,Q0,b0,b2,c2,cands,vcands,scands,U,cur_sol)
% In the root node one solves the standard TRS, ignoring the linear
% inequalities. ---------------------------------------------------
node=[0 0];
[x,val,ca,valc]=Lin_Reduced_dim_TRS(Q0,b0,b2,c2,[],[]);
%Utilizes the general procedure, but without linear constraints
s=size(x,2);
k=2^node(1)+node(2);
cands{k}=x;
if min(size(val))>0
    vcands{k}=val*ones(1,s);
else
    vcands{k}=inf;
end
scands(k)=s;
if size(A,1)>0
    if s>0&A*x(:,1)<b+1e-4
        cur_sol=x(:,1);
        U=val;
    end
    if s>1&A*x(:,2)<b+1e-4
        cur_sol=x(:,2);
        U=val;
    elseif size(ca,2)>0
        for i=1:size(ca,2)
            if A*ca(:,i)<b+1e-4
                [U,ind]=min([U,valc(i)]);
                if ind==2
                    cur_sol=ca(:,i);
                end
            end
            cands{k}=[cands{k},ca(:,i)];
            vcands{k}=[vcands{k},valc(i)];
            s=s+1;
        end
        scands(k)=s;
    end
else
    cur_sol=x(:,1);
    U=val;
end
% 'cands' contains the local and global minimizers of the TRS
% 'vcands' for their values and 'scands' for their numbers.
% 'cur_sol' and 'U' are updated when a better solution which is feasible
% for the original problem is detected.
%--------------------------------------------------------------