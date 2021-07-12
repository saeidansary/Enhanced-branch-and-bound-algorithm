function [U,cur_sol,nn,Lb,status]=BB_QBL_Heur(A,b,Q0,b0,b2,c2,N,gap)
% Solves an instance of problem (QBL):
% minimize    0.5*x'*Q0*x-b0'*x
% subject to  0.5*x'*x-b2'*x+0.5*c2<=0,
%                               A*x<=b.
% ----------------------------------------------------------
% The input should contain A,b,Q0,b0,b2,c2, where
% Q0 is a symmetric matrix in S^n,
% b0,b2 in R^n,
% c2 in R,
% A in R^{mxn},
% b in R^m.
% ----------------------------------------------------------
% n is the dimension, i.e., the variable x is in R^n,
% m is the number of linear constraints of the form ai'*x<=bi.
% ----------------------------------------------------------
% Optional input - determining the two input parameters N and gap:
% N - number of nodes allowed to develop (the algorithm will stop after
%     evaluating N-1 nodes if N is even, or N if odd). default: inf.
% gap - the difference allowed between U (upper bound on the optimal value)
%       and Lb (lower bound) in the final output solution. default: 0.
% The algorithm will stop when the best solution so far reaches 
% a value U with U-Lb<gap, or when N nodes are evaluated. 
% Otherwise it will stop when no open nodes remain.
% ----------------------------------------------------------
% To run Algorithm BB with an ordering heuristics (H):
% 1. Load the input into the workspace
% 2. Run the command 
%       >> [U,cur_sol,nn,Lb,status]=BB_QBL_Heur(A,b,Q0,b0,b2,c2,N,gap)
% Omitting the last two inputs will determine default values N=inf, gap=0.
%       >> [U,cur_sol,nn,Lb,status]=BB_QBL_Heur(A,b,Q0,b0,b2,c2,N,gap)
%       >> [U,cur_sol,nn,Lb,status]=BB_QBL_Heur(A,b,Q0,b0,b2,c2,[],gap)
%       >> [U,cur_sol,nn,Lb,status]=BB_QBL_Heur(A,b,Q0,b0,b2,c2,N)
% ----------------------------------------------------------
% Output:
% U - optimal value (or the best found upper bound if not fully solved)
% cur_sol - optimal solution (or best solution found so far)
% nn - number of nodes evaluated
% Lb - the best lower bound found (min of L values on the deepest layer of
% nodes opened so far).
% status -  
% '0'   == 'solved to optimum'
% 'inf' == 'infeasible' 
% '1'   == 'solved up to the given gap'
% '2'   == 'solved up to a larger gap (due to reaching N developed nodes)'
% '3'   == 'did not find a feasible solution till reached N nodes'. 
% ----------------------------------------------------------

if nargin==6
    N=inf;
    gap=0;
elseif nargin==7
    gap=0;
end
if min(size(N))<1
    N=inf;
end

U=inf;%Upper bound
n=length(b0);
n2=length(b2);
m=size(A,1);
nq=size(Q0,1);
if ~(size (Q0,2)==nq)
    error('The matrix Q0 is not square.')
end
if min(size(A))>0
    if ~(n==nq&&n==n2&&size(A,2)==n&&length(b)==m)
        error('Some dimensions of input vectors or matrices do not match.')
    end
    elseif ~(length(b)==m)
    error('Some dimensions of input vectors or matrices do not match.')
end
if max(abs(Q0-Q0'))>1e-6
    fprintf('Q0 was changed to be symmetric\n')
    Q0=0.5*(Q0+Q0');
end
if max([size(b0,2),size(b,2),size(b2,2)])>1
    error('b,b0,b2 must be column vectors.')
end
if ~(size(c2,1)==1&&size(c2,2)==1)
    error('c2 needs to be a scalar.')
end
if N<1||gap<0
    error('N should be at least 1, and the allowed gap must be nonnegative.')
end
cur_sol=[];
cands=cell(1,1);
vcands=cell(1,1);
scands=zeros(1,1);
L=inf*ones(1,1);%Lower bounds on nodes
status=inf;
Lb=inf;
%Solving the root node [0 0] ----------------------------------------------------
%node=[0 0]
nn=1;
[cands,vcands,scands,U,cur_sol]=Solve_Root_Node(A,b,Q0,b0,b2,c2,cands,vcands,scands,U,cur_sol);
L(1)=min(vcands{1});
if ~(L(1)<U)
    Lb=L(1);
    if U<inf
        status=0;
    end
    fprintf('The problem has been solved. \nOptimal value = %f\n',U);
    fprintf('The optimal solution is saved in cur_sol\n')
    return;
end
j=0;
%-----------------------------------------------------------
%The nodes [i 2k] and [i 2k+1] for k=0,1,2,...------------------------------

for i=1:m
    %--------------
    %Applying the heuristics: finds an index of a constraint
    %with maximal violating candidates and swaps with the upcoming one.
    ind=Max_of_vio(A,b,i-1,cands);
    A([i,ind],:)=A([ind,i],:);
    b([i,ind])=b([ind,i]);
    %--------------
    cands_new=cell(1,2*nnz(L<U));
    vcands_new=cell(1,2*nnz(L<U));
    scands_new=zeros(1,2*nnz(L<U));
    L_new=inf*ones(1,2*nnz(L<U));%Lower bounds on nodes
    nodes_new=sort([2*j,2*j+1]); %"Child-nodes"
            %Notice that j is a vector whose indices are 
            %the open nodes remaining from the previous layer.
    for k=j %0:2^(i-1)-1 
        if L(j==k)<U 
            if nn>N-2
               if U<inf&&U-Lb<gap
                   status=1;
               elseif U<inf
                   status=2;
               else
                   status=3;
               end
               break;
            end
            node=[i 2*k];
% Optional to display output
%             U0=U;% tic;
            [cands_new,vcands_new,scands_new,U,cur_sol]=...
                Solve_i_in_E_Node(A,b,Q0,b0,b2,c2,cands_new,vcands_new,scands_new,U,cur_sol,node,nodes_new);
%            if U<U0
%                node
%                U
%                 cur_sol
%                 max(A*cur_sol-b)
%            end% toc;
            if min(size(vcands_new{nodes_new==2*k}))>0
                L_new(nodes_new==2*k)=min(vcands_new{nodes_new==2*k});
            end
            node=[i 2*k+1];
           %U0=U; % tic;
            [cands_new,vcands_new,scands_new,U,cur_sol]=Solve_i_not_in_E_Node(A,b,Q0,b0,b2,c2,cands,vcands,scands,cands_new,vcands_new,scands_new,U,cur_sol,node,nodes_new,j);
%            if U<U0
%                U
%            end % toc;
            nn=nn+2;
            if min(size(vcands_new{nodes_new==2*k+1}))>0
                L_new(nodes_new==2*k+1)=min(vcands_new{nodes_new==2*k+1});
            end
        end
    end
    if status<inf&&status>0
        break;
    end
    cands=cands_new(L_new<U);
    vcands=vcands_new(L_new<U);
    scands=scands_new(L_new<U);    
    j=nodes_new(L_new<U);%any node for which L_new>=U is being fathomed.
    L=L_new(L_new<U);
    Lb=min(L_new);
    if size(j,2)<1
       if U<inf
            status=0;
       end
       break;
    elseif U-Lb<gap
       status=1;
       break;
    end
end
if status==0
    fprintf('The problem has been solved. \nOptimal value = %f\n',U);
    fprintf('The optimal solution is saved in cur_sol\n')
elseif status==inf
    fprintf('The problem is infeasible.\n')
elseif status==1
    fprintf('Solved within the given gap.\n')
    fprintf('An apporximate solution with value %f is saved in cur_sol\n',U)
elseif status==2
    fprintf('Run till %d developed nodes.\n',nn)
    fprintf('An apporximate solution with value %f is saved in cur_sol\n',U)
elseif status==3
    fprintf('No feasible solution found till %d developed nodes.\n',nn)
end