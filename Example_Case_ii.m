%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Solves an instance of problem (QBL) for Case (ii):
%%%%% minimize       x'*Q*x-2 b0'*x
%%%%% subject to   ||x||^2 <= delta^2
%%%%%                   A*x<=b.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%           n: dimension of x
%%%%%           m_lin: m is number of linear constraint
%%%%%           r_hard_Case: Dimention of optimal solution set        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function   [Q,a,A,b,delta,OUT]=...
            Example_Case_ii(n,m_lin,r_hard_Case,density,delta)

fpl=@(xu)  0.5*(abs(xu)+xu);
if m_lin<2^12
m=2*n;
else
 m=2^n
end
%%
tic
flagfesi=1;
Itrexam=1;
% while   flagfesi>0
    flagmli=0;
    alphamu=1;
    rc=[-10; 3*randn(n-1,1)];
    rc=sort(rc);
    rc(1:r_hard_Case)=min(rc);
    Q=sprandsym(n,density,rc);
    d=min(rc);
    x3=randn(n,1);
    x3=delta*x3/norm(x3);
    a=-(Q-d*speye(n))*x3;
    a=-0.5*a;
    %%  Optimal of TRS
    %%%%%%%%%%%%%%%%%%%% Function of Japon %%%%%%%%%%%%%%%%%%%%%%
    %%    Optimal solution for Japon and Beck
    % %%%%%%%%%%%%%%%%%%%% Function of Japon %%%%%%%%%%%%%%%%%%%%%%
    [x_trsJ,mu_glo,~] = TRSgeph(...
        2*Q,-2*a,speye(n),delta,0,0,0);
    %%%%%%%%%%%%%%%%%%%% Function of Beck %%%%%%%%%%%%%%%%%%%%%%
    Q0=full(Q);
    [x_trsB,mult,lam,hard]=TRS_Glo(2*Q0,2*a,...
        zeros(n,1),0.5*(delta^2));
    x_trs=[x_trsB x_trsJ];
    F_glob_Trs=diag(x_trs'*Q*x_trs-2*a'*x_trs);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%   Find Other Optimal Solution For TRS
    P2=Q+0.5*mu_glo*eye(n);
    w=null(P2);
    options = optimoptions('quadprog');
    options = optimoptions(options,'Display', 'off');
    q = quadprog(eye(n), zeros(n,1),[],[],P2,a,[],[],x_trsJ,options);
    delbar=sqrt(delta^2-q'*q);
    x_r=randn(r_hard_Case,1);
    x_r=delbar*x_r/norm(x_r);
    x_str_feas=q+w*x_r;
    F_optimal=x_str_feas'*Q*x_str_feas-2*a'*x_str_feas;
    %%  construct  linear constraint
    if n<20
        A=example1(n,m);
    elseif m<200
        A=example2(n,m);
    else
        A=example3(n,m);
    end
    stt=A*x_trs(:,1:2);
    for i=1:size(stt(:,1),1)
        b(i,1)=sum(stt(i,:))+0.5;
    end
%%
x_trs=[x_trs x_str_feas];
A*x_trs-b;
ind=find((A*x_trs(:,4)-b)>0);
Ap=A;bp=b;
A(ind,:)=[];
b(ind)=[];
m1=size(A,1);
y=x_trs;
index=[];
 for i=1:3
    ind=[];
ind1=find((A*x_trs(:,i)-b)>0);
if size(ind1,1)==0
    index=[index;i];
end
 end   
 k=1;
S=[];sb=[];
if size(index,1)
for i=1:3
    ind1=[];
ind1=find((A*x_trs(:,i)-b)>0);
if size(ind1,1)>0   
S(k,:)=A(ind1(1),:);
sb(k,1)=b(ind1(1));
k=k+1;
end
end
A=[A(1:m_lin-3,:);S];
b=[b(1:m_lin-3);sb];
end 
%%
m2=size(A,1);
if  m2<m_lin
AA=zeros(2*n,n);
AA=sparse(A);
k=1;
for i=1:n
    AA(k,i)=1;
    AA(k+1,i)=-1;
    k=k+2;
end
bb=0.8*delta*ones(2*n,1);
ind=find((AA*x_trs(:,4)-bb)>0);
sizind=size(ind);
AAp=AA;bbp=bb;
AA(ind,:)=[];
bb(ind)=[];
m3=size(AA,1);
t=m_lin-m2;
A=[A;AA(1:t,:)];
b=[b;bb(1:t,1)];
else
A=A(m2-m_lin+1:m2,:);
b=b(m2-m_lin+1:end);
end

% [A,b,flagcase2,out]=FindRedun_Case(A,b,delta);

% OUT.Msg=out.Msg
OUT.Msg='Problem is Case 2';
OUT.Dimension=n;
OUT.m_linear=size(A,1);
OUT.r_hard_Case=r_hard_Case;
OUT.Cpu_Exam=toc;
OUT.F_Optimal= F_optimal;
OUT;
end






