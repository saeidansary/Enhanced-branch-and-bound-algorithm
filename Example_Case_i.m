%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Solves an instance of problem (QBL) for Case (i):
%%%%% minimize       x'*Q*x-2 b0'*x
%%%%% subject to   ||x||^2 <= delta^2
%%%%%                   A*x<=b.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%           n: dimension of x
%%%%%           m_lin: m is number of linear constraint
%%%%%           b=alp*onese(m_lin,1)
%%%%%           delta <= alp*sqrt(2)     
%%%%%           r_hard_Case: Dimention of optimal solution set 
%%%%%           flagtest=1  :  Test Class 1  of Case i
%%%%%           flagtest=ecept 1  :  Test Class 2  of Case i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function   [Q,a,A,b,delta,OUT]=...
            Example_Case_i(n,m_lin,r_hard_Case,density,alp,flagtest)
fpl=@(xu)  0.5*(abs(xu)+xu);
%%  Example :      ||x||_{inf} <= alp 
if m_lin<n/2
m=2*n;
else
 m=m_lin+n/2;
end
flagcase1=0;
tic
while(flagcase1==0)
k=1;
A=zeros(2*m,n);
A=sparse(A);
for i=1:n
    A(k,i)=1;
    A(k+1,i)=-1;
    k=k+2;
end
A;
b=alp*ones(2*m,1);
delta=0.8*alp*sqrt(2);
t=delta;
m=2*m;
% fprintf(' Elapsed time of Example is %2.2f seconds.   \n',toc )
%%  Object Function (Hard case 2)
rc=[-1; 3*randn(n-1,1)];
rc=sort(rc);
rc(1:r_hard_Case)=min(rc);
Q=sprandsym(n,density,rc);
d=min(rc);
x3=randn(n,1);
x3=delta*x3/norm(x3);
a=-(Q-d*speye(n))*x3;
a=-0.5*a;
%%    Optimal solution for Japon and Beck
% %%%%%%%%%%%%%%%%%%%% Function of Japon %%%%%%%%%%%%%%%%%%%%%%
  [x_trsJ,mu_glo,~] = TRSgeph(...
                        2*Q,-2*a,speye(n),t,0,0,0);        
%%%%%%%%%%%%%%%%%%%% Function of Beck %%%%%%%%%%%%%%%%%%%%%%  
 Q0=full(Q);
[x_trsB,mult,lam,hard]=TRS_Glo(2*Q0,2*a,...
                                  zeros(n,1),0.5*(t^2)); 
x_trs=[x_trsB x_trsJ];                               
F_glob_Trs=diag(x_trs'*Q*x_trs-2*a'*x_trs);                             
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Find Other Optimal Solution For TRS
P2=Q+0.5*mu_glo*eye(n);
w=null(P2);
options = optimoptions('quadprog');
options = optimoptions(options,'Display', 'off');
% q = quadprog(eye(n), -2*cc,[],[],P2,a,[],[],y,options);
q = quadprog(eye(n), zeros(n,1),[],[],P2,a,[],[],x_trsJ,options);
delbar=sqrt(t^2-q'*q);
x_r=randn(r_hard_Case,1);
x_r=delbar*x_r/norm(x_r);
x_str_feas=q+w*x_r; 
x_str_feas'*Q*x_str_feas-2*a'*x_str_feas;
%% 
x_trs=[x_trs x_str_feas];
ind=find((A*x_trs(:,4)-b)>0);
A(ind,:)=[];
b(ind)=[];
[m1,~]=size(A);
y=x_trs;
mm=3;
x_trs1=[];
for i=1:3
    ind=[];
    ind1=find((A*x_trs(:,i)-b)>0);
    if size(ind1,1)>0
        mm=mm-1;
    else
        x_trs1=[ x_trs1  x_trs(:,i)];
    end
end
if flagtest==1
    %%%%%%%%%%%%%%   infeasible optimal Trs_glo  %%%%%%%%%%%%%%
if mm>0
    %%%%%%%%%%%   2 x_trs * x <= 2 ||x_trs||^2- \epsilon
    P=[];px=[];
    for i=1:3
        P(i,:)=2*x_trs1(:,i)';
        px(i)=2*norm(x_trs1(:,i))^2-0.05;
    end
    b=[b;px'];
    A=[A;P];
end
else
%%%%%%%%%%%%%%   Parallel linear Constraints  %%%%%%%%%%%%%% 
I=eye(n);
if m_lin<100
    iterLin=floor(m_lin/3.5)
else
    iterLin=floor(m_lin/4)
end
ind=[];
for i=1:3
    ind1=find(x_trs1(:,i)>0);
    if isempty(ind1)==1
        ind(i)=i;
    else
        ind(i)=ind1(1);
    end
end
for j=1:iterLin
    for i=1:3
        P(i,:)=(2)*x_trs1(:,i)'+(1e-3)*I(:,ind(i))'*j;
        px(i)=2*norm(x_trs1(:,i))^2-(1e-5)*j ;
    end
    b=[b;px'];
    A=[A;P];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=size(A,1);
A=A(m-m_lin+1:end,:);
b=b(m-m_lin+1:end);
% [A,b,flagcase2,out]=FindRedun_Case(A,b,delta);
if n<1000 && m_lin<700
[A25,b25,flagcase225,out]=Find_Redun_Case(A,b,delta);
if sum(out.Msg)==1466
    flagcase1=1;
end
else
flagcase1=1;
end
end
OUT.dimension=size(A,2);
OUT.m_linear=size(A,1);
OUT.r_hard_Case=r_hard_Case;
% OUT.Msg=out.Msg
OUT.Msg='Problem is Case 1';
OUT.Cpe_Exam=toc;
OUT.F=F_glob_Trs(1);

%%  Figure 3
if n==3
%     A=Ap;
%     b=bp;
for i=1:m1+mm
Hyperplane(A(i,:),-2*delta,2*delta,-2*delta,2*delta,b(i),1)
hold on
end

[x,y,z] = ellipsoid(0,0,0,delta,delta,delta,200);
surf(x,y,z)

for i=1:3
 plot3(x_trs(1,i),x_trs(2,i),x_trs(3,i),'r*')   
 hold on
end

axis equal
end



