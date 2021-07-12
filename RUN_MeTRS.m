%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Solves an instance of problem (QBL):
%%%%% minimize       x'*Q*x-2 b0'*x
%%%%% subject to   ||x||^2 <= delta^2
%%%%%                   A*x<=b.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all
% warning off
fpl=@(xu)  0.5*(abs(xu)+xu);
%%    Free Problem of M-eTRS
n=50;
m=10;
Q=sprandsym(n,1);
a=randn(n,1);
delta=2;
x0=randn(n,1);
x0=delta*x0/norm(x0);
A=randn(m,n);
b=A*x0;

%%    Case i:  Problem of M-eTRS
% n=100;
% m_lin=200;
% density=1;
% r_hard_Case=20;
% alp=10;
% [Q,a,A,b,delta,OUT_Example_Casei]=...
%             Example_Case_i(n,m_lin,r_hard_Case,density,alp,1);
%         OUT_Example_Casei
%         Ap=A;
%%   Example :      Case ii:  Problem of M-eTRS 
% n=100;
% m_lin=200;
% density=1;
% r_hard_Case=20;
% delta=2;
% [Q,a,A,b,delta,OUT_Example_Caseii]=...
%             Example_Case_ii(n,m_lin,r_hard_Case,density,delta);
%         OUT_Example_Caseii
%%  
%%      Algorithm 1  (Find redundant of linear constraints)
tic
[A2,b2,flagcase23,OUT_new_AK2]=Find_Redun_Case(A,b,delta);
OUT_new_AK2.CPU=toc
OUT_new_AK2

%%     EBB Algorithm
disp('**************************************')
disp('*****     AF Algorithm Algorithm        ****')
% disp('**************************************')
OUT_AFA=Solve_MeTRS(Q,a,A,b,delta)
% 
% 
%%   BB Algorithm
disp('**************************************')
disp('*****         BB Algorithm        ****')
disp('**************************************')
tic
Q0=full(Q);
A=full(A);
[U1,cur_sol1,nn1,Lb1,status1]=...
    BB_QBL_Heur(A,b,2*Q0,2*a,zeros(n,1),-delta^2 ,5000,1e-3);
timebeck1=toc;
xb1=cur_sol1;
 OUTBeck1.Msg=status1;
 OUTBeck1.Dimention=n;
 OUTBeck1.NumLin=size(A,1); 
 OUTBeck1.CPU=timebeck1;
 OUTBeck1.Node=nn1;
 OUTBeck1.U=U1;
 OUTBeck1.Lb=Lb1;
OUTBeck1
OUT_AFA




 
 
 
 

