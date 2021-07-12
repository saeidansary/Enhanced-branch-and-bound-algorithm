function  OUT=Solve_MeTRS(Q,a,A,b,delta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Solves an instance of problem (QBL):
%%%%% minimize       x'*Q*x-2 b0'*x
%%%%% subject to   ||x||^2 <= delta^2
%%%%%                   A*x<=b.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fpl=@(xu)  0.5*(abs(xu)+xu);

%% Finde redundant linear constraint and determain Case (i) or (ii)
tic
mp=size(A,1);
[A,b,flagcase2,out]=Find_Redun_Case(A,b,delta);
[m,n]=size(A);
Timecase=toc;
flagfeasible=out.feasible;
 if flagfeasible==1
%% Case i
if flagcase2==0
 [tim,xm,Msg]=solv_metrsnGnL(A,b,Q,a,delta,1); 
OUTDRS.Case='Problem is Case (i)';
OUTDRS.Lineer_Con=mp;
OUTDRS.Redund_Lineer=mp-m;
OUTDRS.CpuCase=Timecase;
OUTDRS.Dimention=n;
OUTDRS.LinearCon=m;
OUTDRS.Msg=Msg;
OUTDRS.CPU=tim;
OUTDRS.x=xm;
fDrs=xm'*Q*xm-2*a'*xm;
OUTDRS.f=fDrs;
OUTDRS.fesDR=norm(xm)^2-delta^2;
OUTDRS.maxlin=norm(fpl(A*xm-b));
OUTDRS;
OUT=OUTDRS;
%% Case ii    
else
 rrr=1;
 Q0=full(Q);
 A=full(A);
 AA1=A;bb1=b;del1=delta;QQ=Q;aa=a;
 save('example','QQ','aa','AA1','bb1','del1')
save('chinexa','rrr');
[U,cur_sol,nn,Lb,status]=...
    BB_QBL_Heur_EBB(A,b,2*Q0,2*a,zeros(n,1),-delta^2 ,1000,1e-8);
timebeck=toc;
xb=cur_sol;
% U,nn,Lb,status
OUTBeck.Case='Problem is Case (ii)';
OUTBeck.Lineer_Con=mp;
OUTBeck.Redund_Lineer=mp-m;
OUTBeck.CpuCase=Timecase;
OUTBeck.Dimention=n;
OUTBeck.LinearCon=m;
OUTBeck.Msg=status;
OUTBeck.Dimention=n;
OUTBeck.NumLin=size(A,1);
OUTBeck.CPU=timebeck;
OUTBeck.Node=nn;
OUTBeck.x=cur_sol(:,1);
if size(cur_sol(:,1),1)>1
OUTBeck.f=cur_sol'*Q*cur_sol-2*a'*cur_sol;
% OUTBeck.fTRS=fglobstr;
OUTBeck.fesbeck=norm(cur_sol)^2-delta^2;
OUTBeck.maxlin=norm(fpl(A*xb-b));
end
 OUT=OUTBeck;
end
%%
 else
  OUT=out; 
 end%%%  flagfeasible==1

