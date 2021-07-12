function [A,b,flagcase2,OUT]=Find_Redun_Case(A,b,delta)
%%%%%%%%%%%%%%%%%%%% Algorithm 1 %%%%%%%%%%%%%%%%%%%%
%% This algorithm finds some redundant linear constraints of m-eTRS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   ||x|| <=delta
%%%%%%%%%%%    Ax <= b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Case i :  flagcase2=0
%%%%%%%%%%  Case ii:  flagcase2=1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n]=size(A);
flagfeasible=1;
flagcase2=0;
redun1=[];
Xtild=zeros(n,m);nrmXtild=zeros(m,1);
%%%%%%%%%%%%% distance of center ball %%%%%%%%%%%%%%%%%
for i=1:m
  nrmbi=norm(A(i,:));
  Xtild(:,i)=b(i)*A(i,:)'/nrmbi^2;
  nrmXtild(i)=norm(Xtild(:,i));
  if  nrmXtild(i)>delta
      if b(i) <0
          disp('problem is infeasible')
          flagfeasible=0;
          OUT.feasible=0;
          break
      end
   redun1=[redun1 i];
  end
end
A(redun1,:)=[];b(redun1)=[];Xtild(:,redun1)=[];
m1=size(A,1);
if flagfeasible==1
OUT.feasible=1;
%%%%%%% determining parallel linear Constraint %%%%%%%%% 
e=0.1*exp(1);
P1=1:m1;P2=P1;
  for i=P1
      P2(P2<=i)=[];
     for j=P2
        % par=A(i,:)-(max(abs(A(i,:)))/max(abs(A(j,:))))*A(j,:);
        par=(A(i,:)+e)./(A(j,:)+e);
        if max(par)-min(par)<1e-6
            tji=A(j,:)*Xtild(:,i)-b(j);
            tij=A(i,:)*Xtild(:,j)-b(i);
            if   tji>=0 && tij<=0
                P2(P2==j)=[];
                P1(P1==j)=[];
            elseif    tji<=0 && tij>=0
                P1(P1==i)=[];
                break
            end
        end
    end
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 OUT.num_Parallel_Constraint=m1-length(P1);
m1=size(A,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
P=A*A';
if m1>0
%%%%%%%%%%%%% Chech intersection inside ball %%%%%%%%%%%%%%%%%    
    delta2=delta-1e-7;
    P1=1:m1;P2=P1;
    for i=P1
      P2(P2<=i)=[];
     for j=P2
            tt= P(j,j)*P(i,i)-P(i,j)^2;
            teta1=(b(i)*P(j,j)-b(j)*P(i,j));
            teta2=(b(j)*P(i,i)-b(i)*P(i,j));
            %  nrmXhat2=((teta1^2)*P(i,i)...
            %     +(teta2^2)*P(j,j)+2*teta1*teta2*P(i,j))/tt.^2;
            %  ttt=(teta1*A(i,:)+teta2*A(j,:))/tt;
            % nrmXhat2=norm(ttt);
            nrmXhat2=norm(teta1*A(i,:)+teta2*A(j,:))/tt;
            if nrmXhat2<delta2
                flagcase2=1;
            else
       %%%%%%%%%%%%% Chech Step 5 %%%%%%%%%%%%%%%%%
                tji=A(j,:)*Xtild(:,i)-b(j);
                tij=A(i,:)*Xtild(:,j)-b(i);
            if   tji>=0 && tij<=0
                P2(P2==j)=[];
                P1(P1==j)=[];
            elseif    tji<=0 && tij>=0
                P1(P1==i)=[];
                break
            end
            end
        end
    end
OUT.Step3_redun=m1-length(P1);  
A=A(P1,:);b=b(P1);Xtild=Xtild(:,P1);
m1=size(A,1);
if flagcase2==1
  OUT.Msg='Problem is Case 2';  
else
    m1=size(A,1);
    P1=1:m1;P2=P1;
    for i=P1
      P2(P2<=i)=[];
     for j=P2
     if A(j,:)*Xtild(:,i)-b(j)>=0
        P2(P2==j)=[];
        P1(P1==j)=[];
         break
     end
  end
  end
   OUT.Msg='Problem is Case 1';
   A=A(P1,:);
   b=b(P1);   
end
end%%% m1>0
m3=size(A,1);
OUT.Linear_Constraints=m3;
OUT.Number_noparalel_redun=m1-length(P1);
end%%%flagfeasible


