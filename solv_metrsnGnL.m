function  [tim,xm,Msg]=solv_metrsnGnL(A,b,Q0,a,t,flagsolv)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Solves an instance of problem (QBL):
%%%%% minimize       x'*Q*x-2 b0'*x
%%%%% subject to   ||x||^2 <= delta^2
%%%%%                   A*x<=b.
%%%%% A fast algorithm for solving the extended trust region
%%%%% subproblem with non-intersecting linear constraints
%%%%%   if  flagsolv=1, then for solving we use  Generalized eigenvalue
%%%%%   if  flagsolv>1, then for solving we use
%%%%%                   1-eTRS, TRS and LNG_TRS of Beck Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n]=size(A);
tic
%%%%%%%%%%  Find global solution of TRS
if   flagsolv==1
    [x_trs,mu_glo,hard] = TRSgeph(...
                        2*Q0,-2*a,speye(n),t,0,0,1);
    mult=hard;
    if hard==1

        P2=Q0+0.5*mu_glo*eye(n);
        w=null(P2);
        options = optimoptions('quadprog');
        options = optimoptions(options,'Display', 'off');
    q = quadprog(eye(n),zeros(n,1),[],[],P2,a,[],[],[],options);
    delbar=sqrt(t^2-q'*q);
    x_trs(:,1)=q+delbar*w;
    x_trs(:,2)=q-delbar*w;
    end
else
    Q0=full(Q0);
    [x_trs,mult,mu_glo,hard]=TRS_Glo(2*Q0,2*a,zeros(n,1),0.5*t^2);
end
f_trs=x_trs'*Q0*x_trs-2*a'*x_trs;
%%
flag=0;
m2=size(x_trs,2);
for i=1:m2
    if   sum(A*x_trs(:,i)-b<=1e-8) == m  %% global of trs is global of m-etrs
        xm=x_trs(:,i);
        flag=1;
        Msg='Optimal solution of TRS is optimal.';
        break
    end
end
%%%%%%%%  If global of TRS is not global for m-eTRS
if flag==0
    xm=solv_metrsnGnL_DrS(x_trs,A,b,Q0,a,t,mult,mu_glo,hard,flagsolv);
    f_xm=xm'*Q0*xm-2*a'*xm;
    if   abs(f_trs-f_xm)<=1e-6
        Msg='Optimal solution of TRS is optimal.';
    else
        Msg='Optimal solution of TRS is not optimal.';
    end
end
tim=toc;
end
    function  xm=solv_metrsnGnL_DrS(x_trs,A,b,Q0,a,t...
                                         ,mult,mu_glo,hard,flagsolv)
        [m,n]=size(A);
         m3=1;
         if hard==1
             m3=size(x_trs,2);
             for ii=1:m3
                 [~,inmaxG(ii,1)]=max(A*x_trs(:,ii)-b);
             end
         else
             [~,inmaxG]=max(A*x_trs(:,1)-b);
         end
        %eigs(sparse(Q0),1,'sa')<=-1e-8
        %         [~,x_loc,~,~,valoc,~,~]=...
        %             TRS_Loc(2*Q0,2*a,zeros(n,1),0.5*t^2);
        %%%%%  Find local nonglobal of TRS
        flageigeq=0;%%% if flageigeq=0,then hav not  local nonglobal
        x_loc=[];
        [v ei]=eigs(Q0,2,'sa');
        ei=diag(ei);
        if abs(ei(1)-ei(2))<1e-7 || abs(v(:,1)'*a)<=1e-7
            flageigeq=1;
        end
        if flageigeq==0
        if   flagsolv==1 
            [x_loc,laglocal,valoc,fla,kkt1,kkt2]=localminimizer(Q0,...
                                                            -a,t^2,n);
        else
            Q0=full(Q0);
            [x_loc,valoc,mu]=TRS_Loc144(2*Q0,2*a,zeros(n,1),0.5*t^2....
                ,x_trs,mult,mu_glo,hard);
        end
        end
        %%%%%%    Step 3
        if  size(x_loc,1)>0  &&  flageigeq==0
            [~,inmaxL]=max(A*x_loc-b);
            if sum(A*x_loc-b<=1e-8) == m
                %             [U,cur_sol,nn,Lb,status]=...
                %                 BB_QBL_Heur(A(inmaxG,:),b(inmaxG)...
                %                 ,2*Q0,2*a,zeros(n,1),-t^2);%,N,gap);
                for ii=1:m3
         [cur_sol1(:,ii),f3(ii,1)]=eqTRS(Q0,a,...
                        A(inmaxG(ii,1),:)',b(inmaxG(ii,1)),t,flagsolv);              
                end
           [f,indm]=min(f3);
            xm=cur_sol1(:,indm);     
                if valoc<f
                    xm=x_loc;
                else
                    xm=cur_sol1(:,1);
                end
         %%%%%%    Step 4
%             elseif inmaxL==inmaxG  
%                 %             [U,cur_sol,nn,Lb,status]=...
%                 %                 BB_QBL_Heur(A(inmaxG,:),b(inmaxG),...
%                 %                 2*Q0,2*a,zeros(n,1),-t^2);%,N,gap);
%                 [cur_sol,~]=eqTRS(Q0,a,A(inmaxG,:)',b(inmaxG),t,flagsolv);
%                 xm=cur_sol(:,1);
            else
                %             A1=A([inmaxG,inmaxL],:);
                %             b1=b([inmaxG,inmaxL]);
                %             [U,cur_sol,nn,Lb,status]=...
                %                 BB_QBL_Heur(A1,b1,...
                %                 2*Q0,2*a,zeros(n,1),-t^2);%,N,gap);
                for ii=1:m3
         [cur_sol1(:,ii),f3(ii)]=eqTRS(Q0,a,A(inmaxG(ii),:)'...
                                         ,b(inmaxG(ii)),t,flagsolv);
                end
           [f,indm]=min(f3);
            cur_sol1=cur_sol1(:,indm);
                [cur_sol2,f1]=eqTRS(Q0,a,A(inmaxL,:)',b(inmaxL),t,flagsolv);
                if f<f1
                    xm=cur_sol1(:,1);
                else
                    xm=cur_sol2(:,1);
                end
            end
        else
                for ii=1:m3
         [cur_sol1(:,ii),f3(ii,1)]=eqTRS(Q0,a,...
                        A(inmaxG(ii,1),:)',b(inmaxG(ii,1)),t,flagsolv);              
                end
           [f,indm]=min(f3);
            xm=cur_sol1(:,indm);
       end
        
    end
 