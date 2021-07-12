function  [x_s,f]=eqTRS(Q,c,a,b,t,flagsolv)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Solves an instance of problem (QBL):
%%%%% minimize       x'*Q*x-2 c'*x
%%%%% subject to   ||x||^2 <= t^2
%%%%%                   a'*x = b.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=size(Q,1);

% [w,xb]=Extract_lin(a',b,1e-8);
% xb
%%
w=sparse(null(a'));
if abs(b)<1e-8
    xb=zeros(n,1);
else
    xb=b*a/(norm(a)^2);
end
%%  TRS Beck

% % % %Computes the global solution(s) of a standard TRS
% % % % [x_trs,mult,lam,hard]=TRS_Glo(Q0,b0,b1,alpha)
% % % %min     0.5*x'*Q0*x-b0'*x
% % % %s.t.    0.5*x'*x-b1'*x<=alpha;
if flagsolv>1
[x_trs1,~,~,~]=TRS_Glo(2*w'*(Q*w),2*w'*(c-Q*xb),...
                                  zeros(n-1,1),0.5*(t^2-xb'*xb));
  y1=x_trs1(:,1);
 x_s = xb + w*y1; 
 f=x_s'*Q*x_s - 2*c'*x_s;
end
%%  Japon TRS
% % % % minimize (x^TAx)/2+ ax
% % % % subject to 
% % % %if  fequal1==0 then   x^TBx <= Del^2
% % % %if  fequal1==1 then   x^TBx == Del^2
 if flagsolv==1    
[x_trs,~,~] = TRSgeph(...
           2*w'*(Q*w),-2*w'*(c-Q*xb),speye(n-1),sqrt(t^2-xb'*xb),0,0,0);
 y=x_trs(:,1);
 x_s = xb + w*y;
 f=x_s'*Q*x_s - 2*c'*x_s;
 end
