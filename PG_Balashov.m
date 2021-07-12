function     out=PG_Balashov(A,b,delta)
%------------------------------------------------------------
%--------min ||(Ax-b)+||^2--------------------------------
%-------- s.t.  ||x||=delta------------------------------------
%------------------------------------------------------------
%------------------------------------------------------------

%%  Main algorithm
[m,n]=size(A);
pl=@(x) max(x,0);
PB=@(x) delta*x/norm(x);

Fx=@(A,b,x) norm(pl(A*x-b))^2 ;
gradq=@(A,b,x) A'*pl(A*x-b) ;
tolx=1e-8;
tolF=1e-16;
gamma=0.4;
% eta=1.4;
eta=norm(A,2)*norm(A',2)+1;
t=0.2;
 maxitr=5000;
% x0=delta*orth(rand(n,1));
x0=rand(n,1);
if norm(gradq(A,b,x0))>1e-5
x0=delta*gradq(A,b,x0)/norm(gradq(A,b,x0));
else 
    x0=delta*orth(x0);
end

h=1;
ff(h)=Fx(A,b,x0);
for i=1:5000
    
    %% back tracking
    g=gradq(A,b,x0);
    g=g;
%     L=1/t;
%     for j=1:100
%         L=eta*L;
% %         if Fx(A,b,x0)-Fx(A,b,PB(x0-(1/L)*g)) >= gamma*L*norm(x0-PB(x0-(1/L)*g))^2
% %             break
% %         end
%     end
    x00=x0;
    x0=x0-(1/eta)*g; 
    x0=PB(x0);
  h=h+1;
ff(h)=Fx(A,b,x0);  
    if abs(Fx(A,b,x0) - Fx(A,b,x00))<tolF  || norm(x0-x00)<tolx
        break
    end
    
end
out.itr=i;
out.x=x0;
out.Fval=Fx(A,b,x0);
out.feasi=norm(x0)-delta;
out;

 


