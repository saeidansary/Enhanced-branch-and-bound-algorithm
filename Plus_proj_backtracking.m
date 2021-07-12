function     out= Plus_proj_backtracking(A,b,delta)
%------------------------------------------------------------
%--------min ||(Ax-b)+||^2--------------------------------
%-------- s.t.  ||x||=delta------------------------------------
%------------------------------------------------------------
%------------------------------------------------------------
n=size(A,2);
plus=@(x)(max(0,x));
func=@(x)(norm(plus(A*x-b))^2);
gradF=@(x)(2*A'*plus(A*x-b));
P=@(x)delta*(x/norm(x));
%------------------------------------------------------------
eps=1e-16; maxIter=2000;
%------------------------------------------------------------
%%
x0=rand(n,1);x0=P(x0); iter=0;
nfe=0;
% L1=sort(eig(A'*A));
% L=0.5*(L1(end)+L1(1));
% tic
while iter<maxIter
    iter=iter+1;
    %................................................
    g=gradF(x0);
%     alpha=1/L;
    alpha=1;
    beta=0.5; flg=1;
    while flg==1
        nfe=nfe+1;
        z1=x0-alpha*g;
        z2=x0-P(z1);
        G_alpha=(1/alpha)*z2;
        f1=func(x0-alpha*G_alpha)-func(x0);
        f2=-alpha*g'*G_alpha+(alpha/2)*norm(G_alpha)^2;
        if f1>f2
            alpha=beta*alpha;
        else
            flg=0;
        end
    end
    x1=x0-alpha*G_alpha;
%     clc
%     disp(['iter: ',num2str(iter)])
%     disp(' ')
%     disp(['f(x): ',num2str(func(x1))])
    %     hold on
    %     plot(iter,func(x1),'or')
    %---------------------------------------------------------
    if norm(x1-x0,inf)<=eps
        break
    end
    x0=x1;
end

% time=toc;
time='turn off';
out.itr=iter;
out.nfe=nfe;
out.x=x1;
out.Fval=func(x1);
out.feasi=norm(x1)-delta;
out.CPU=time;

