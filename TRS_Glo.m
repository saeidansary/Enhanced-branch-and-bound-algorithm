function [x_trs,mult,lam,hard]=TRS_Glo(Q0,b0,b1,alpha)
%Computes the global solution(s) of a standard TRS
%min q0(x):=0.5*x'*Q0*x-b0'*x
%s.t.q1(x)<=0, where
%q1=@(x) 0.5*x'*x-b1'*x-alpha;
hard=0;
n=length(b0);
c1=-2*alpha;
DEL2=b1'*b1-c1;
[P,lamQ]=eig(Q0);
c=P'*(b0-Q0*b1);
phi=@(t)sum((c.^2)./(diag(lamQ)+t).^2);
%q1((Q0+t*eye(n))\(b0+t*b1))<=0   <===>  phi(t)<=DEL2
%This can be shown by changing variables to homog. q1 into ||z||<=r
f=@(t)0.5*phi(t)-0.5*DEL2;
phi_pr=@(t) -2*sum((c.^2)./(diag(lamQ)+t).^3);
% phi_pr2=@(t) 6*sum((c.^2)./(diag(lamQ)+t).^4);
gamma=@(t)(phi(t))^(-0.5);
gamma_pr=@(t)-0.5*(phi(t))^(-1.5)*phi_pr(t);
% gamma_pr2=@(t)0.75*(phi(t)^(-2.5)*phi_pr(t)^2)-0.5*(phi(t))^(-1.5)*phi_pr2(t);
f_ga=@(t)gamma(t)-DEL2^(-0.5);
x_trs=[];
mult=1;
if -0.5*(b1'*b1)-alpha>0%-0.5*(b1'*b1)+c1/2>0
    lam=-min(diag(lamQ));
    return;
end
lamQ=diag(lamQ);
for i=2:n
if lamQ(i)-lamQ(1)<1e-7
    mult=mult+1;
end
% mult is the multiplicity of the minimal eigenvalue of Q0
end
%f=@(lam) 0.5*norm((Q0+lam*eye(n))\(b0+lam*b1))^2-b1'*((Q0+lam*eye(n))\(b0+lam*b1))-alpha;
[L,p]=chol(Q0,'lower');
if(p==0)
    x_naive=L'\(L\b0);
    if(0.5*norm(x_naive)^2-b1'*x_naive<=alpha)
        x_trs=x_naive;
        lam=0;
    else
        u=1;
        while(phi(u)>DEL2)
            u=2*u;
        end
        if u>1
            u=u/2;
        end
        if abs(f(u))<1e-8
            lam=u;
        else
            lam=Newton_on_IPD(f_ga,gamma_pr,u,-min(lamQ),1e-7);
        end
        x_trs=(Q0+lam*eye(n))\(b0+lam*b1);
    end
else
% if abs(c(1))>1e-7
%     %c(1)
%     tolep=abs(c(1))/sqrt(DEL2);
% else
%     tolep=(sum((c(abs(c)>1e-7).^2)./((lamQ(abs(c)>1e-7))-min(lamQ)).^2)-DEL2)/(2*sum((c(abs(c)>1e-7).^2)./((lamQ(abs(c)>1e-7))-min(lamQ)).^3));
% end
u=-min(lamQ)+1e-4;%max(tolep,1e-4);
% if tolep<1e-4
%     fprintf('nearly hard case, %1.14f\n',tolep)
% end
if(f(u)<0)%the hard case, where b0-min(eig(Q0))*b1 belongs to
    %Range(Q0-min(eig(Q0))*eye(n))
    hard=mult;
    lam=-min(lamQ);
    [P,D]=eig(Q0);
    Ds=D+lam*eye(n);
    g=P'*b0+lam*P'*b1;
    if n<1+mult
        y1=zeros(mult,1);
    else
        y1=[zeros(mult,1);g(1+mult:n)./diag(Ds(1+mult:n,1+mult:n))];
    end
    x1=P*y1;
    if mult==1||0.5*(x1'*x1)-b1'*x1+c1/2<=0
      t=-(x1-b1)'*P(:,1)+sqrt(((x1-b1)'*P(:,1))^2-(x1-b1)'*(x1-b1)+b1'*b1+2*alpha)*[1,-1];
      x_trs=x1*[1,1]+P(:,1)*t;
    else
        rti=sqrt(2*alpha+b1'*b1-(x1-b1)'*(x1-b1)+((x1-b1)'*P(:,1:mult)*P(:,1:mult)'*(x1-b1)));
        w=-P(:,1:mult)'*(x1-b1)+rti*[1;zeros(mult-1,1)];
        x_trs=x1+P(:,1:mult)*w;
    end
    %This gives all the solutions if mult==1, but if mult>1 there can be
    %infinitely many, as [1;zeros(mult-1,1)] can be replaced with any unit
    %vector in R^mult
else
    while(f(u)>0)
        u=2*u;
    end
    if u>-min(lamQ)+1e-7
        u=u/2;
    end
        if abs(f(u))<1e-8
            lam=u;
        else
            lam=Newton_on_IPD(f_ga,gamma_pr,u,-min(lamQ),1e-7);
        end
    x_trs=(Q0+lam*eye(n))\(b0+lam*b1);
end
end
if max(abs(imag(x_trs)))>1e-4
    fprintf('inaccurate\n');
end
x_trs=real(x_trs);
