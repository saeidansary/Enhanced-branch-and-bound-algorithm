function [glo,loc,mult,val,valoc,mu,mu_glo]=TRS_Loc(Q0,b0,b1,alpha)
%Computes local minimizers of a TRS
%(the global by calling TRS_Glo, the non-global by the procedure with S)
q0=@(x) 0.5*x'*Q0*x-b0'*x;
%q1=@(x) 0.5*x'*x-b1'*x-alpha;
n=length(b0);
mu=[];
if n==1&&Q0(1,1)<0
    % In the 1-dimensional case no need for the second-small eigenvalue of
    % Q0, but the local-non-global-solution might exist
    % the procedure should be simpler
    discr=b1^2+2*alpha;
    if discr<0
        glo=[];
        loc=[];
        mult=1;
        val=[];
        valoc=[];
    else
        x1=b1+sqrt(discr);
        x2=b1-sqrt(discr);
        if q0(x1)<q0(x2)
            glo=x1;
            val=q0(x1);
            mu_glo=(b0-Q0*x1)/(x1-b1);
            loc=x2;
            valoc=q0(x2);
            mu=(b0-Q0*x2)/(x2-b1);
        else
            glo=x2;
            val=q0(x2);
            mu_glo=(b0-Q0*x2)/(x2-b1);
            loc=x1;
            valoc=q0(x1);
            mu=(b0-Q0*x1)/(x1-b1);
        end
        mult=1;
    end
    return;
end
c1=-2*alpha;
DEL2=b1'*b1-c1;
loc=[];
valoc=inf;
Q0=0.5*(Q0+Q0');
[P,lamQ]=eig(Q0);
c=P'*(b0-Q0*b1);
phi=@(t)sum((c.^2)./(diag(lamQ)+t).^2)-DEL2;
phi_pr=@(t) -2*sum((c.^2)./(diag(lamQ)+t).^3);
phi_pr2=@(t) 6*sum((c.^2)./(diag(lamQ)+t).^4);
[glo,mult,mu_glo,hard]=TRS_Glo(Q0,b0,b1,alpha);
if size(glo,2)<1
    val=[];
    return;
end
val=q0(glo(:,1));
if hard>1 %a hard case with mult>1
    fprintf('A hard case of dim>1. There can be other solutions for the subproblem in trs_trans\n')
end
if (mult==1)&&(size(glo,2)==1)&&(n>1)
    L_IPD=-lamQ(1,1);
    if L_IPD<1e-8
        return;
    end
    L_I1=-lamQ(2,2);
    muL=L_I1+1e-4;
    muU=L_IPD-1e-4;
    mu=(muL+muU)/2;
    while 1
        if abs(phi(mu))<1e-4&&phi_pr(mu)>=0
            if mu>-1e-8
                loc=(Q0+mu*eye(n))\(b0+mu*b1);
                valoc=q0(loc);
            end
            break;
        end
        if phi(mu)<=0||phi_pr(mu)<=0
            mu1L=mu;
            mu1U=muU;
        else
            mu1L=muL;
            mu1U=mu;
        end
        muNe=mu-phi(mu)/phi_pr(mu);
        if (phi(mu)>0&&abs(phi_pr(mu))<1e-4)||(phi(mu)>0&&(muNe<mu1L||muNe>mu1U))
           break;             
        end
        %S=mu-(phi_pr(mu)-sqrt(phi_pr(mu)^2-2*phi_pr2(mu)*phi(mu))/phi_pr2(mu));%A second order polynomial approximation
        S=muNe;%alternative option (Standard Newton - usually better)
        %S=(mu1L+mu1U)/2;%alternative option (Bisection)
        if (~isreal(S))||S>L_IPD||S<L_I1
            S=inf;       
        end
        if S>=mu1L+0.01*(mu1U-mu1L)&&S<=mu1L+0.99*(mu1U-mu1L)%Safeguard
            mu1=S;
        else
            mu1=(mu1L+mu1U)/2;
        end
        muL=mu1L;
        muU=mu1U;
        mu=mu1;
    end
end

    