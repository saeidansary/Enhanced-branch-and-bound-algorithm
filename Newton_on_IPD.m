function mu=Newton_on_IPD(f,fpr,lam1,lo,eps)
%Performs a safeguarded Newton algorithm to compute mu^* in IPD
mu=lam1;
while f(mu)<0
    mu=mu+1;
end
if mu>lam1
    mu=mu-1;
end
muU=mu+1;
while abs(f(mu))>eps
    mu_New=mu-f(mu)/fpr(mu);
    if (mu_New>lo)&&(mu_New<muU)&&abs(f(mu_New))<0.99*abs(f(mu))
        mu=mu_New;
    else
        mu=(mu+lo)/2;
    end
end
