function [xlocal,laglocal,fvallocal,fla,kkt1,kkt2]=localminimizer(A,a,delta,n)
% This function computes the candidat for the local non-global minimizer of the follwoing TRS:
%   min x'Ax+2*a'x
%       x'x<=delta
%=========================================================================
opts.v0=full(sum(A))';
opts.isreal=1;
% opts.issym=1;
A=0.5*(A+A');
if size(A,1)<=1000
  lam124=sort(eigs(A)); 
 lambdamin=[lam124(1)  0;0  lam124(2)];
else
[~,lambdamin]=eigs(A,2,'sa',opts);
end
lambdamin2=lambdamin(2,2);
M1=[sparse(n,n) speye(n);speye(n) sparse(n,n)];
opts.v0=randn(2*n,1);
fla=0; % fla will changed to one if there is no LNGM for TRS
%===================================================
[v,d]=eigs(@M0,2*n,-M1,2,'lr',opts);
if (abs(imag(d(2,2))))>1e-8
%     fprintf('Second eigenvalue of M(lambda) is complex \n')
    if real(d(2,2))<=-lambdamin2
        fla=1;
%         fprintf('fla=1: there is no local nonglobal minimizer\n')
    else
        [~,d,flagsecon]=realM0n1(v(:,1),d(1,1),lambdamin2);
        if flagsecon==0,
            fla=1;
%             fprintf('fla=1: there is no local nonglonal minimizer \n');
        elseif (d>=0)&&((d+lambdamin2)>=1e-8)&&((d+lambdamin1)<-1e-8)
            df=d;
            opts.v0=randn(2*n,1);
            opts.issym=0;
            opts.isreal=1;
            %========= added for dense problme
            %             opts.tol=1e-5;
            % %             opts.maxit=5;
            %             opts.p=5;
            %================
            [v,d]=eigs(@Asigma,2*n,1,df,opts);
            xopt1=-sign(a'*v(n+1:end))*sqrt(delta)*v(1:n)/norm(v(1:n));
            fval1=xopt1'*A*xopt1+2*a'*xopt1;
%             fprintf('eigenvector for the second real eige of M(lambda) is computed\n')
            lambda_1LNGM=d;
            kkt1=norm((A+lambda_1LNGM*speye(n))*xopt1+a,inf);
            kkt2=lambda_1LNGM*(norm(xopt1)^2-delta);
%             fprintf('kkt1 and kkt2 for LNGM are: %g  %g \n',kkt1,kkt2);
        else
            fla=1;
%             fprintf('fla=1: there is no local nonglonal minimizer \n');
            
        end
    end
elseif (d(2,2)>1e-12)&&((d(2,2)+lambdamin2)>=1e-8)
    xopt1=-sign(a'*v(n+1:end,2))*sqrt(delta)*v(1:n,2)/norm(v(1:n,2));
    fval1=xopt1'*A*xopt1+2*a'*xopt1;
    lambda_1LNGM=d(end,end);
    kkt1=norm((A+lambda_1LNGM*speye(n))*xopt1+a,inf);
    kkt2=lambda_1LNGM*(norm(xopt1)^2-delta);
%     fprintf('kkt1 and kkt2 for LNGM are: %g  %g \n',kkt1,kkt2);
elseif (abs(d(2,2))<=1e-12)&&((d(2,2)+lambdamin2)>=1e-8)
    xopt1=-sign(a'*v(n+1:end,2))*sqrt(delta)*v(1:n,2)/norm(v(1:n,2));
    fval1=xopt1'*A*xopt1+2*a'*xopt1;
    lambda_1LNGM=d(end,end);
    kkt1=norm((A+lambda_1LNGM*speye(n))*xopt1+a,inf);
    kkt2=lambda_1LNGM*(norm(xopt1)^2-delta);
%     fprintf('kkt1 and kkt2 for LNGM are: %g  %g \n',kkt1,kkt2);
else
%     fprintf('fla=1: there is no  local minimzer\n')
    fla=1;
end
if fla==1
    xlocal=[];laglocal=[];fvallocal=[];kkt1=[];kkt2=[];
else
    xlocal=xopt1; laglocal=lambda_1LNGM; fvallocal=fval1;
end
%==============================================
    function f=M0(x)
        f=[-x(1:n)+A*x(n+1:end);A*x(1:n)-a*(a'*x(n+1:end))/delta];
    end
%==============================================
    function y=Asigma(x)
        y = ([-A  (a*a')/delta;speye(n) -A]-df*speye(2*n))\x;
    end

end
% function [v,d,flagsecon]=realM0n1(v,d,lambdamin2,A,n,a,delta)
% %===============================================================
% alpha1=[];base1=[];
% flag_shd=1;
% base1=[base1,v];
% alpha1=[alpha1;d/(v'*v)];
% flagsecon=1;
% while flag_shd
%     opts.v0=randn(2*n,1);
%     opts.isreal=0;
%     opts.issym=0;
%     [v,d]=eigs(@Areal,2*n,1,'lr',opts);
%     if (real(d))<=-lambdamin2
%         flagsecon=0;
%         break;
%     elseif (abs(imag(d)))<=1e-8
%         break;
%     end
%     alpha1=[alpha1;d/(v'*v)];
%     base1=[base1,v];
%     
% end
% %============================================
%     function [v] = Areal(x)
%         if flag_shd==1
%             v=[-A*x(1:n)+a*(a'*x(n+1:end))/delta;x(1:n)-A*x(n+1:end)]-(base1*(alpha1.*(base1'*x)));
%         end
%     end
% end

