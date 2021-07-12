function     [x p]=hardoptfind(y,K,xo,mult,mu_glo,Q,a,cc,t,maxith)
n=size(Q,2);
load example
A=AA1;b=bb1;
% Q=QQ;a=aa;t=del1;
P2=Q+mu_glo*eye(n);
a=a+mu_glo*cc;
w=null(P2);
% options = optimoptions('quadprog');
% options = optimoptions(options,'Display', 'off');
% q = quadprog(eye(n), -2*cc,[],[],P2,a,[],[],y,options);
q = quadprog(eye(n), -2*cc,[],[],P2,a,[],[],y);
delbar=sqrt(t^2-q'*q+2*cc'*q);
delbar2=sqrt(delbar^2+(w'*(q-cc))'*(w'*(q-cc)));
% y1=randn(size(w,2),1);
% y3=delbar2*y1/norm(y1)-w'*(q-cc);
% x2=K*(q+w*y3)+xo;
% x1=K*y(:,1)+xo;
A1=A*K*w;
b1=b-A*K*q-A*xo;
b2=b1+A1*(w'*(q-cc));
% n2=size(A1,2);

p=1;
% x=linprog(zeros(n2,1),A1,b2,[],[],...
%           -2*delbar*ones(n2,1),2*delbar*ones(n2,1));
% if size(x,1)>=1
% if max(A1*x-b2)<=1e-6
for i=1:1
%   [out1]=optim_PLTRS(A1,b2,n2,delbar2);
%   out2=fun_PGbrack(A1,b2,delbar2)
% out= Plus_proj_backtracking(A1,b2,delbar2);
out=PG_Balashov(A1,b2,delbar2);
       yy=out.x-w'*(q-cc);
%     sss(i)=max(A1*yy-b1);
    if max(A1*yy-b1)<=1e-5
%       [out.line  A1*out.x-b2];
        p=max(A1*yy-b1);
        break
    end
end
x=q+w*yy;
% else
%   x=y(:,1);  
% end
% else
%   x=y(:,1);  
% end
% out;
% min(sss);









