function [xq,D,con]=v_f(Q,c,A,b,xq)% obtaining a local vertex optimal
  
   n=size(Q, 1); 
   xq=cplexlp(Q*xq+c,A,b);
%    xq=linprog(Q*xq+c,A,b)
 while(1)

  D=zeros(n);    % directiones
  idx= A*xq-b>=-10^(-7);
  JJ=A(idx,:); % active constraints
  
  if size(JJ,1)~=n% degenerated case or less active ...
       con=0;
       return; 
  end

  for i=1:n
      GH=JJ;
      GH(i,:)=[];
      [~, ~, Vy] = svd(GH);
       d = Vy(:,end);
     if JJ(i,:)*d>0
         d=-d;
     end
     d=norm(d)^(-1)*d;
     D(:, i)=d; 
     
     idx1= A*d>=10^(-7);
     tt=(b(idx1)-A(idx1,:)*xq)./(A(idx1,:)*d);
     ax=xq-max(-tt)*d;
     
      if ax'*Q*ax+2*ax'*c<xq'*Q*xq+2*xq'*c-10^(-5)
          xq=ax;
          con=0;
          break;
      end
  end
  if i==n
    con=1;
    return;
  end

end

