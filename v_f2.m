function [xq,D,con]=v_f2(Q,c,A,b,xq)% obtaining a local vertex optimal
  
   n=size(Q, 1); 
   xq=cplexlp(Q*xq+c,A,b);
%    xq=linprog(Q*xq+c,A,b);


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
  end
    con=1;
   


