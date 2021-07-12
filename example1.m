 function A=example1(n,m)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%  ||x||_1 <=1
 %%%%%%   A*x <=1
 %%%%%%   A \in R^{m,n}
 %%%%%%   m<= 2^n
 %%%%%%  for example should n<22 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clc
% clear all
% n=3;
% m=3;
% s=1;

A=zeros(m,n);
tic
for t=1:n
p=example_find(n,m,n-t+1);
A(:,t)=p;
end

function p=example_find(n,m,s)
    k=1;
for i=1:2^(n+1)/2^s
    if mod(i,2)==1  
  for j=1:2^(s-1) 
      if k<=m
    p(k,1)=1;
    k=k+1;
      end
  end
    end
   if mod(i,2)==0  
     for j=1:2^(s-1) 
         if k<=m
    p(k,1)=-1;
    k=k+1;
         end
     end
    end

end
p;

end
 end











