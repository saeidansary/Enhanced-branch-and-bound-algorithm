 function A=example3(n,m)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%  ||x||_1 <=1
 %%%%%%   A*x <=1
 %%%%%%   A \in R^{m,n}
 %%%%%%   m<= 2^n
 %%%%%%  for example should m << 2^n 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=zeros(m,n);
k=1;
while(k<m+1)
for i=1:n    
t=rand;
if t>0.5
   Q(1,i)=2*randn;
else
     Q(1,i)=-1;
end
end
flag=0;
% for j1=1:k-1    
%     for j2=1:k-1
%         if j1<j2 && sum(abs(A(j1,:)-A(j2,:)))==0
%         flag=1;   
%         end
%     end
% end
if flag==0
    A(k,:)=Q; 
    k=k+1;
end
end
A;
 end 













