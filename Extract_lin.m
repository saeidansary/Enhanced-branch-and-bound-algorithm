function [K,xo]=Extract_lin(A,b,tol)
%Computes K and xo for the affine substitution x=K*y+xo
%where the latter is the general solution (for all y in R^k, k=dim(null(A)))
%n=size(A,2);
m=size(A,1);
r=rank(A,1e-5);
rb=rank([A,b],1e-5);
if r<rb
   xo=[];    
   K=[];
elseif r<m
    R=rref([A,b],1e-5);fprintf('r<m\n')
    Rb=R(1:rb,1:end-1);
    bb=R(1:rb,end);
    xo=Rb\bb;
    if norm(A*xo-b)>tol
        xo=A\b; %possible warning appears here
    end
    K=null(A);
else
    xo=A\b;
    K=null(A);
end

