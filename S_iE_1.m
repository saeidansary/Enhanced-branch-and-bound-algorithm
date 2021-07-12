function [is,js1]=S_iE_1(node,nodes_new)
%  computes the indices of nodes in S_{i,E}^1 that are open
is=node(1);
js1=[];
bx=de2bi(node(2));
s=length(bx);
for i=s:-1:1
    if bx(i)==1
        ch=bx;
        ch(i)=0;
        n1=bi2de(ch);
       if min(abs(nodes_new-n1))==0
            %checks whether n1 belongs to nodes_new
            js1=[js1,n1];
        end
    end
end
    
    