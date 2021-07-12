function [x,val,ca,valc,flag_GeSy,erg]=Lin_Reduced_dim_TRS(Q0,b0,b1,c1,Aeq,beq)
%Solves a TRS obtained by extracting variables from the linear system of
%equations Aeq*x=beq and plugging them into the original TRS with
%Q0,b0,b1,c1
x=[];
val=inf;
ca=[];
valc=[];
flag_GeSy=0;
erg=0;
if(size(Aeq,1)>0)
    [K,xo]=Extract_lin(Aeq,beq,1e-4);%x=K*y+xo
    if min(size(xo))<1
        return;
    end     
    if min(size(K))>0
        [y,yloc]=TRS_Loc(K'*Q0*K,K'*(b0-Q0*xo),K'*(b1-xo),-c1/2-0.5*(xo'*xo)+b1'*xo);
    else
        if 0.5*(xo'*xo)-b1'*xo+c1/2<1e-4
            y=zeros(0,1);
        else
            y=[];
        end
        yloc=[];
    end
    if size(y,2)==1
        x=K*y+xo;
        val=0.5*x'*Q0*x-b0'*x;
    elseif size(y,2)>1
        %y
        x=K*y+[xo,xo];
        val=0.5*x(:,1)'*Q0*x(:,1)-b0'*x(:,1);
        %v2=0.5*x(:,2)'*Q0*x(:,2)-b0'*x(:,2)
    end
    if size(yloc,2)>0
        ca=K*yloc+xo;
        valc=0.5*ca'*Q0*ca-b0'*ca;
    end    
else
    [x,ca,~,val,valc]=TRS_Loc(Q0,b0,b1,-c1/2);
end
if size(Aeq,1)>0
    if size(x,2)>1
        erg=max(erg,max(abs(Aeq*x(:,1)-beq)));
        if erg>1e-4
            flag_GeSy=1;
            %x(:,1)=x(:,1)-Aeq'*((Aeq*Aeq')\(Aeq*x(:,1)-beq));
        end
        erg=max(erg,max(abs(Aeq*x(:,2)-beq)));
        if erg>1e-4
            flag_GeSy=1;
            %x(:,2)=x(:,2)-Aeq'*((Aeq*Aeq')\(Aeq*x(:,2)-beq));
        end
    elseif size(x,2)==1
        erg=max(erg,max(abs(Aeq*x-beq)));
        if erg>1e-4
            flag_GeSy=1;
            %x(:,1)=x(:,1)-Aeq'*((Aeq*Aeq')\(Aeq*x(:,1)-beq));
        end
        if size(ca,2)>0
            erg=max(erg,max(abs(Aeq*ca-beq)));
            if erg>1e-4
                flag_GeSy=1;
                %ca=ca-Aeq'*((Aeq*Aeq')\(Aeq*ca-beq));
            end
        end
    end
end
    
    