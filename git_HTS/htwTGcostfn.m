function [Tcost,F]=htwTGcostfn(NG2,Pg,a,b,c)   
%Thermal generation cost calculation
%   Fuel cost
    Tcost=0;
    F=zeros(1,NG2); 
    for g=1:NG2
        F(g)=a(g)+b(g)*Pg(g)+c(g)*(Pg(g)^2);
        Tcost=Tcost+F(g);
    end

end


