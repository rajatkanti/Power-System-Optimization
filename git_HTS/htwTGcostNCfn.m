function [Tcost,F]=htwTGcostNCfn(NG2,Pg,a,b,c,d1,e1,Pmin)   
% Thermal generation cost calculation
%   Fuel cost
    F=zeros(1,NG2); Tcost=0;
    for i=1:NG2
        F(i)=a(i)+b(i)*Pg(i)+c(i)*(Pg(i)^2)+abs(d1(i)*sin(e1(i)*(Pmin(i)-Pg(i))));
        Tcost=Tcost+F(i);
    end

end


