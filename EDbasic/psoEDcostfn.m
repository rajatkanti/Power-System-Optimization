function [Tcost]=psoEDcostfn(ng,a,b,c,Pg)
% This function calculates cost or emissions for a generation schedule
% Cost function is quadratic a+b*Pgi+c*Pgi^2
% $Author: Dr. Rajat Kanti Samal$ $Date: 24-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%
%@Author: Rajat Kanti Samal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


F=zeros(1,ng); Tcost=0;
for i=1:ng
    F(i)=a(i)+b(i)*Pg(i)+c(i)*(Pg(i)^2);
    Tcost=Tcost+F(i);
end