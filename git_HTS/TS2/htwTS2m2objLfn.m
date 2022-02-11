function [Tcost]=htwTS2m2objLfn(PgH)
% objective function hydrothermal TS-2
%   Refer the main function for details
% $Author: Dr. Rajat Kanti Samal$ $Date: 11-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$


NG1=2;% Number of hydro generators
NG2=1;% Number of thermal generators
NG=NG1+NG2;
M=24;
totPdTS2=[30	33	35	38	40	45	50	59	61	58	56	57	60	61	65	68	71	62	55	50	43	33	31	30];

B= [0.0 0.000 0.0
    0.0 0.001 0.0
    0.0 0.000 0.0005];

%% Cost Calculation
%a+bPgi+cPgi2; Thermal generation coefficients
a=15;
b=3;
c=0.001;

PgT=zeros(M,1);F=zeros(M,1);
for m=1:M % for 24 hours    
    %Thermal Generation
    PgT(m)=totPdTS2(m)-sum(PgH(m,:));  
    % Loss
    Pg=[PgT(m) PgH(m,:)];
    [ Ploss ] = htwBClossfn( NG,Pg,B );
    PgT(m)=PgT(m)+Ploss;
    % Cost 
    F(m)=a+b*PgT(m)+c*(PgT(m)^2); 
end
Tcost=sum(F);

