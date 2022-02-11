function [Tcost]=htwTS2objLfn(Qp)
% objective function for hydrothermal TS2
%   Refer the main function for details
% $Author: Dr. Rajat Kanti Samal$ $Date: 11-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$


M=24; NG1=2;NG2=1;NG=NG1+NG2;
totPd=[30	33	35	38	40	45	50	59	61	58	56	57	60	61	65	68	71	62	55	50	43	33	31	30];

B= [0.0 0.000 0.0
    0.0 0.001 0.0
    0.0 0.000 0.0005];

%% Cost Calculation

%ah+bhPgi+chPgi2; Hydro discharge coefficients
ah=[0.2 0.4];
bh=[0.03 0.06];
ch=[0.00005 0.0001];

%a+bPgi+cPgi2; Thermal generation coefficients
a=15;
b=3;
c=0.01;
PgT=zeros(M,1);PgH=zeros(M,NG1); F=zeros(M,1);
q=Qp;
for m=1:M % for 24 hours    
    PgH(m,1)=(-bh(1)+sqrt(bh(1)^2-4*ch(1)*(ah(1)-q(m))))/(2*ch(1));
    PgH(m,2)=(-bh(2)+sqrt(bh(2)^2-4*ch(2)*(ah(2)-q(M+m))))/(2*ch(2));
    PgT(m)=totPd(m)-sum(PgH(m,:));  
    % Loss calculation
    Pg=[PgT(m) PgH(m,:)];
    [ Ploss ] = htwBClossfn( NG,Pg,B );
    PgT(m)=PgT(m)+Ploss;
    % Cost 
    F(m)=a+b*PgT(m)+c*(PgT(m)^2); 
end
Tcost=sum(F);


% B coefficients

