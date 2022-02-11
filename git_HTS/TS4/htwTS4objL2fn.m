function [Tcost]=htwTS4objL2fn(Qp)

% objective function for hydrothermal TS-4
%   Refer main function for details
% $Author: Dr. Rajat Kanti Samal$ $Date: 11-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$


global totPd;

M=24; NG1=1;NG2=1;NG=NG1+NG2;
% totPd=[1200 1500];

B= [0.00 0.00000
    0.00 0.00008];

%% Cost Calculation
%a+bPgi+cPgi2; Thermal generation coefficients
a=575;
b=9.2;
c=0.00184;
%ah+bhPgi+chPgi2; Hydro discharge coefficients
ah=330;
bh=4.97;
ch=0;

PgT=zeros(M,1);PgH=zeros(M,1); F=zeros(M,1);
for m=1:M % for M intervals    
    PgH(m)=(Qp(m)-ah)/bh;  % since ch=0
    PgT(m)=totPd(m)-PgH(m);  
    % Loss calculation
    Pg=[PgT(m) PgH(m)];
    [ Ploss ] = htwBClossfn( NG,Pg,B );
    PgT(m)=PgT(m)+Ploss;
    % Cost 
    F(m)=a+b*PgT(m)+c*(PgT(m)^2); 
end
Tcost=sum(F);


% B coefficients

