function [Tcost]=htwTS1objLfn(Qp)
% hydrothermal TS-1 objective function
% Refer main function for test system details
% $Author: Dr. Rajat Kanti Samal$ $Date: 11-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$

global totPd


M=24; NG1=1;NG2=1;NG=NG1+NG2;
% totPd=[455 425	415	407	400	420	487	604	665	675	695	705	580	605	616	653	721	740	700	678	630	585	540	503];

B= [0.00005 0.00001
    0.00001 0.00015];

%% Cost Calculation
%a+bPgi+cPgi2; Thermal generation coefficients
a=373.7;
b=9.606;
c=0.001991;

%ah+bhPgi+chPgi2; Hydro discharge coefficients
ah=61.53;
bh=-0.009079;
ch=0.0007749;

PgT=zeros(M,1);PgH=zeros(M,1); F=zeros(M,1);
q=Qp;
for m=1:M % for 24 hours    
    PgH(m)=(-bh+sqrt(bh^2-4*ch*(ah-q(m))))/(2*ch);    
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

