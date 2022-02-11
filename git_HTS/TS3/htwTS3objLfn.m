function [Tcost]=htwTS3objLfn(qT)
% objective function for hydrothermal TS-3 
%   Refer the reference provided in main function for details
% $Author: Dr. Rajat Kanti Samal$ $Date: 11-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$


M=24; NG1=2; NG2=2; NGht=NG1+NG2;

B= [0.00014 0.000010 0.000015 0.000015
    0.000010 0.00006 0.000010 0.000013
    0.000015 0.000010 0.000068 0.000065
    0.000015 0.000013 0.000065 0.00007];

%% Cost Calculation
%ah+bhPgi+chPgi2; Hydro discharge coefficients
ah=[1.98 0.936];
bh=[0.306	0.612];
ch=[0.000216	0.00036];
%a+bPgi+cPgi2; Thermal generation coefficients
a=[25	30];
b=[3.2	3.4];
c=[0.0025	0.0008];
PgH=zeros(M,NG1); PgT=zeros(M,NG2); F=zeros(M,1); Fgen=zeros(M,NG2);
%Ploss=zeros(M,1);
for m=1:M % for 24 hours    
    % Hydro Generation
    PgH(m,1)=(-bh(1)+sqrt(bh(1)^2-4*ch(1)*(ah(1)-qT(m))))/(2*ch(1));
    PgH(m,2)=(-bh(2)+sqrt(bh(2)^2-4*ch(2)*(ah(2)-qT(M+m))))/(2*ch(2));
    % Thermal Generation
    PgT(m,:)=[qT(2*M+m) qT(3*M+m)];
    %Calculate Loss
    PgHT=[PgT(m,:) PgH(m,:)];
    [ Ploss ] = htwBClossfn( NGht,PgHT,B );
    PgT(m,1)=PgT(m,1)+Ploss;% Assuming generator 1 to be slack generator
    [F(m),Fgen(m,:)]=htwTGcostfn(NG2,PgT(m,:),a,b,c);        
end
Tcost=sum(F);

%% Inline function for minimum cost of thermal units  
    
end% End of main function


% B coefficients

