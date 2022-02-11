function [Tcost]=htwTS7objLfn(qT)
% objective function for hydrothermal TS-7 
%   Refer the main function for details
% $Author: Dr. Rajat Kanti Samal$ $Date: 11-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$

M=4;  Msize=12;
NG1=2; NG2=4; NGht=NG1+NG2;

B= [0.000049 0.000014 0.000015 0.000015 0.000020 0.000017
    0.000014 0.000045 0.000016 0.000020 0.000018 0.000015
    0.000015 0.000016 0.000039 0.000010 0.000012 0.000012
    0.000015 0.000020 0.000010 0.000040 0.000014 0.000010
    0.000020 0.000018 0.000012 0.000014 0.000035 0.000011
    0.000017 0.000015 0.000012 0.000010 0.000011 0.000036];

%% Cost Calculation
%ah+bhPgi+chPgi2; Hydro discharge coefficients
ah=[260 250];
bh=[8.5 9.8];
ch=[0.00986 0.01140];
%a+bPgi+cPgi2; Thermal generation coefficients
a=[10 10 20 20];
b=[3.25	2 1.75 1];
c=[0.0083 0.0037 0.0175 0.0625];
d1=[12 18 16 14];
e1=[0.045 0.037	0.038 0.04];
PgTmin1=[20 30 40 50];
PgH=zeros(M,NG1); PgT=zeros(M,NG2); Fgen=zeros(M,NG2); F=zeros(M,1); 
%Ploss=zeros(M,1);
for m=1:M % for M intervals   
    % Hydro Generation
    PgH(m,1)=(-bh(1)+sqrt(bh(1)^2-4*ch(1)*(ah(1)-qT(m))))/(2*ch(1));
    PgH(m,2)=(-bh(2)+sqrt(bh(2)^2-4*ch(2)*(ah(2)-qT(M+m))))/(2*ch(2));
    % Thermal Generation
    PgT(m,:)=[qT(2*M+m) qT(3*M+m) qT(4*M+m) qT(5*M+m)];
    %Calculate Loss
    PgHT=[PgT(m,:) PgH(m,:)];
    [ Ploss ] = htwBClossfn( NGht,PgHT,B );
    PgT(m,1)=PgT(m,1)+Ploss;% Assuming generator 1 to be slack generator
    [F(m),Fgen(m,:)]=htwTGcostNCfn(NG2,PgT(m,:),a,b,c,d1,e1,PgTmin1);        
end
Tcost=sum(F)*Msize;

%% Inline function for minimum cost of thermal units  
    
end% End of main function


% B coefficients

