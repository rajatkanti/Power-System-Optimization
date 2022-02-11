function [Tcost]=htwTS6objLfn(qT)
% objective function for hydrothermal TS-6 
% Refer the main function for details
% $Author: Dr. Rajat Kanti Samal$ $Date: 11-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$

M=3;  Msize=8;
NG1=2; NG2=2; NGht=NG1+NG2;

B= [14 1.0 1.5 1.5
    1.0 6.0 1.0 1.3
    1.5 1.0 6.8 6.5
    1.5 1.3 6.5 7.0]*(10^(-5));

%% Cost Calculation
%ah+bhPgi+chPgi2; Hydro discharge coefficients
ah=[1.980 0.936];
bh=[0.306 0.612];
ch=[0.000216 0.000360];
%a+bPgi+cPgi2; Thermal generation coefficients
a=[25 30];
b=[3.2 3.4];
c=[0.0025	0.0008];
d1=[12 14];
e1=[0.0550 0.0450];
PgTmin1=[50 50];
PgH=zeros(M,NG1); PgT=zeros(M,NG2); Fgen=zeros(M,NG2); F=zeros(M,1); 
%Ploss=zeros(M,1);
for m=1:M % for M intervals   
    % Hydro Generation
    PgH(m,1)=(-bh(1)+sqrt(bh(1)^2-4*ch(1)*(ah(1)-qT(m))))/(2*ch(1));
    PgH(m,2)=(-bh(2)+sqrt(bh(2)^2-4*ch(2)*(ah(2)-qT(M+m))))/(2*ch(2));
    % Thermal Generation
    PgT(m,:)=[qT(2*M+m) qT(3*M+m)];
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

