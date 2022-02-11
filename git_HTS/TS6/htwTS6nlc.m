function [cnl,ceqnl]=htwTS6nlc(qT)
% non linear constraint (power balance) for hydrothermal TS-6 
% Refer the main function for details
% $Author: Dr. Rajat Kanti Samal$ $Date: 11-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$


NG1=2;NG2=2;NGht=NG1+NG2;M=3;
totPd=[900 1200 1100];

B= [14 1.0 1.5 1.5
    1.0 6.0 1.0 1.3
    1.5 1.0 6.8 6.5
    1.5 1.3 6.5 7.0]*(10^(-5));

%a+bPgi+cPgi2
ah=[1.980 0.936];
bh=[0.306 0.612];
ch=[0.000216 0.000360];

PgH=zeros(M,NG1);PgT=zeros(M,NG2);ceqnl=zeros(M,1);
for m=1:M    
    PgH(m,1)=(-bh(1)+sqrt(bh(1)^2-4*ch(1)*(ah(1)-qT(m))))/(2*ch(1));
    PgH(m,2)=(-bh(2)+sqrt(bh(2)^2-4*ch(2)*(ah(2)-qT(M+m))))/(2*ch(2));
    % Thermal Generation
    PgT(m,:)=[qT(2*M+m) qT(3*M+m)];
    %Calculate Loss
    PgHT=[PgT(m,:) PgH(m,:)];
    [ Ploss ] = htwBClossfn( NGht,PgHT,B );
    PgT(m,1)=PgT(m,1)+Ploss;% Assuming generator 1 to be slack generator
    % equality constraint
    ceqnl(m)=sum(PgH(m,:))+sum(PgT(m,:))-(totPd(m)+Ploss);
end
cnl=[];