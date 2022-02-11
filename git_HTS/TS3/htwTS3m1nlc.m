function [cnl,ceqnl]=htwTS3m1nlc(qT)
% non linear constraint function for hydrotherma TS-3 
%   Refer the reference provided in main function for details
% $Author: Dr. Rajat Kanti Samal$ $Date: 11-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$


NG1=2;NG2=2;NGht=NG1+NG2;M=24;
totPd=[400	300	250	250	250	300	450	900	1230 1250 1350 1400 1200 1250 1250 1270 1350 1470 1330 1250 1170 1050 900 600];

B= [0.00014 0.000010 0.000015 0.000015
    0.000010 0.00006 0.000010 0.000013
    0.000015 0.000010 0.000068 0.000065
    0.000015 0.000013 0.000065 0.00007];

%a+bPgi+cPgi2
ah=[1.98 0.936];
bh=[0.306	0.612];
ch=[0.000216	0.00036];

PgH=zeros(M,NG1);PgT=zeros(M,NG2);ceqnl=zeros(M,1);
for m=1:M    
    PgH(m,1)=(-bh(1)+sqrt(bh(1)^2-4*ch(1)*(ah(1)-qT(m))))/(2*ch(1));
    PgH(m,2)=(-bh(2)+sqrt(bh(2)^2-4*ch(2)*(ah(2)-qT(M+m))))/(2*ch(2));
    PgT(m,:)=[qT(2*M+m) qT(3*M+m)];
    %Calculate Loss
    PgHT=[PgT(m,:) PgH(m,:)];
    [ Ploss ] = htwBClossfn( NGht,PgHT,B );
    PgT(m,1)=PgT(m,1)+Ploss;% Assuming generator 1 to be slack generator
    % equality constraint
    ceqnl(m)=sum(PgH(m,:))+sum(PgT(m,:))-(totPd(m)+Ploss);
end
cnl=[];