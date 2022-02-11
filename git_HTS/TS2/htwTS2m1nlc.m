function [cnl,ceqnl]=htwTS2m1nlc(Qp)
% non-linear constraint function for hydrothermal TS-2 
%   Refer the main function for details
% $Author: Dr. Rajat Kanti Samal$ $Date: 11-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$

NG1=2;M=24;
totPd=[30	33	35	38	40	45	50	59	61	58	56	57	60	61	65	68	71	62	55	50	43	33	31	30];

%a+bPgi+cPgi2
ah=[0.2 0.4];
bh=[0.03 0.06];
ch=[0.00005 0.0001];

PgH=zeros(M,NG1);cnl=zeros(M,1);
q=Qp; 
for m=1:M    
    PgH(m,1)=(-bh(1)+sqrt(bh(1)^2-4*ch(1)*(ah(1)-q(m))))/(2*ch(1));
    PgH(m,2)=(-bh(2)+sqrt(bh(2)^2-4*ch(2)*(ah(2)-q(M+m))))/(2*ch(2));
    cnl(m)=sum(PgH(m,:))-totPd(m);
end

ceqnl=[];