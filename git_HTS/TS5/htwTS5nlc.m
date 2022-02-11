function [cnl,ceqnl]=htwTS5nlc(qT)
% non linear constraint for hydrothermal test system-5
%   Refer main function for details
% $Author: Dr. Rajat Kanti Samal$ $Date: 11-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$


NG1=2;NG2=2;NGht=NG1+NG2;M=4;
totPd=[1200 1500 1400 1700];

B= [4.0 1.0 1.5 1.5
    1.0 3.5 1.0 1.2
    1.5 1.0 3.9 2.0
    1.5 1.2 2.0 4.9]*(10^(-5));

%a+bPgi+cPgi2
ah=[260 250];
bh=[8.5 9.8];
ch=[0.00986 0.0114];

PgH=zeros(M,NG1);PgT=zeros(M,NG2);ceqnl=zeros(M,1);
for m=1:M    
    PgH(m,1)=(-bh(1)+sqrt(bh(1)^2-4*ch(1)*(ah(1)-qT(m))))/(2*ch(1));
    PgH(m,2)=(-bh(2)+sqrt(bh(2)^2-4*ch(2)*(ah(2)-qT(M+m))))/(2*ch(2));
    PgT(m,:)=[qT(2*M+m) qT(3*M+m)];
    % Calculate Loss
    PgHT=[PgT(m,:) PgH(m,:)];
    [ Ploss ] = htwBClossfn( NGht,PgHT,B );
    PgT(m,1)=PgT(m,1)+Ploss;% Assuming generator 1 to be slack generator
    % equality constraint
    ceqnl(m)=sum(PgH(m,:))+sum(PgT(m,:))-(totPd(m)+Ploss);
end
cnl=[];