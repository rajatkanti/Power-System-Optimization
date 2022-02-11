%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This program performs hydrothermal coordination 
% Hydrothermal TS1 
% Test System-1 in 
% "Thang Trung Nguyen, Dieu Ngoc Vo, Anh Viet Truong, Cuckoo search algorithm
% for short-term hydrothermal scheduling, Applied Energy, Volume 132, 2014, 
% Pages 276-287, ISSN 0306-2619, https://doi.org/10.1016/j.apenergy.2014.07.017"
% $Author: Dr. Rajat Kanti Samal$ $Date: 11-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$

clc;
clear;

global totPd;

%% Hydro discharge coefficients
%ah+bhPgi+chPgi2
ah=61.53;
bh=-0.009079;
ch=0.0007749;
Wj=2559.6;M=24;
totPd=[455 425	415	407	400	420	487	604	665	675	695	705	580	605	616	653	721	740	700	678	630	585	540	503];
Qmin=zeros(M,1);%lower bound
Qmax=zeros(M,1);%upper bound
for m=1:24
    Qmin(m)=ah;
    Qmax(m)=ah+bh*totPd(m)+ch*(totPd(m)^2);
end
%% Intial solution
disp('Now creating initial solution')
Qp0=zeros(M,1);
tmp=0;
for h=1:(M-1)      
    Qp0(h)=Qmin(h)+rand(1,1)*(Qmax(h)-Qmin(h));
    tmp=tmp+Qp0(h);      
end   
Qp0(M)=Wj-tmp;
disp('End of creation of initial solution...')

%% Linear Constraints
A=[];b=[]; Aeq=zeros(1,M);
for h=1:M
    Aeq(h)=1;
end
beq=Wj;
options = optimset('Algorithm', 'interior-point');


%% Optimal solution
startTime=clock;
[QpMC,OBJmc] = fmincon(@htwTS1objLfn,Qp0,A,b,Aeq,beq,Qmin,Qmax,[],options);
endTime=clock;
disp(endTime-startTime)

%% Hydro and Thermal Generation
NG1=1;% Number of hydro generators
NG2=1;% Number of thermal generators
NG=NG1+NG2;
B= [0.00005 0.00001
    0.00001 0.00015];
a=373.7;
b=9.606;
c=0.001991;
PgH=zeros(M,1);PgT=zeros(M,1);F=zeros(M,1); 
q=QpMC;
for m=1:M
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

res=[QpMC PgH PgT F totPd'];

% disp('Now writing results...');
% xlswrite('C:\Users\RKSAMAL\Documents\MATLAB\IOFiles\HTW\HTW.xlsx', res, 'TS1', 'D6:H29');

disp('Now writing results...');
xlswrite('C:\Users\RKSAMAL\Documents\MATLAB\IOFiles\gitIO\PSO-HTS.xlsx', res, 'TS1', 'D6:H29');



