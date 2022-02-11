%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This program performs hydrothermal coordination 
% Hydrothermal TS2
% Test System-2 (two hydro and one thermal unit) in 
% "Thang Trung Nguyen, Dieu Ngoc Vo, Anh Viet Truong, Cuckoo search algorithm
% for short-term hydrothermal scheduling, Applied Energy, Volume 132, 2014, 
% Pages 276-287, ISSN 0306-2619, https://doi.org/10.1016/j.apenergy.2014.07.017"
% $Author: Dr. Rajat Kanti Samal$ $Date: 11-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$

clc;
clear;

NG1=2;% Number of hydro generators
NG2=1;% Number of thermal generators
NG=NG1+NG2;
%% Hydro discharge coefficients
%ah+bhPgi+chPgi2
ah=[0.2 0.4];
bh=[0.03 0.06];
ch=[0.00005 0.0001];
Wj=[25 35];M=24;
totPd=[30	33	35	38	40	45	50	59	61	58	56	57	60	61	65	68	71	62	55	50	43	33	31	30];
%totPd1=totPd/2;
totPd1=totPd;
Qmin=zeros(M,2);%lower bound
Qmax=zeros(M,2);%upper bound
for m=1:24
    for g=1:NG1
        Qmin(m,g)=ah(g);
        Qmax(m,g)=ah(g)+bh(g)*totPd1(m)+ch(g)*(totPd1(m)^2);
    end
end
%% Intial solution
disp('Now creating initial solution')
Qp0=zeros(M,2);
tmp=zeros(2,1);
for m=1:M%(M-1) 
    for g=1:NG1
        Qp0(m,g)=Qmin(m,g)+rand(1,1)*(Qmax(m,g)-Qmin(m,g));
        tmp(g)=tmp(g)+Qp0(m,g);  
    end
end   
% Qp0(M,1)=Wj(1)-tmp(1);
% Qp0(M,2)=Wj(2)-tmp(2);
Qp01=zeros(2*M,1);
Qp01(1:M)=Qp0(:,1);
Qp01((M+1):2*M)=Qp0(:,2);
disp('End of creation of initial solution...')

%% Linear Constraints
A=[];b=[]; Aeq=zeros(1,2*M);
for m=1:M    
    Aeq(1,m)=1;  
    Aeq(2,M+m)=1;
end
beq=Wj';
options = optimset('Algorithm', 'interior-point','MaxFunEvals',3000);


%% Optimal solution
% [QpMC,OBJmc] = fmincon(@htwTS2objLfn,Qp01,A,b,Aeq,beq,Qmin,Qmax,[],options);
startTime=clock;
[QpMC,OBJmc] = fmincon(@htwTS2objLfn,Qp01,A,b,Aeq,beq,Qmin,Qmax,@htwTS2m1nlc,options);
endTime=clock;
disp(endTime-startTime)

%% Hydro and Thermal Generation
B= [0.0 0.000 0.0
    0.0 0.001 0.0
    0.0 0.000 0.0005];
%a+bPgi+cPgi2; Thermal generation coefficients
a=15;
b=3;
c=0.01;
PgT=zeros(M,NG2);PgH=zeros(M,NG1); F=zeros(M,1);
q=QpMC;
for m=1:M % for 24 hours    
    PgH(m,1)=(-bh(1)+sqrt(bh(1)^2-4*ch(1)*(ah(1)-q(m))))/(2*ch(1));
    PgH(m,2)=(-bh(2)+sqrt(bh(2)^2-4*ch(2)*(ah(2)-q(M+m))))/(2*ch(2));
    PgT(m)=totPd(m)-sum(PgH(m,:));  
    % Loss calculation
    Pg=[PgT(m) PgH(m,:)];
    [ Ploss ] = htwBClossfn( NG,Pg,B );
    PgT(m)=PgT(m)+Ploss;
    % Cost 
    F(m)=a+b*PgT(m)+c*(PgT(m)^2); 
end
Tcost=sum(F)
size(QpMC)
res=[QpMC(1:M) QpMC((M+1):2*M) PgH(:,1) PgH(:,2) PgT F totPd'];

% disp('Now writing results...');
% xlswrite('C:\Users\RKSAMAL\Documents\MATLAB\IOFiles\HTW\HTW.xlsx', res, 'TS2', 'D6:J29');

disp('Now writing results...');
xlswrite('C:\Users\RAJAT\Documents\MATLAB\IOFiles\gitIO\PSO-HTS.xlsx', res, 'TS2', 'D6:J29');



