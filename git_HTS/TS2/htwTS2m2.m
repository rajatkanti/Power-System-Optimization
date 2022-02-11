%% This program performs hydrothermal coordination (Method-2)
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
M=24;
%% Hydro discharge coefficients
%ah+bhPgi+chPgi2
ah=[0.2 0.4];
bh=[0.03 0.06];
ch=[0.00005 0.0001];
Wj=[25 35];
totPdTS2=[30 33	35	38	40	45	50	59	61	58	56	57	60	61	65	68	71	62	55	50	43	33	31	30];
PgHmin=zeros(M,NG1);%lower bound
PgHmax=zeros(M,NG1);%upper bound
for m=1:M
    for g=1:NG1
        PgHmin(m,g)=0;
        PgHmax(m,g)=totPdTS2(m);
    end
end
%% Intial solution
disp('Now creating initial solution')
PgH0=zeros(M,1);
tmp=zeros(2,1);
for m=1:M  
    for g=1:NG1
        PgH0(m,g)=PgHmin(m,g)+rand(1,1)*(PgHmax(m,g)-PgHmin(m,g));
        tmp(g)=tmp(g)+PgH0(m,g); 
    end
end   
disp('End of creation of initial solution...')

%% Linear Constraints
A=[];b=[]; Aeq=[];beq=[];
options = optimset('Algorithm', 'interior-point');


%% Optimal solution
[PgHMC,OBJmc] = fmincon(@htwTS2m2objLfn,PgH0,A,b,Aeq,beq,PgHmin,PgHmax,@htwTS2nlc,options);

%% Hydro and Thermal Generation
B= [0.0 0.000 0.0
    0.0 0.001 0.0
    0.0 0.000 0.0005];
%a+bPgi+cPgi2; Thermal generation coefficients
a=15;
b=3;
c=0.001;
PgH=PgHMC;PgT=zeros(M,1);F=zeros(M,1);q=zeros(M,NG2);
for m=1:M
    % Discharge
    for g=1:NG1
        q(m,g)=ah(g)+bh(g)*PgH(m,g)+ch(g)*(PgH(m,g)^2);
    end
    % Thermal Generation
    PgT(m)=totPdTS2(m)-sum(PgH(m,:));        
    % Loss
    Pg=[PgT(m) PgH(m,:)];
    [ Ploss ] = htwBClossfn( NG,Pg,B );
    PgT(m)=PgT(m)+Ploss;
    % Cost 
    F(m)=a+b*PgT(m)+c*(PgT(m)^2);         
end
Tcost=sum(F)

res=[q(:,1) q(:,2) PgH(:,1) PgH(:,2) PgT F totPdTS2'];

disp('Now writing results...');
xlswrite('C:\Users\RAJAT\Documents\MATLAB\IOFiles\gitIO\PSO-HTS.xlsx', res, 'TS2', 'M6:S29');



