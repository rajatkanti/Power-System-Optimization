%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This program performs hydrothermal coordination 
% Hydrothermal TS3
% Test System-3 (two hydro and two thermal units) in 
% "Thang Trung Nguyen, Dieu Ngoc Vo, Anh Viet Truong, Cuckoo search algorithm
% for short-term hydrothermal scheduling, Applied Energy, Volume 132, 2014, 
% Pages 276-287, ISSN 0306-2619, https://doi.org/10.1016/j.apenergy.2014.07.017"
% $Author: Dr. Rajat Kanti Samal$ $Date: 11-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$

clc;
clear;

NG1=2;% Number of hydro generators
NG2=2;% Number of thermal generators
NGht=NG1+NG2;
M=24;
%% Hydro and Thermal Generation Limits
% Hydro Generation Coefficients
ah=[1.98 0.936];
bh=[0.306	0.612];
ch=[0.000216	0.00036];
Wj=[2500	2100];
totPd=[400	300	250	250	250	300	450	900	1230 1250 1350 1400 1200 1250 1250 1270 1350 1470 1330 1250 1170 1050 900 600];
Qmin=zeros(M,NG1);Qmax=zeros(M,NG1);
PgTmin=zeros(M,NG2);PgTmax=zeros(M,NG2);
for m=1:24
    for g=1:NG1
        Qmin(m,g)=ah(g);
        Qmax(m,g)=ah(g)+bh(g)*totPd(m)+ch(g)*(totPd(m)^2);
    end
    for g=1:NG2
        PgTmin(m,g)=0;
        PgTmax(m,g)=totPd(m);
    end
end
xMin=[Qmin(:,1); Qmin(:,2); PgTmin(:,1); PgTmin(:,2)];
xMax=[Qmax(:,1); Qmax(:,2); PgTmax(:,1); PgTmax(:,2)];
%% Intial solution
disp('Now creating initial solution')
Qp0=zeros(M,NG1);PgT0=zeros(M,NG2);
tmp=zeros(2,1);tmpT=zeros(2,1);
for m=1:M%(M-1) 
    for g=1:NG1
        Qp0(m,g)=Qmin(m,g)+rand(1,1)*(Qmax(m,g)-Qmin(m,g));
        tmp(g)=tmp(g)+Qp0(m,g);  
    end
    for g=1:NG2
        PgT0(m,g)=PgTmin(m,g)+rand(1,1)*(PgTmax(m,g)-PgTmin(m,g));
        tmpT(g)=tmpT(g)+PgT0(m,g);  
    end
end   
decVar0=zeros(NGht*M,1);
decVar0(1:M)=Qp0(:,1);
decVar0((M+1):2*M)=Qp0(:,2);
decVar0((2*M+1):3*M)=PgT0(:,1);
decVar0((3*M+1):4*M)=PgT0(:,2);
disp('End of creation of initial solution...')
% decVar0
%% Linear Constraints
A=[];b1=[]; Aeq=zeros(1,NGht*M);
for m=1:M    
    Aeq(1,m)=1;  
    Aeq(2,M+m)=1;
end
beq=Wj';
options = optimset('Algorithm', 'interior-point','MaxFunEvals',50000);


%% Optimal solution
startTime=clock;
[qT,OBJmc] = fmincon(@htwTS3objLfn,decVar0,A,b1,Aeq,beq,xMin,xMax,@htwTS3m1nlc,options);
endTime=clock;
disp(endTime-startTime)

%% Hydro and Thermal Generation
B= [0.00014 0.000010 0.000015 0.000015
    0.000010 0.00006 0.000010 0.000013
    0.000015 0.000010 0.000068 0.000065
    0.000015 0.000013 0.000065 0.00007];
%a+bPgi+cPgi2; Thermal generation coefficients
a=[25	30];
b=[3.2	3.4];
c=[0.0025	0.0008];

PgH=zeros(M,NG1); PgT=zeros(M,NG2); F=zeros(M,1); Fgen=zeros(M,NG2);
for m=1:M % for 24 hours    
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
% Tcost=sum(F)
% size(qT)
res=[qT(1:M) qT((M+1):2*M) PgH(:,1) PgH(:,2) PgT(:,1) PgT(:,2) Fgen(:,1) Fgen(:,2) totPd'];

% disp('Now writing results...');
% xlswrite('C:\Users\RKSAMAL\Documents\MATLAB\IOFiles\HTW\HTW.xlsx', res, 'TS3', 'D6:L29');

disp('Now writing results...');
xlswrite('C:\Users\RAJAT\Documents\MATLAB\IOFiles\gitIO\PSO-HTS.xlsx', res, 'TS3', 'D6:L29');



