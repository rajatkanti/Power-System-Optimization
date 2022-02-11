%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This program performs hydrothermal coordination 
% Hydrothermal TS7 (non-convex)
% Test System-7 (two hydro and four thermal units) in 
% "Thang Trung Nguyen, Dieu Ngoc Vo, Anh Viet Truong, Cuckoo search algorithm
% for short-term hydrothermal scheduling, Applied Energy, Volume 132, 2014, 
% Pages 276-287, ISSN 0306-2619, https://doi.org/10.1016/j.apenergy.2014.07.017"
% $Author: Dr. Rajat Kanti Samal$ $Date: 11-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$


clc;
clear;

NG1=2;% Number of hydro generators
NG2=4;% Number of thermal generators
NGht=NG1+NG2;
M=4;Msize=12;
%% Hydro and Thermal Generation Limits
% Hydro Generation Coefficients
ah=[260 250];
bh=[8.5 9.8];
ch=[0.00986 0.01140];
Wj=[125000 286000];
% Thermal generation coefficients (a+bPgi+cPgi2)
a=[10 10 20 20];
b=[3.25	2 1.75 1];
c=[0.0083 0.0037 0.0175 0.0625];
d1=[12 18 16 14];
e1=[0.045 0.037	0.038 0.04];
% Load demand
totPd=[900 1100 1000 1300];
% Limits
Qmin=zeros(M,NG1);Qmax=zeros(M,NG1);
Qmax1=zeros(M,NG1);Qmax2=zeros(M,NG1);PgHmax=[250 500];
PgTmin=zeros(M,NG2);PgTmax=zeros(M,NG2);
PgTmin1=[20 30 40 50];PgTmax1=[125 175 250 300];
for m=1:M
    for g=1:NG1
        Qmin(m,g)=ah(g);
        Qmax1(m,g)=ah(g)+bh(g)*totPd(m)+ch(g)*(totPd(m)^2);
        Qmax2(m,g)=ah(g)+bh(g)*PgHmax(g)+ch(g)*(PgHmax(g)^2);
        Qmax(m,g)=min(Qmax1(m,g),Qmax2(m,g));
        
    end
    for g=1:NG2
        PgTmin(m,g)=PgTmin1(g);
        PgTmax(m,g)=PgTmax1(g);
    end
end
xMin=[Qmin(:,1); Qmin(:,2); PgTmin(:,1); PgTmin(:,2); PgTmin(:,3); PgTmin(:,4)];
xMax=[Qmax(:,1); Qmax(:,2); PgTmax(:,1); PgTmax(:,2); PgTmax(:,3); PgTmax(:,4)];
%% Intial solution
disp('Now creating initial solution')
Qp0=zeros(M,NG1);PgT0=zeros(M,NG2);
tmp=zeros(NG1,1);tmpT=zeros(NG2,1);
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
decVar0((4*M+1):5*M)=PgT0(:,3);
decVar0((5*M+1):6*M)=PgT0(:,4);
disp('End of creation of initial solution...')
% decVar0
%% Linear Constraints
A=[];b1=[]; Aeq=zeros(NG1,NGht*M);
for m=1:M    
    Aeq(1,m)=Msize;  
    Aeq(2,M+m)=Msize;
end
beq=Wj';
options = optimset('Algorithm', 'interior-point','MaxFunEvals',6000);


%% Optimal solution
startTime=clock;
[qT,OBJmc] = fmincon(@htwTS7objLfn,decVar0,A,b1,Aeq,beq,xMin,xMax,@htwTS7nlc,options);
endTime=clock;
disp(endTime-startTime)

%% Hydro and Thermal Generation
B= [0.000049 0.000014 0.000015 0.000015 0.000020 0.000017
    0.000014 0.000045 0.000016 0.000020 0.000018 0.000015
    0.000015 0.000016 0.000039 0.000010 0.000012 0.000012
    0.000015 0.000020 0.000010 0.000040 0.000014 0.000010
    0.000020 0.000018 0.000012 0.000014 0.000035 0.000011
    0.000017 0.000015 0.000012 0.000010 0.000011 0.000036];
PgH=zeros(M,NG1); PgT=zeros(M,NG2); Fgen=zeros(M,NG2); F=zeros(M,1); 
for m=1:M % for M intervals  
    PgH(m,1)=(-bh(1)+sqrt(bh(1)^2-4*ch(1)*(ah(1)-qT(m))))/(2*ch(1));
    PgH(m,2)=(-bh(2)+sqrt(bh(2)^2-4*ch(2)*(ah(2)-qT(M+m))))/(2*ch(2));
    % Thermal Generation
    PgT(m,:)=[qT(2*M+m) qT(3*M+m) qT(4*M+m) qT(5*M+m)];
    %Calculate Loss
    PgHT=[PgT(m,:) PgH(m,:)];
    [ Ploss ] = htwBClossfn( NGht,PgHT,B );
    PgT(m,1)=PgT(m,1)+Ploss;% Assuming generator 1 to be slack generator
    [F(m),Fgen(m,:)]=htwTGcostNCfn(NG2,PgT(m,:),a,b,c,d1,e1,PgTmin(m,:)); 
end
Tcost=sum(F)*Msize
size(qT)
res=[qT(1:M)*Msize qT((M+1):2*M)*Msize PgH(:,1) PgH(:,2) PgT(:,1) PgT(:,2) ...
    PgT(:,3) PgT(:,4) Fgen(:,1)*Msize Fgen(:,2)*Msize Fgen(:,3)*Msize Fgen(:,4)*Msize totPd'];

% disp('Now writing results...');
% xlswrite('C:\Users\RKSAMAL\Documents\MATLAB\IOFiles\HTW\HTW.xlsx', res, 'TS7', 'D6:P9');

disp('Now writing results...');
xlswrite('C:\Users\RAJAT\Documents\MATLAB\IOFiles\gitIO\PSO-HTS.xlsx', res, 'TS7', 'D6:P9');




