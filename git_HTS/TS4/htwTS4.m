%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This program performs hydrothermal coordination 
% Hydrothermal TS4
% Test System-4 (one hydro and one thermal unit) in 
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
ah=330;
bh=4.97;
ch=0;
Wj=100000;
choice=2;
switch(choice)
    case 1
        M=2; Msize=12;
        totPd=[1200 1500];
        outCol='D6:H7';
    case 2
        M=24; Msize=1;
        totPd=[1200 1200 1200 1200 1200 1200 1200 1200 1200 1200 1200 1200 1500 1500 1500 1500 1500 1500 1500 1500 1500 1500 1500 1500];
        outCol='L6:P29';
end
Qmin=zeros(M,1);%lower bound
Qmax=zeros(M,1);%upper bound
for m=1:M
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
Qp0(M)=(Wj/Msize)-tmp;
% To ensure feasibility of initial solution
while Qp0(M)<Qmin(M) || Qp0(M)>Qmax(M)
    tmp=0;
    for h=1:(M-1)      
        Qp0(h)=Qmin(h)+rand(1,1)*(Qmax(h)-Qmin(h));
        tmp=tmp+Qp0(h);      
    end   
    Qp0(M)=(Wj/Msize)-tmp;
end
disp('End of creation of initial solution...')

%% Linear Constraints
A=[];b1=[];Aeq=zeros(1,M);
for h=1:M
    Aeq(h)=Msize;
end
beq=Wj;
options = optimset('Algorithm', 'interior-point');


%% Optimal solution
startTime=clock;
switch(choice)
    case 1
        [QpMC,OBJmc] = fmincon(@htwTS4objLfn,Qp0,A,b1,Aeq,beq,Qmin,Qmax,[],options);
    case 2
        [QpMC,OBJmc] = fmincon(@htwTS4objL2fn,Qp0,A,b1,Aeq,beq,Qmin,Qmax,[],options);
end

endTime=clock;
disp(endTime-startTime)

%% Hydro and Thermal Generation
NG1=1;% Number of hydro generators
NG2=1;% Number of thermal generators
NG=NG1+NG2;
B= [0.00 0.00000
    0.00 0.00008];
% Thermal generator coefficients
a=575;
b=9.2;
c=0.00184;
PgH=zeros(M,1);PgT=zeros(M,1);F=zeros(M,1); 
Qp=QpMC;
for m=1:M
    PgH(m)=(Qp(m)-ah)/bh;  % since ch=0
    PgT(m)=totPd(m)-PgH(m);  
    % Loss calculation
    Pg=[PgT(m) PgH(m)];
    [ Ploss ] = htwBClossfn( NG,Pg,B );
    PgT(m)=PgT(m)+Ploss;
    % Cost 
    F(m)=a+b*PgT(m)+c*(PgT(m)^2);   
end
Tcost=sum(F)*Msize;

res=[QpMC*Msize PgH PgT F*Msize totPd'];

% disp('Now writing results...');
% xlswrite('C:\Users\RAJAT\Documents\MATLAB\IOFiles\gitIO\HTW.xlsx', res, 'TS4', outCol);

disp('Now writing results...');
xlswrite('C:\Users\RAJAT\Documents\MATLAB\IOFiles\gitIO\PSO-HTS.xlsx', res, 'TS4', outCol);



