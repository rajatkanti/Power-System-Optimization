%% This program performs hydrothermal coordination (Method 2)
% Hydrothermal TS1 
% Test System-1 in 
% "Thang Trung Nguyen, Dieu Ngoc Vo, Anh Viet Truong, Cuckoo search algorithm
% for short-term hydrothermal scheduling, Applied Energy, Volume 132, 2014, 
% Pages 276-287, ISSN 0306-2619, https://doi.org/10.1016/j.apenergy.2014.07.017"
% Method-2: In this method, hydro generations are variables
% $Author: Dr. Rajat Kanti Samal$ $Date: 11-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$


clc;
clear;

%% Hydro discharge coefficients
%ah+bhPgi+chPgi2
ah=61.53;
bh=-0.009079;
ch=0.0007749;
Wj=2559.6;M=24;
totPd=[455 425	415	407	400	420	487	604	665	675	695	705	580	605	616	653	721	740	700	678	630	585	540	503];
PgHmin=zeros(M,1);%lower bound
PgHmax=zeros(M,1);%upper bound
for m=1:24
    PgHmin(m)=0;
    PgHmax(m)=totPd(m);
end
%% Intial solution
disp('Now creating initial solution')
PgH0=zeros(M,1);
tmp=0;
for h=1:M      
    PgH0(h)=PgHmin(h)+rand(1,1)*(PgHmax(h)-PgHmin(h));
    tmp=tmp+PgH0(h);      
end   
disp('End of creation of initial solution...')

%% Linear Constraints
A=[];b=[]; Aeq=[];beq=[];
options = optimset('Algorithm', 'interior-point');


%% Optimal solution
[PgHMC,OBJmc] = fmincon(@htwTS1m2objLfn,PgH0,A,b,Aeq,beq,PgHmin,PgHmax,@htwTS1nlc,options);

%% Hydro and Thermal Generation
NG1=1;% Number of hydro generators
NG2=1;% Number of thermal generators
NG=NG1+NG2;
B= [0.00005 0.00001
    0.00001 0.00015];
a=373.7;
b=9.606;
c=0.001991;
PgH=PgHMC;PgT=zeros(M,1);F=zeros(M,1);q=zeros(m,1);
for m=1:M
    % Discharge
    q(m)=ah+bh*PgH(m)+ch*(PgH(m)^2);
    % Thermal Generation
    PgT(m)=totPd(m)-PgH(m);        
    % Loss
    Pg=[PgT PgH];
    [ Ploss ] = htwBClossfn( NG,Pg,B );
    PgT(m)=PgT(m)+Ploss;
    % Cost 
    F(m)=a+b*PgT(m)+c*(PgT(m)^2);         
end
Tcost=sum(F)

res=[q PgH PgT F totPd'];

% disp('Now writing results...');
% xlswrite('C:\Users\Rajat\Documents\MATLAB\IOFiles\HTW\HTW.xlsx', res, 'TS1', 'L6:P29');



