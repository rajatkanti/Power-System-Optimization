%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program solves minimum cost and emission dispatch for the problem in
% Power System Optimization by Kothari,Dhillon, P-346 
% using Weighting method (P-354)
% $Author: Dr. Rajat Kanti Samal$ $Date: 24-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$
%@Author: Rajat Kanti Samal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;

totPd=1800;% load demand in MW
delP=totPd; % intial error 
ng=6; % number of generators
nb=6; % number of buses
% Fuel cost curves of the generators a+bPgi+cPgi2
maxOption=20; 
alpha=0.001; % Convergece speed of iterative method

%% Cost and emission coefficients
disp('Now reading cost and emission coefficients...');
coeffAll=xlsread('C:\Users\RAJAT\Documents\MATLAB\IOFiles\OPTio\PSO-MO.xlsx', 'System', 'D6:Q11');
% Cost coefficients
a1=coeffAll(:,1); b1=coeffAll(:,2); c1=coeffAll(:,3);
% NOx emission coefficients
a2=coeffAll(:,4); b2=coeffAll(:,5); c2=coeffAll(:,6);
% SO2 emission coefficients
a3=coeffAll(:,7); b3=coeffAll(:,8); c3=coeffAll(:,9);
% CO2 emission coefficients
a4=coeffAll(:,10); b4=coeffAll(:,11); c4=coeffAll(:,12);
% Generator limits
genLimMin=coeffAll(:,13);
genLimMax=coeffAll(:,14);

%% Weights
disp('Now reading weights...');
allWeights=xlsread('C:\Users\RAJAT\Documents\MATLAB\IOFiles\OPTio\PSO-MO.xlsx', 'Hourly', 'B4:E23');
w1=allWeights(:,1); w2=allWeights(:,2); w3=allWeights(:,3); w4=allWeights(:,4);

%% B-coefficients
B=[0.000200 0.000010 0.000015 0.000005 0.000000 -0.00003
   0.000010 0.000300 -0.00002 0.000001 0.000012 0.000010
   0.000015 -0.00002 0.000100 -0.00001 0.000010 0.000008
   0.000005 0.000001 -0.00001 0.000150 0.000006 0.000050
   0.000000 0.000012 0.000010 0.000006 0.000250 0.000020
   -0.00003 0.000010 0.000008 0.000050 0.000020 0.000210];

nfSolution=zeros(maxOption,12);
Fmin=zeros(4,1);Fmax=zeros(4,1);

for dmOption=1:maxOption
    
    %% Obtain weighted cost/emission coefficients
    a=w1(dmOption)*a1+w2(dmOption)*a2+w3(dmOption)*a3+w4(dmOption)*a4;
    b=w1(dmOption)*b1+w2(dmOption)*b2+w3(dmOption)*b3+w4(dmOption)*b4;
    c=w1(dmOption)*c1+w2(dmOption)*c2+w3(dmOption)*c3+w4(dmOption)*c4;
    
    
    %% Economic schedule without losses
    [Pg,lambda]=psoEDnoLossfn(ng,b,c,totPd); % to be used as input for with loss    
    [Pg,lambda,Pl]=psoEDwithLossfn(totPd,ng,nb,b,c,B,Pg,lambda,alpha);
    
    %% Compute objective function values
    F=zeros(4,1); % Values of four objective functions
    a=a1;b=b1;c=c1;
    F(1)=psoEDcostfn(ng,a,b,c,Pg);% Total cost in Rs/h 
    a=a2;b=b2;c=c2;% Total NOx emissions in kg/h
    F(2)=psoEDcostfn(ng,a,b,c,Pg);
    a=a3;b=b3;c=c3;% Total SO2 emissions in kg/h
    F(3)=psoEDcostfn(ng,a,b,c,Pg);
    a=a4;b=b4;c=c4;% Total CO2 emissions in kg/h
    F(4)=psoEDcostfn(ng,a,b,c,Pg);
    
    %% Solution vector
    nfSolution(dmOption,:)=[Pg(1) Pg(2) Pg(3) Pg(4) Pg(5) Pg(6) F(1) F(2) F(3) F(4)...
        lambda Pl];
    
    %% Compute membership function of solution
    %% the following code segment calculates min and max value of objectives    
    if dmOption==1
        %initialize min and max values to that of first option
        for objNo=1:4
            Fmin(objNo)=F(objNo);
            Fmax(objNo)=F(objNo);
        end
    else
        %if objective function is less than min or more than max change...
        for objNo=1:4
            if F(objNo)<Fmin(objNo)
                Fmin(objNo)=F(objNo);
            end
            if F(objNo)>Fmax(objNo)
                Fmax(objNo)=F(objNo);
            end
        end
    end    
end

%% membership function of each objective function
mu=zeros(maxOption,4);sumOptions=0;
for dmOption=1:maxOption
    for objNo=1:4
        if nfSolution(dmOption,objNo+6)<=Fmin(objNo)
            mu(dmOption,objNo)=1;
        elseif nfSolution(dmOption,objNo+6)>=Fmax(objNo)
            mu(dmOption,objNo)=0;
        else
            mu(dmOption,objNo)=(Fmax(objNo)-nfSolution(dmOption,objNo+6))/(Fmax(objNo)-Fmin(objNo));     
        end 
        sumOptions=sumOptions+mu(dmOption,objNo);
    end
end


%% membership function of solution
mfSolution=zeros(maxOption,1);bO=1;
for dmOption=1:maxOption
    sumObjectivesMF=0;
    for objNo=1:4
        sumObjectivesMF=sumObjectivesMF+mu(dmOption,objNo);        
    end
    mfSolution(dmOption)=sumObjectivesMF/sumOptions;
    if mfSolution(dmOption)>=mfSolution(bO)
        bO=dmOption;
    end
end

results=[nfSolution mu mfSolution]; 

%% Write results
disp('Now writing results...')
xlswrite('C:\Users\RAJAT\Documents\MATLAB\IOFiles\OPTio\PSO-MO.xlsx', results, 'Hourly', 'F4:V4');

%% Write best solution
hour=1;
bestResultW(hour,:)=[bO w1(bO) w2(bO) w3(bO) w4(bO) nfSolution(bO,1:12) mu(bO,1:4) mfSolution(bO)];
disp('Writing the best solution');
xlswrite('C:\Users\RAJAT\Documents\MATLAB\IOFiles\OPTio\PSO-MO.xlsx', bestResultW, 'BestH', 'A3:V3');







