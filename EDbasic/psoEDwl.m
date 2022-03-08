%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program solves economic dispatch problem with losses
% $Author: Dr. Rajat Kanti Samal$ $Date: 24-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$
%@Author: Rajat Kanti Samal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; 
clear;

example=2;

switch(example)
    case 1
        %% Kothari-Dhillon P-139, Example 3.1
        % Cost coefficients (a+b*Pg+c*Pg^2)
        a=[200 240];
        b=[10.333 10.833];
        c=[0.00889 0.00741];
        Pgmin=[0 0];
        Pgmax=[500 500];
        totPd=150; 
        ng=2; nb=2;
        alpha=0.05;
        B= [0.001 -0.0002
            -0.0002 0.002];
            
    case 2
        %% K. Uma Rao, P-327, Ex. 8.12
        % Cost coefficients (a+b*Pg+c*Pg^2)
        a=[320 200];
        b=[6.2 6.0];
        c=[0.004 0.003];
        Pgmin=[0 0];
        Pgmax=[500 500];
        totPd=412.35; 
        ng=2; nb=2; 
        alpha=0.001;
        B= (0.01)*[0.0125 0
            0      0.00625];
        
end

[Pg,lambda]=psoEDnoLossfn(ng,b,c,totPd);
[Pg2,lambda2,Ploss]=psoEDwithLossfn(totPd,ng,nb,b,c,B,Pg,lambda,alpha);
[totCost] = psoEDcostfn(ng,a,b,c,Pg2);

