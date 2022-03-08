%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program solves economic dispatch problem without losses
% $Author: Dr. Rajat Kanti Samal$ $Date: 24-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$
%@Author: Rajat Kanti Samal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; 
clear;

example=4;

switch(example)
    case 1
        %% Kothari-Dhillon P-139, Example 3.1
        % Cost coefficients (a+b*Pg+c*Pg^2)
        a=[120 120];
        b=[22 16];
        c=[0.05 0.06];
        Pgmin=[20 20];
        Pgmax=[100 100];
        totPd=180; 
        ng=2; 
        
    case 2
        %% Kothari-Dhillon, P-174, Example 3.3
        % Cost coefficients (a+b*Pg+c*Pg^2)
        a=[200 240];
        b=[10.333 10.833];
        c=[0.00889 0.00741];
        Pgmin=[0 0];
        Pgmax=[500 500];
        totPd=150; 
        ng=2; 
        
    case 3
        %% K. Uma Rao, P-310, Ex. 8.1
        % Cost coefficients (a+b*Pg+c*Pg^2)
        a=[1.5 1.9];
        b=[20 30];
        c=[0.1 0.1];
        Pgmin=[0 0];
        Pgmax=[500 500];
        totPd=200; 
        ng=2; 
        
    case 4
        %% K. Uma Rao, P-312, Ex. 8.3
        % Cost coefficients (a+b*Pg+c*Pg^2)
        a=[350 500 600];
        b=[7.2 7.3 6.74];
        c=[0.004 0.0025 0.003];
        Pgmin=[0 0];
        Pgmax=[500 500];
        totPd=800; 
        ng=3; 
        
end

[Pg,lambda]=optEDnoLossfn(ng,b,c,totPd)
[F,totCost] = optEDcostfn(ng,a,b,c,Pg)

