%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution of Linear Programming problem using linprog of MATLAB
% Optmization toolbox
% Example 3.4, P-156, Engineering Optimization, S. S. Rao
% $Author: Dr. Rajat Kanti Samal$ $Date: 2022/01/25$    $Version: 1.0 $
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;

% Objective function
f = [-1 -2 -1]; 

% Inequality Constraints
A = [2 1 -1; 2 -1 5; 4 1 1]; b=[2 6 6]' ;

% Solution 
x = linprog(f,A,b)
profit = f*x