%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution of Linear Programming problem using linprog of MATLAB
% Optmization toolbox
% Example 2.1-1 (Reddy Mikks compnay) in Operations Research (Taha)
% $Author: Dr. Rajat Kanti Samal$ $Date: 2022/01/25 18:23:52 $    $Version: 1.0 $
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;

% Objective function
f = [-5 -4]; % Since the problem is of maximization, we multiply by -1

% Inequality Constraints
A = [6 4; 1 2; -1 1; 0 1]; b=[24 6 1 2]'; 

% Solution 
x = linprog(f,A,b)
profit = -f*x