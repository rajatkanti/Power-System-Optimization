%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demonstration of use of fmincon for solution of NLP problems
% Exercise 7.16 of S. S. Rao, Engineering Optimization
% $Author: Dr. Rajat Kanti Samal$ $Date: 2022/01/25$    $Version: 1.0 $
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear; 

% Constraints written in standard <= form
A = [1 2; 4 3; 6 1];
b = [5; 10; 7];

% Starting guess at the solution
x0 = [1;1];   

% Solve using fmincon
% x provides the variable values adn fval is the optimum function value
[x,fval] = fmincon(@psoNLP2fn,x0,A,b);

disp('The opitmal values are ')
x
num2str(fval)