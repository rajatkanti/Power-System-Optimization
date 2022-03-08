%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demonstration of use of fmincon for solution of NLP problems
% Find values of x that minimize f(x) = –x1x2x3, 
% starting at the point x = [10;10;10], subject to the constraints:
% 0 ? x1 + 2x2 + 2x3 ? 72.
% $Author: Dr. Rajat Kanti Samal$ $Date: 2022/01/25$    $Version: 1.0 $
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear; 

% Constraints written in standard <= form
A = [-1 -2 -2; ...
      1  2  2];
b = [0;72];

% Starting guess at the solution
x0 = [10;10;10];   

% Solve using fmincon
% x provides the variable values adn fval is the optimum function value
[x,fval] = fmincon(@psoNLPfn,x0,A,b);

disp('The opitmal values are ')
x
num2str(fval)