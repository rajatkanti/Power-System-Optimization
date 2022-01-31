%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demonstration of use of fmincon for solution of NLP problems
% Exercise 7.15 of S. S. Rao, Engineering Optimization
% $Author: Dr. Rajat Kanti Samal$ $Date: 2022/01/25$    $Version: 1.0 $
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear; 

% Constraints written in standard <= form
A = []; b = [];
Aeq = []; beq = [];
lb=[-Inf; -Inf]; ub=[Inf, Inf];

% Starting guess at the solution
x0 = [1;1];   


% Solve using fmincon
% x provides the variable values adn fval is the optimum function value
% x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon)
[x,fval] = fmincon(@psoNLP3fn,x0,A,b,Aeq,beq,lb,ub,@psoNLP3nlc);

disp('The opitmal values are ')
x
num2str(fval)