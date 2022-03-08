%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Economic Dispatch including losses
% The losses are computed by B-coefficients.
% Example 3.3., P-147, Power System Optimization by Kothari, Dhillon
% $Author: Dr. Rajat Kanti Samal$ $Date: 2022/02/01$    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$

clc;
clear;

%% Data; Quadratic cost function aPg2+bPg+c
NG = 2; % number of generators
a = [0.00889 0.00741];
b = [10.333 10.833]; 
c = [200 240]; 
Pd = 150; % MW
B = [0.001   -0.0002
     -0.0002   0.002]; Bi0=0;  
 
%% Initial solution (without losses)
sum1=0; sum2=0;
for i=1:NG
    sum1 = sum1 + b(i)/(2*a(i));
    sum2 = sum2 + 1/(2*a(i)); 
end

lambda = (Pd + sum1)/sum2;

Pg=zeros(NG,1); % ITERMAX number of rows, NG number of columns
for i=1:NG
    Pg(i) = (lambda - b(i))/(2*a(i)); 
end

Pg;

%% Iterative Process
ITERMAX=15; % Maximum number of iterations
stepSize = 0.05; tolerance=0.001; iterCount=zeros(ITERMAX,1); 
deltaP = zeros(ITERMAX,1); Ploss = zeros(ITERMAX,1); 
PgN = zeros(ITERMAX,NG); % ITERMAX number of rows and NG number of columns

% Initial generation
lambdaN=zeros(ITERMAX,1);
lambdaN(1)=lambda; % obtained from the lambda without losses
for i=1:NG
    PgN(1,i)=Pg(i); 
end
disp('Initial generation without considering losses')
PgN;
% ITERATIVE LOOP
for iter=1:ITERMAX
    iterCount(iter)=iter; 
    sum3=zeros(NG,1); % this contains sum_j=1^NG 2*B(i,j)*Pg(j)
    for i=1:NG
        for j=1:NG
            if j ~= i
                sum3(i) = sum3(i)+2*B(i,j)*PgN(iter,j);
            end
        end
    end
    
    for i=1:NG
        PgN(iter+1,i) = (lambdaN(iter)*(1-Bi0-sum3(i))-b(i))/(2*(a(i)+lambdaN(iter)*B(i,i)));
    end
    
    % Compute losses
    Ploss(iter)= psoBcoeffLfn(PgN(iter+1,:), NG, B);
    
    % Compute power imbalance
    deltaP(iter)=(Pd+Ploss(iter)) - sum(PgN(iter+1,:));
    
    % Modify lambda
    lambdaN(iter+1) = lambdaN(iter) + stepSize*(deltaP(iter));
    
    if abs(deltaP(iter)) < tolerance
        break; 
    end
    
end % end of for loop

%% display output
% PgN
% Ploss
% deltaP

%% Write to ouput file
results = [iterCount PgN lambdaN deltaP Ploss];
% disp('Now writing iterative process to output file...')
% xlswrite('C:\Users\RAJAT\Documents\MATLAB\Learning\PSO-DD\PSO-IO.xlsx', results, 'EDwithLoss', 'C7:H20');    



% Iterative loop










