function [Pg,lambda]=psoEDnoLossfn(ng,b,c,totPd)
% Economic generator scheduling without losses (iterative technique)
% Reads data from excel. Change load data and no of generators in program
% Generator limits included
% $Author: Dr. Rajat Kanti Samal$ $Date: 24-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$

Pd=totPd;
delP=Pd; 

Pg=zeros(1,ng);
lambda=1; %this can be changed to initial guess 
itermax=30; iter=0;

%% Solution using iterative method
% a,b,c may be cost coefficients or emission coefficients
tmp=0;
for i=1:ng
    tmp=tmp+(1/(2*c(i)));
end

while (abs(delP)>0.001 && iter<itermax)
    sumPg=0; 
    iter=iter+1;     
    delLAM=delP/tmp; %increment in lambda
    lambda=lambda+delLAM;
    for i=1:ng
        Pg(i)=(lambda-b(i))/(2*c(i));    
        sumPg=sumPg+Pg(i);
    end   
    delP=Pd-sumPg;      
end
             
        
end 


