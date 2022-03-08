function [Pg,lambda]=psoEDnoLossLimfn(ng,a,b,c,Pgmin,Pgmax,totPd)
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

%% Incorporate generator limits
fixed=zeros(ng,1);
sumPg=0;nf=1; 
for i=1:ng    
    if Pg(i) < Pgmin(i)
            fixed(nf)=i;
            nf=nf+1;
            Pg(i)=Pgmin(i);
    elseif Pg(i)>Pgmax(i)
                Pg(i)=Pgmax(i);
                fixed(nf)=i;
                nf=nf+1;
    end    
    sumPg=sumPg+Pg(i);
end
delP=Pd-sumPg; 
%% Allocate remaining load in non-limit-violated generators
if abs(delP)>0.001 %only if constraints are violated     
    
    while (abs(delP)> 0.0001 && iter<itermax)
        sumPg=0; 
        iter=iter+1;     
        delLAM=delP/tmp; %increment in lambda
        lambda=lambda+delLAM;        
        for i=1:ng
            %fixed power generators no longer participate in iteration,
            %however, their power is included in sumPg
            %What if limit violation in this iteration also??????
            if (isempty(find(fixed==i, 1)))%no i in the array 'fixed'
                Pg(i)=(lambda-b(i))/(2*c(i));               
            end
            sumPg=sumPg+Pg(i);
        end   
        delP=Pd-sumPg; 
    end    
    
end % End of 'if abs(delP)>0.001 loop'
             
        
end 


