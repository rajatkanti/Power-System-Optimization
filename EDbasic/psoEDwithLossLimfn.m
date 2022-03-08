function [Pg,lambda,Ploss]=psoEDwithLossLimfn(Pgmin,Pgmax,Pd,nb,ng,b,c,B,Pg,lambda)
% Economic generator scheduling with losses (iterative technique)
% Reads data from excel. Change load data and no of generators in program
% Generator limits included
% $Author: Dr. Rajat Kanti Samal$ $Date: 24-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$

delP=Pd; 
n=ng; %number of generators

itermax=50; iter=0;iter1=0;iter2=0; Bio=0;

%% Solution using iterative technique
% a,b,c may be cost coefficients or emission coefficients
tmp=0;
for i=1:n
    tmp=tmp+(1/(2*c(i)));
end

while (abs(delP)>0.001 && iter<itermax)
    sumPg=0; 
    iter=iter+1;     
    delLAM=delP/tmp; %increment in lambda
    lambda=lambda+delLAM;
    for i=1:n
    Pg(i)=(lambda-b(i))/(2*c(i));    
    sumPg=sumPg+Pg(i);
    end   
    delP=Pd-sumPg;      
end

%% Incorporate generator limits 
fixed=zeros(n,1);
sumPg=0;nf=1; 
for i=1:n    
    if Pg(i) < Pgmin(i)
            fixed(nf)=i;
            nf=nf+1;
            Pg(i)=Pgmin(i);
        else if Pg(i)>Pgmax(i)
                Pg(i)=Pgmax(i);
                fixed(nf)=i;
                nf=nf+1;
            end
    end
    sumPg=sumPg+Pg(i);
end
% fixed %to see for which generators constraint is violated

delP=Pd-sumPg; 
if abs(delP)>0.001 %only if constraints are violated
    ouput('Generator Limits Violated')
    while (abs(delP)> 0.0001 && iter1<itermax)
        sumPg=0; 
        iter1=iter1+1;     
        delLAM=delP/tmp; %increment in lambda
        lambda=lambda+delLAM;
        for i=1:n
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
    F=zeros(1,n); Tcost=0;
    for i=1:n
        F(i)=a(i)+b(i)*Pg(i)+c(i)*(Pg(i)^2);
        Tcost=Tcost+F(i);
    end

    
end %End of generator limit violation if loop

%% Incorporate lossess
Ploss=0;
for i=1:nb
    for j=1:nb
        Ploss=Ploss+B(i,j)*Pg(i)*Pg(j);
    end
end
Ploss;
sumPg=0;
for i=1:n
    sumPg=sumPg+Pg(i);
end
delP=Pd+Ploss-sumPg;
alpha=0.001;

%fprintf('Ite Pg1(MW)   Pg2(MW) Pl(MW) lambda(Rs/MWh)  delP(MW)\n');

while (abs(delP)>0.001 && iter2<itermax)
    iter2=iter2+1;
    sumPg=0;
    sumBP=0;    
    % dPgdLam=0;
    
    for i=1:n
        for j=1:n
            if (i~=j)
            sumBP=sumBP+2*B(i,j)*Pg(j);
            end
        end
        Pg(i)=(lambda*(1-Bio-sumBP)-b(i))/(2*(c(i)+lambda*B(i,i)));        
        % dPgdLam=dPgdLam+(c(i)+b(i)*B(i,i))/(2*(c(i)+lambda*B(i,i))^2);
    end
    Ploss=0;
    for i=1:nb
        for j=1:nb
            Ploss=Ploss+B(i,j)*Pg(i)*Pg(j);
        end
    end    
    
    for i=1:n
        sumPg=sumPg+Pg(i);
    end
    delP=Pd+Ploss-sumPg;
    %delLAM=delP/dPgdLam;
    %lambda=lambda+delLAM;
    lambda=lambda+alpha*delP;    
end
% End of Loss inclusion

Ploss=0;
for i=1:nb
    for j=1:nb
        Ploss=Ploss+B(i,j)*Pg(i)*Pg(j);
    end
end



