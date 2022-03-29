function [Pg,lambda,Ploss]=psoEDwithLoss2fn(Pd,ng,nb,b,c,B,Pg,lambda)
% Economic generator scheduling with losses (iterative technique)
% It implements algorithm of P-327 in K. Uma Rao
% Generator limit checks not incorporated
% $Author: Dr. Rajat Kanti Samal$ $Date: 24-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$

n=ng; %number of generators

itermax=100; iter=0; Bio=0; 

%% Incorporate lossess
Ploss=psoBcoeffLfn( Pg,ng,B );
sumPg=0; 
for i=1:ng
    sumPg=sumPg+Pg(i);
end
delP=Pd+Ploss-sumPg;
% alpha=0.001;
% Step (i): initial estimate of lambda is provided in function call

while (abs(delP)>0.001 && iter<itermax)
    iter=iter+1;
    sumPg=0;
    sumBP=0;    
    
    %% Obtain PGi(r)using Eq. 8.42, P-326 KUR (Step (ii))
    dPgdLamSum =0; 
    for i=1:ng
        for j=1:ng
            if (i~=j)
                sumBP=sumBP+2*B(i,j)*Pg(j);
            end
        end
        sumBP2=2*c(i)*sumBP;
        % Eq. 8.42, P-326, KUR
        Pg(i)=(lambda*(1-Bio-sumBP)-b(i))/(2*(c(i)+lambda*B(i,i))); 
        % Eq. 8.45, P-326, KUR
        dPgdLam=(c(i)+b(i)*B(i,i)-sumBP2)/(2*(c(i)+lambda*B(i,i))^2); 
        dPgdLamSum = dPgdLamSum + dPgdLam; 
    end
    %% Compute losses (Step (iii))
    Ploss=psoBcoeffLfn( Pg,nb,B ); 
    
    %% Compute DeltaP(r) from Eq. 8.46 (Step (v))
    for i=1:n
        sumPg=sumPg+Pg(i);
    end
    delP=Pd+Ploss-sumPg;
    
    %% STEP-VI: Compute increment in lambda (Eq. 8.47)
    delLAM=delP/dPgdLamSum;
    % delLAM
    lambda=lambda+delLAM;
    %lambda=lambda+alpha*delP;    
end
% End of Loss inclusion
%iter
%delP
Ploss=psoBcoeffLfn( Pg,ng,B );   



