function [Pg,lambda,Ploss]=psoEDwithLossfn(Pd,ng,nb,b,c,B,Pg,lambda,alpha)
% Economic generator scheduling with losses (iterative technique)
% It implements algorithm of Eq. 8.42, P-326 in K. Uma Rao
% and also the algorithm of Eq. 3.21, P-145, Kothari-Dhillon
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

while (abs(delP)>0.001 && iter<itermax)
    iter=iter+1;
    sumPg=0;
    sumBP=0;    
    % dPgdLam=0;    
    for i=1:ng
        for j=1:ng
            if (i~=j)
                sumBP=sumBP+2*B(i,j)*Pg(j);
            end
        end
        Pg(i)=(lambda*(1-Bio-sumBP)-b(i))/(2*(c(i)+lambda*B(i,i)));     
        % dPgdLam=dPgdLam+(c(i)+b(i)*B(i,i))/(2*(c(i)+lambda*B(i,i))^2);
    end
    Ploss=psoBcoeffLfn( Pg,nb,B );    
    for i=1:n
        sumPg=sumPg+Pg(i);
    end
    delP=Pd+Ploss-sumPg;
    %delLAM=delP/dPgdLam;
    %lambda=lambda+delLAM;
    lambda=lambda+alpha*delP;    
end
% End of Loss inclusion
%iter
%delP
Ploss=psoBcoeffLfn( Pg,ng,B );   



