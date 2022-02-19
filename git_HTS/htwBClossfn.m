function [ Ploss ] = htwBClossfn( NG,Pg,B )
%Calculation of loss using B coefficients
%   This function calculates loss using B coefficients.



Ploss=0;
for i=1:NG
    for j=1:NG
        Ploss=Ploss+Pg(i)*B(i,j)*Pg(j);
    end
end



end% End of function

