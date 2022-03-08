function [ Ploss ] = psoBcoeffLfn( Pg,NG,B )
%computation of loss using B-coefficients
%   this function computes loss for economic dispatch
% $Author: Dr. Rajat Kanti Samal$ $Date: 2022/02/01$    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$


Ploss=0;

for i=1:NG
    for j=1:NG
        Ploss = Ploss + Pg(i)*B(i,j)*Pg(j);
    end
end

end

