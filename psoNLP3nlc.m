function [ nlc1, nlc2 ] = psoNLP3nlc( x )
% non-linear constraints for psoNLP3.m
%   demonstration of handling non-linear constraints using fmincon
% $Author: Dr. Rajat Kanti Samal$ $Date: 2022/01/25$    $Version: 1.0 $
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$


nlc1=-x(1)^2 + x(2) -4;
nlc2 = -(x(1)-2)^2 + x(2)-3; 

end

