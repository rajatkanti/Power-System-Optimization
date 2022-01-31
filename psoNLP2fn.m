function f = psoNLP2fn( x )
% objective function for fmincon
%   Exercise 7.16 of S. S. Rao, Engineering Optimization
% $Author: Dr. Rajat Kanti Samal$ $Date: 2022/01/25$    $Version: 1.0 $
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$


f = (x(1)-1)^2 + (x(2)-2)^2 - 4; % function equation

end

