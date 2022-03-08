% function [ output_args ] = psoNLPfn( input_args )
function f = psoNLPfn( x )
% objective function for psoNLP.m
%   demonstration of non-linear optimization using fmincon
% $Author: Dr. Rajat Kanti Samal$ $Date: 2022/01/25$    $Version: 1.0 $
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$

f = -x(1) * x(2) * x(3); % function equation

end

