% function [ output_args ] = psoNLPfn( input_args )
function f = psoNLP3fn( x )
% objective function for fmincon
%   main program is psoNLP3.m
% $Author: Dr. Rajat Kanti Samal$ $Date: 2022/01/25$    $Version: 1.0 $
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$

f = (x(1)-1)^2 + (x(2)-5)^2; % function equation

end

