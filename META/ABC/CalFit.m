function fit = CalFit(f)
% This function computes the fitness value for 
% Artificial Bee Colony Optimization in OCTAVE
% $Author: Dr. Rajat Kanti Samal$ $Date: 20-Jul-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$
% Reference: SWAYAM course on Computer Aided Applied Single Objective Optimization

if f >= 0
    fit = 1 / (1+f);
else
    fit = 1 + abs(f)
end