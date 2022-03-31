function F = Rastrigin( x )
%Rastrigin function
%   for evaluating optimization techniques

F = sum(x.^2 - 10.*cos(2.*pi.*x) + 10);


end

