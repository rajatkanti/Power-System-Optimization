function offspring = PolyMutation( offspring,Pm,etam,lb,ub )
% Polynomial Mutation for Genetic Algorihtm
%   This function is used for mutation
% $Author: Dr. Rajat Kanti Samal$ $Date: 21-Jul-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$
% Reference: SWAYAM course on Computer Aided Applied Single Objective Optimization

[Np, D] = size(offspring);

for i=1:Np
    r = rand; 
    if r<Pm
        for j=1:D
            r = rand;
            if r<0.5
                delta = (2*r)^(1/(etam+1)) -1; 
            else
                delta = (1 - 2*(1-r))^(1/(etam+1)) - 1; 
            end
            offspring(i,j) = offspring(i,j) + (ub(j) - lb(j))*delta; 
        end
        offspring(i,:) = max(offspring(i,:),lb); 
        offspring(i,:) = min(offspring(i,:),ub); 
    end
end


end

