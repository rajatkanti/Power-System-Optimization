function offspring = CrossoverSBX( Parent,Pc,etac,lb,ub )
%Crossover function for genetic algorithm
%   SBX crossover is performed

[Np, D] = size(Parent); % No of population and decision variables
indx = randperm(Np); % Permutating numbers from 1 to Np
Parent = Parent(indx,:); % Randomly suffling parent solutions
offspring = NaN(Np,D); % Matrix to store offspring solutions

for i=1:2:Np % Selecting parents in pairs for crossover
    
    r = rand; % To check if crossover is to be performed
    if r < Pc
        
        for j=1:D
            
            r = rand;
            
            if r <= 0.5
                beta = (2*r)^(1/(etac+1));
            else
                beta = (1/(2*(1-r)))^(1/(etac+1)); 
            end
            
            offspring(i,j) = 0.5*(((1+beta)*Parent(i,j))+((1-beta)*Parent(i+1,j)));
            offspring(i+1,j) = 0.5*(((1-beta)*Parent(i,j))+((1+beta)*Parent(i+1,j)));
            
        end
        
        offspring(i,:) = max(offspring(i,:),lb);
        offspring(i+1,:) = min(offspring(i+1,:),ub);
        
    else
        
        offspring(i,:) = Parent(i,:);
        offspring(i+1,:) = Parent(i+1,:);
        
        
    end
    
    
    
end

end

