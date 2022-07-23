%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program implements Genetic Algorithm in OCTAVE
% $Author: Dr. Rajat Kanti Samal$ $Date: 20-Jul-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$
% Reference: SWAYAM course on Computer Aided Applied Single Objective Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;                    % To clear the command window
clear;                  % To clear the workspace

disp('################### GENETIC ALGORITHM ###################')

% rng(2,'twister')        % Fixing the random number generator and seed

%% Problem Settings
lb = [0 0];  % Lower Bound
ub = [10 10]; % Upper Bound
prob = @SphereNew; % Fitness Function

%% Algorithm Parameters
Np = 6; % Population Size
T = 50; % Number of iterations
etac = 20; % Distribution index for crossover
etam = 20; % Distribution index for mutation
Pc = 0.8; % Crossover probability
Pm = 0.2; % Mutation probability

%% Genetic Algorithm
f = NaN(Np,1); % Vector to store the fitness function
OffspringObj = NaN(Np,1); % Vector to store fitness function of the offspring
D = length(lb);  % Number of decision variables
P = repmat(lb,Np,1) + repmat((ub-lb),Np,1).*rand(Np,D); % Generation of initial population

for p = 1:Np
    f(p) = prob(P(p,:));
end

%% Iteration Loop
for t=1:T
    
    %% Tournament Selection
    MatingPool = TournamentSelection(f,Np); 
    Parent = P(MatingPool,:); 
    
    %% Crossover
    offspring = CrossoverSBX(Parent,Pc,etac,lb,ub); 
    
    %% Mutation
    offspring = PolyMutation(offspring,Pm,etam,lb,ub);
    
    for i=1:Np
        OffspringObj(i) = prob(offspring(i,:)); 
    end
    
    CombinedPopulation = [P; offspring];
    [f,ind] = sort([f; OffspringObj]); % mu+lambda selection
    
    f = f(1:Np);
    P = CombinedPopulation(ind(1:Np),:);    
    
end

[bestfitness, ind] = min(f) 
bestsol = P(ind,:) 




