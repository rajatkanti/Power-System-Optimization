%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program implements Differential Evolution in OCTAVE
% $Author: Dr. Rajat Kanti Samal$ $Date: 20-Jul-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$
% Reference: SWAYAM course on Computer Aided Applied Single Objective Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clc                            %To clear the command window
clear                          %To clear the workspace

disp('###################### DIFFERENTIAL EVOLUTION #########################')

%% problem setting
lb = [0 0];                     %lower bound          
ub = [10 10];                  %upper bound
prob = @SphereNew;             %fitness function


%% Parameters for Differential Evolution
Np = 5;                         %Population size
T = 100;                         %No. of iteration
Pc = 0.8;
F = 0.85;


%% Starting of DE
f = NaN(Np,1);

fu = NaN(Np,1);

D = length(lb);

U = NaN(Np,D);

P = repmat(lb,Np,1) + repmat((ub-lb),Np,1).*rand(Np,D);

for p = 1:Np
    f(p) = prob(P(p,:));
end

%% Iteration loop
for t = 1:T
    
    for i = 1:Np
        
        %% Mutation
        Candidates = [1:i-1 i+1:Np];
        idx = Candidates(randperm(Np-1,3));
        
        X1 = P(idx(1),:);
        X2 = P(idx(2),:);
        X3 = P(idx(3),:);
        
        V = X1 + F*(X2 - X3);
        
        
        %% Crossover
        
        del = randi(D,1);
        for j = 1:D
            
            if (rand <= Pc) || del == j
                U(i,j) = V(j);
            else
                U(i,j) = P(i,j);
            end
            
        end
    end
    
    %% Bounding and Greedy Selection
    for j = 1:Np
        
        U(j,:) = min(ub,U(j,:));
        U(j,:) = max(lb,U(j,:));
        
        fu(j) = prob(U(j,:));
        
        if fu(j) < f(j)
            P(j,:) = U(j,:);
            f(j) = fu(j);
        end
    end
    
    
end
  
  


 [bestfitness,ind] = min(f)
 bestsol = P(ind,:)
 
 
 
 
 
 
 
