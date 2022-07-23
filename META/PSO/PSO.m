%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program implements Particle Swarm Optimization in OCTAVE
% $Author: Dr. Rajat Kanti Samal$ $Date: 20-Jul-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$
% Reference: SWAYAM course on Computer Aided Applied Single Objective Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
clear

disp('################### PARTICLE SWARM OPTIMIZATION ###################')

%% problem settings
lb = [0 0];
ub = [10 10];
prob = @SphereNew;

%% Algorithm parameters
Np = 10;
T = 50;
w = 0.8;
c1 = 1.5;
c2 = 1.5;

 
%% Particle Swarm Optimization
f = NaN(Np,1);

D = length(lb);

P = repmat(lb,Np,1) + repmat((ub-lb),Np,1).*rand(Np,D);
v = repmat(lb,Np,1) + repmat((ub-lb),Np,1).*rand(Np,D);

for p = 1:Np
    f(p) = prob(P(p,:));
end

pbest = P;
f_pbest = f;


[f_gbest,ind] = min(f_pbest);
gbest = P(ind,:);

for t = 1:T
    
    for p = 1:Np
        
        v(p,:) = w*v(p,:) + c1*rand(1,D).*(pbest(p,:)-P(p,:)) + c2*rand(1,D).*(gbest - P(p,:));
       
        P(p,:) = P(p,:) + v(p,:);
        
        P(p,:) = max(P(p,:),lb);
        P(p,:) = min(P(p,:),ub);
        
        f(p) = prob(P(p,:));
        
        if f(p) < f_pbest(p)
            
            f_pbest(p) = f(p);
            pbest(p,:) = P(p,:);
            
            if f_pbest(p) < f_gbest
                
                f_gbest= f_pbest(p);
                gbest = pbest(p,:);
                
            end
            
        end
    end
end

bestfitness = f_gbest
bestsol = gbest
            
        
