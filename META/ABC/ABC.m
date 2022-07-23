%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program implements Artificial Bee Colony Optimization in OCTAVE
% $Author: Dr. Rajat Kanti Samal$ $Date: 20-Jul-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$
% Reference: SWAYAM course on Computer Aided Applied Single Objective Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
clear


disp('################# ARTIFICIAL BEE COLONY OPTIMIZATION #################')


%% Problem settings
lb = [0 0];
ub = [10 10];
prob = @SphereNew;

%% Algorithm parameters
Np = 10;
T = 50;
limit = 5;

%% Starting of ABC
f = NaN(Np,1);
fit = NaN(Np,1);
trial = NaN(Np,1);

D = length(lb);

P = repmat(lb,Np,1) + repmat((ub-lb),Np,1).*rand(Np,D);

for p = 1:Np
    f(p) = prob(P(p,:));
    fit(p) = CalFit(f(p));
end

[bestobj, ind] = min(f);
bestsol = P(ind,:);

for t = 1:T
    
    %% Employed Bee Phase
    for i = 1:Np
    [trial,P,fit,f] = GenNewSol(prob, lb, ub, Np, i, P, fit, trial, f, D);
    end
    
    %% onlooker Bee Phase
    probability = 0.9 * (fit/max(fit)) + 0.1;
    
    m = 0; n = 1;
    
    while(m < Np)
        n;
        prob(n);
        if(rand < prob(n))
            [trial,P,fit,f] = GenNewSol(prob, lb, ub, Np, n, P, fit, trial, f, D);
            m = m + 1;
        end
        n = mod(n,Np) + 1;
    end
    
    [bestobj,ind] = min([f;bestobj]);
    CombinedSol = [P;bestsol];
    bestsol = CombinedSol(ind,:);
    
     %% Scout Bee Phase
     [val,ind] = max(trial);
     
     if (val > limit)
         trial(ind) = 0;
         P(ind,:) = lb + (ub-lb).*rand(1,D);
         f(ind) = prob(P(ind,:));
         fit(ind) = CalFit(f(ind));
     end
end

[bestobj,ind] = min([f;bestobj])
CominedSol = [P ;bestsol];
bestsol = CombinedSol(ind,:)
         
    
