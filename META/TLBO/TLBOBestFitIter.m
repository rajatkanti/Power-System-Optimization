clc                            %To clear the command window
clear                          %To clear the workspace

%% problem setting
lb = [0 0];                     %lower bound          
ub = [10 10];                  %upper bound
prob = @SphereNew;             %fitness function

%% Algorithm parameters
Np = 10;                         %Population size
T = 50;                         %No. of iteration

%% starting of TLBO
f = NaN(Np,1);
BestFitIter = NaN(T+1,1);

D = length(lb);

P = repmat(lb,Np,1) + repmat((ub-lb),Np,1).*rand(Np,D);

for p = 1:Np
  f(p) = prob(P(p,:));
end

BestFitIter(1) = min(f)

%% Iteration loop
for t = 1: T
    
    for i = 1:Np
        %% teacher phase
        Xmean = mean(P);             %Determining mean of the population
        
        [~,ind] = min(f);            %Determining the location of the teacher
        Xbest = P(ind,:);            %copying the solution acting as teacher
        
        TF = randi([1 2],1,1);       %Generating either 1 or 2 randomly for teaching factor
        
        Xnew = P(i,:) + rand(1,D).*(Xbest - TF*Xmean);   %Generating the new solution
        
        Xnew = min(ub, Xnew);        %Bounding the violating variables to their upper bound
        Xnew = max(lb, Xnew);        %Bounding the violating variables to their lower bound
        
        fnew = prob(Xnew);           %Evaluating the fitness of the newly generated solution
        
        if (fnew<f(i))               %Gready selection
            P(i,:) = Xnew;           %Include the new solution in population
            f(i) = fnew;             %Include the fitness function value of the new solution in population
        end
        
        
%% learner phase

p = randi([1 Np], 1,1);

%% Ensuring that the current member is not the partner
while i == p
    p = randi([1 Np],1,1);
end

if f(i)<f(p)
    Xnew = P(i,:) + rand(1, D).*(P(i,:) - P(p,:));
else
    Xnew = P(i,:) - rand(1, D).*(P(i,:) - P(p,:));
end

Xnew = min(ub, Xnew);
Xnew = min(lb, Xnew);

fnew = prob(Xnew);

if (fnew<f(i))
    P(i,:) = Xnew;
    f(i) = fnew;
end

    end
    
 BestFitIter(T+1) = min(f);
    
    disp(['Iteration' num2str(t) ': Best fitness = ' num2str(BestFitIter(t+1))])    
end

[bestfitness,ind] = min(f)
bestsol = P(ind,:)
