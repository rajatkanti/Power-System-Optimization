function MatingPool = TournamentSelection( f, Np )
%Tournament Selection
%   Allows each solution to participate exactly twice

MatingPool = NaN(Np,1); 
indx = randperm(Np); % Randomly suffling the index of population


for i = 1:Np-1                              % Poolsize is Np
    Candidate = [indx(i) indx(i+1)]; % Selecting one pair of population
    CandidateObj = f(Candidate);
    [~, ind] = min(CandidateObj); % Selecting winner
    MatingPool(i) = Candidate(ind);  % Selecting index of winner  
end

% Tounament selection between last and first member
Candidate = [indx(Np) indx(1)];
CandidateObj = f(Candidate);
[~, ind] = min(CandidateObj);
MatingPool(Np) = Candidate(ind);

