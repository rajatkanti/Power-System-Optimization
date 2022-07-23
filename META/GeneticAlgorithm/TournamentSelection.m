function MatingPool = TournamentSelection( f, Np )
%Tournament Selection
%   Allows each solution to participate exactly twice
% $Author: Dr. Rajat Kanti Samal$ $Date: 20-Jul-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$
% Reference: SWAYAM course on Computer Aided Applied Single Objective Optimization

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

