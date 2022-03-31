function [trial,P,fit,f] = GenNewSol(prob, lb, ub, Np, n, P, fit, trial, f, D)

j = randi(D,1)
p = randi(Np,1)

while (p == n)
    p = randi(Np,1)
end

Xnew = P(n,:)

phi = -1 + (1-(-1))*rand

Xnew(j) = P(n,j) + phi*(P(n,j) - P(p,j))
Xnew(j) = min(Xnew(j),ub(j))
Xnew(j) = max(Xnew(j),lb(j))

ObjNewSol = prob(Xnew)
FitnessNewSol = CalFit(ObjNewSol)

if (FitnessNewSol > fit(n))
    P(n,:) = Xnew
    fit(n) = FitnessNewSol
    f(n) = ObjNewSol
    trial(n) = 0
else
    trial(n) = trial(n)+1
end