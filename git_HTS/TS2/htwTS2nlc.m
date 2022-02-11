function [cnl,ceqnl]=htwTS2nlc(PgH)
% non linear constraint (Method 2) for hydrotherma TS-2 
%   Refer the main function for details
% $Author: Dr. Rajat Kanti Samal$ $Date: 11-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$


NG1=2;M=24;
Wj=[25 35];
totPdTS2=[30	33	35	38	40	45	50	59	61	58	56	57	60	61	65	68	71	62	55	50	43	33	31	30];

%a+bPgi+cPgi2
ah=[0.2 0.4];
bh=[0.03 0.06];
ch=[0.00005 0.0001];

q=zeros(M,2); cnl=zeros(M,1);
for m=1:M  
    cnl(m)=sum(PgH(m,:))-totPdTS2(m);
    for g=1:NG1
        q(m,g)=ah(g)+bh(g)*PgH(m,g)+ch(g)*(PgH(m,g)^2);
    end
end

Q=zeros(NG1,1);ceqnl=zeros(NG1,1);
for g=1:NG1
    Q(g)=sum(q(:,g));
    ceqnl(g)=Q(g)-Wj(g);
end

