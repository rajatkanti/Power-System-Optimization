function [cnl,ceqnl]=htwTS1nlc(PgH)
% hydrotherma TS-1 non linear constraint (Method 2)
% Refer main function for test system details
% $Author: Dr. Rajat Kanti Samal$ $Date: 11-Feb-2022 $    $Version: 1.0$
% $Veer Surendra Sai University of Technology, Burla, Odisha, India$

NG2=1;
Wj=2559.6;

%a+bPgi+cPgi2
ah=61.53;
bh=-0.009079;
ch=0.0007749;

q=zeros(24,1); 
for m=1:24         
    q(m)=ah+bh*PgH(m)+ch*(PgH(m)^2);
end

Q=sum(q);


cnl=[];
ceqnl=Q-Wj;