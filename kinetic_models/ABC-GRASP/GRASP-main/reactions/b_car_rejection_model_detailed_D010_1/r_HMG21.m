function v = r_HMG21(X,negEff,Kr,KnegEff,L,n) 
% Parameters definition 
[vR,eR] = r_HMG21Catalytic(X,Kr); 
Q = L*eR.^n; 
KnegEff = KnegEff(ones(size(negEff,2),1),:); 
Q = Q.*((1 + sum(negEff'./KnegEff,2)).^n)'; 
% Reaction rate 
v = n*vR./(1 + Q);