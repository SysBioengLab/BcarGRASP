function [v,E1,E2] = r_MVD11(X,K)
% Metabolites definition 
A = X(1,:);
B = X(2,:);
P = X(3,:);
Q = X(4,:);
R = X(5,:);
S = X(6,:);
% Parameters definition K 
k01 = K(1);
k02 = K(2);
k03 = K(3);
k04 = K(4);
k05 = K(5);
k06 = K(6);
k07 = K(7);
k08 = K(8);
k09 = K(9);
k10 = K(10);
k11 = K(11);
k12 = K(12);
k13 = K(13);
k14 = K(14);
%  Numerator terms
E1 = k02*k13*k04*k06*k08.*P*k10.*Q+k02*k13*k04*k11*k06*k08.*P+k02*k13*k04*k11*k06*k09+k02*k13*k04*k11*k09*k07+k02*k13*k11*k09*k07*k05+k02*k04*k06*k08.*P*k10.*Q*k12.*R+k13*k11*k09*k07*k05*k03.*B;
E2 = k01.*A*k04*k13*k06*k08.*P*k10.*Q+k01.*A*k04*k13*k06*k11*k08.*P+k01.*A*k04*k13*k06*k11*k09+k01.*A*k04*k13*k11*k09*k07+k01.*A*k13*k11*k09*k07*k05+k01.*A*k04*k06*k08.*P*k10.*Q*k12.*R+k04*k06*k08.*P*k10.*Q*k12.*R*k14.*S;
E3 = k03.*B*k06*k01.*A*k08.*P*k13*k10.*Q+k03.*B*k06*k01.*A*k08.*P*k13*k11+k03.*B*k06*k01.*A*k13*k11*k09+k03.*B*k01.*A*k13*k11*k09*k07+k06*k08.*P*k10.*Q*k12.*R*k14.*S*k02+k03.*B*k06*k01.*A*k08.*P*k10.*Q*k12.*R+k03.*B*k06*k08.*P*k10.*Q*k12.*R*k14.*S;
E4 = k05*k08.*P*k03.*B*k10.*Q*k01.*A*k13+k05*k08.*P*k03.*B*k01.*A*k13*k11+k05*k03.*B*k01.*A*k13*k11*k09+k08.*P*k10.*Q*k12.*R*k14.*S*k02*k04+k05*k08.*P*k10.*Q*k12.*R*k14.*S*k02+k05*k08.*P*k03.*B*k10.*Q*k01.*A*k12.*R+k05*k08.*P*k03.*B*k10.*Q*k12.*R*k14.*S;
E5 = k07*k10.*Q*k05*k03.*B*k01.*A*k13+k07*k05*k03.*B*k01.*A*k13*k11+k10.*Q*k12.*R*k14.*S*k02*k04*k06+k07*k10.*Q*k12.*R*k14.*S*k02*k04+k07*k10.*Q*k05*k12.*R*k14.*S*k02+k07*k10.*Q*k05*k12.*R*k03.*B*k01.*A+k07*k10.*Q*k05*k12.*R*k03.*B*k14.*S;
E6 = k09*k07*k05*k03.*B*k01.*A*k13+k12.*R*k14.*S*k02*k04*k06*k08.*P+k09*k12.*R*k14.*S*k02*k04*k06+k09*k12.*R*k07*k14.*S*k02*k04+k09*k12.*R*k07*k14.*S*k05*k02+k09*k12.*R*k07*k05*k03.*B*k01.*A+k09*k12.*R*k07*k14.*S*k05*k03.*B;
E7 = k14.*S*k02*k04*k06*k08.*P*k10.*Q+k14.*S*k11*k02*k04*k06*k08.*P+k14.*S*k11*k02*k09*k04*k06+k14.*S*k11*k02*k09*k04*k07+k14.*S*k11*k02*k09*k07*k05+k11*k09*k07*k05*k03.*B*k01.*A+k14.*S*k11*k09*k07*k05*k03.*B;
% Denominator terms
Denominator = E1+E2+E3+E4+E5+E6+E7;
% Enzyme abundances terms
E1 = E1./Denominator;
E2 = E2./Denominator;
E3 = E3./Denominator;
E4 = E4./Denominator;
E5 = E5./Denominator;
E6 = E6./Denominator;
E7 = E7./Denominator;
% Reaction rate 
v = +k07.*E4-k08.*P.*E5;