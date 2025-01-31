function [v,E1,E2] = r_CrtIa1(X,K)
% Metabolites definition 
A = X(1,:);
B = X(2,:);
C = X(3,:);
P = X(4,:);
Q = X(5,:);
R = X(6,:);
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
k15 = K(15);
k16 = K(16);
k17 = K(17);
k18 = K(18);
%  Numerator terms
E1 = k02*k17*k04*k06.*P*k08*k10*k12.*Q*k14+k02*k17*k04*k15*k06.*P*k08*k10*k12.*Q+k02*k17*k04*k15*k06.*P*k13.*C*k08*k10+k02*k17*k04*k15*k06.*P*k13.*C*k08*k11+k02*k17*k04*k15*k06.*P*k13.*C*k11*k09+k02*k17*k04*k15*k13.*C*k11*k09*k07.*B+k02*k17*k15*k13.*C*k11*k09*k07.*B*k05+k02*k04*k06.*P*k08*k10*k12.*Q*k14*k16+k17*k15*k13.*C*k11*k09*k07.*B*k05*k03;
E2 = k01.*A*k04*k17*k06.*P*k08*k10*k12.*Q*k14+k01.*A*k04*k17*k06.*P*k15*k08*k10*k12.*Q+k01.*A*k04*k17*k06.*P*k15*k08*k13.*C*k10+k01.*A*k04*k17*k06.*P*k15*k08*k13.*C*k11+k01.*A*k04*k17*k06.*P*k15*k13.*C*k11*k09+k01.*A*k04*k17*k15*k13.*C*k11*k09*k07.*B+k01.*A*k17*k15*k13.*C*k11*k09*k07.*B*k05+k01.*A*k04*k06.*P*k08*k10*k12.*Q*k14*k16+k04*k06.*P*k08*k10*k12.*Q*k14*k16*k18.*R;
E3 = k03*k06.*P*k01.*A*k08*k17*k10*k12.*Q*k14+k03*k06.*P*k01.*A*k08*k17*k10*k15*k12.*Q+k03*k06.*P*k01.*A*k08*k17*k10*k15*k13.*C+k03*k06.*P*k01.*A*k08*k17*k15*k13.*C*k11+k03*k06.*P*k01.*A*k17*k15*k13.*C*k11*k09+k03*k01.*A*k17*k15*k13.*C*k11*k09*k07.*B+k06.*P*k08*k10*k12.*Q*k14*k16*k18.*R*k02+k03*k06.*P*k01.*A*k08*k10*k12.*Q*k14*k16+k03*k06.*P*k08*k10*k12.*Q*k14*k16*k18.*R;
E4 = k05*k08*k03*k10*k01.*A*k12.*Q*k17*k14+k05*k08*k03*k10*k01.*A*k12.*Q*k17*k15+k05*k08*k03*k10*k01.*A*k17*k15*k13.*C+k05*k08*k03*k01.*A*k17*k15*k13.*C*k11+k05*k03*k01.*A*k17*k15*k13.*C*k11*k09+k08*k10*k12.*Q*k14*k16*k18.*R*k02*k04+k05*k08*k10*k12.*Q*k14*k16*k18.*R*k02+k05*k08*k03*k10*k01.*A*k12.*Q*k14*k16+k05*k08*k03*k10*k12.*Q*k14*k16*k18.*R;
E5 = k07.*B*k10*k05*k12.*Q*k03*k14*k01.*A*k17+k07.*B*k10*k05*k12.*Q*k03*k01.*A*k17*k15+k07.*B*k10*k05*k03*k01.*A*k17*k15*k13.*C+k07.*B*k05*k03*k01.*A*k17*k15*k13.*C*k11+k10*k12.*Q*k14*k16*k18.*R*k02*k04*k06.*P+k07.*B*k10*k12.*Q*k14*k16*k18.*R*k02*k04+k07.*B*k10*k05*k12.*Q*k14*k16*k18.*R*k02+k07.*B*k10*k05*k12.*Q*k03*k14*k01.*A*k16+k07.*B*k10*k05*k12.*Q*k03*k14*k16*k18.*R;
E6 = k09*k12.*Q*k07.*B*k14*k05*k03*k01.*A*k17+k09*k12.*Q*k07.*B*k05*k03*k01.*A*k17*k15+k09*k07.*B*k05*k03*k01.*A*k17*k15*k13.*C+k12.*Q*k14*k16*k18.*R*k02*k04*k06.*P*k08+k09*k12.*Q*k14*k16*k18.*R*k02*k04*k06.*P+k09*k12.*Q*k07.*B*k14*k16*k18.*R*k02*k04+k09*k12.*Q*k07.*B*k14*k05*k16*k18.*R*k02+k09*k12.*Q*k07.*B*k14*k05*k16*k03*k01.*A+k09*k12.*Q*k07.*B*k14*k05*k16*k03*k18.*R;
E7 = k11*k14*k09*k07.*B*k05*k03*k01.*A*k17+k11*k09*k07.*B*k05*k03*k01.*A*k17*k15+k14*k16*k18.*R*k02*k04*k06.*P*k08*k10+k11*k14*k16*k18.*R*k02*k04*k06.*P*k08+k11*k14*k09*k16*k18.*R*k02*k04*k06.*P+k11*k14*k09*k16*k07.*B*k18.*R*k02*k04+k11*k14*k09*k16*k07.*B*k18.*R*k05*k02+k11*k14*k09*k16*k07.*B*k05*k03*k01.*A+k11*k14*k09*k16*k07.*B*k18.*R*k05*k03;
E8 = k13.*C*k11*k09*k07.*B*k05*k03*k01.*A*k17+k16*k18.*R*k02*k04*k06.*P*k08*k10*k12.*Q+k13.*C*k16*k18.*R*k02*k04*k06.*P*k08*k10+k13.*C*k16*k11*k18.*R*k02*k04*k06.*P*k08+k13.*C*k16*k11*k18.*R*k09*k02*k04*k06.*P+k13.*C*k16*k11*k18.*R*k09*k02*k07.*B*k04+k13.*C*k16*k11*k18.*R*k09*k02*k07.*B*k05+k13.*C*k16*k11*k09*k07.*B*k05*k03*k01.*A+k13.*C*k16*k11*k18.*R*k09*k07.*B*k05*k03;
E9 = k18.*R*k02*k04*k06.*P*k08*k10*k12.*Q*k14+k18.*R*k15*k02*k04*k06.*P*k08*k10*k12.*Q+k18.*R*k15*k02*k13.*C*k04*k06.*P*k08*k10+k18.*R*k15*k02*k13.*C*k04*k11*k06.*P*k08+k18.*R*k15*k02*k13.*C*k04*k11*k06.*P*k09+k18.*R*k15*k02*k13.*C*k04*k11*k09*k07.*B+k18.*R*k15*k02*k13.*C*k11*k09*k07.*B*k05+k15*k13.*C*k11*k09*k07.*B*k05*k03*k01.*A+k18.*R*k15*k13.*C*k11*k09*k07.*B*k05*k03;
% Denominator terms
Denominator = E1+E2+E3+E4+E5+E6+E7+E8+E9;
% Enzyme abundances terms
E1 = E1./Denominator;
E2 = E2./Denominator;
E3 = E3./Denominator;
E4 = E4./Denominator;
E5 = E5./Denominator;
E6 = E6./Denominator;
E7 = E7./Denominator;
E8 = E8./Denominator;
E9 = E9./Denominator;
% Reaction rate 
v = +k05.*E3-k06.*P.*E4;