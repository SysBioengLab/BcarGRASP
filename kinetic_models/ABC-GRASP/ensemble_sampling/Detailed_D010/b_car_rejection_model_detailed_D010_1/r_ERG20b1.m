function [v,E1,E2] = r_ERG20b1(X,K)
% Metabolites definition 
A = X(1,:);
B = X(2,:);
C = X(3,:);
D = X(4,:);
P1 = X(5,:);
P2 = X(6,:);
Q = X(7,:);
R = X(8,:);
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
k19 = K(19);
k20 = K(20);
%  Numerator terms
E1 = k02*k09*k12*k19*k04*k14*k06*k16+k02*k09*k12*k19*k04*k14*k17*k06+k02*k09*k12*k19*k04*k17*k06*k15+k02*k09*k12*k19*k04*k07*k14*k16+k02*k09*k12*k19*k04*k07*k14*k17+k02*k09*k12*k19*k04*k07*k17*k15+k02*k09*k12*k19*k07*k14*k05*k16+k02*k09*k12*k19*k07*k14*k17*k05+k02*k09*k12*k19*k07*k17*k05*k15+k02*k09*k12*k04*k14*k06*k16*k18.*P2+k02*k09*k12*k04*k07*k14*k16*k18.*P2+k02*k09*k12*k07*k14*k05*k16*k18.*P2+k02*k09*k19*k04*k17*k06*k15*k13.*D+k02*k09*k19*k04*k07*k17*k15*k13.*D+k02*k09*k19*k07*k17*k05*k15*k13.*D+k02*k12*k19*k04*k14*k06*k16*k08.*P1+k02*k12*k19*k04*k14*k17*k06*k08.*P1+k02*k12*k19*k04*k17*k06*k15*k08.*P1+k02*k12*k04*k14*k06*k16*k08.*P1*k18.*P2+k02*k19*k04*k17*k06*k15*k08.*P1*k13.*D+k09*k12*k19*k07*k14*k05*k16*k03.*B+k09*k12*k19*k07*k14*k17*k05*k03.*B+k09*k12*k19*k07*k17*k05*k15*k03.*B+k09*k12*k07*k14*k05*k16*k03.*B*k18.*P2+k09*k19*k07*k17*k05*k15*k03.*B*k13.*D;
E2 = k01.*A*k04*k09*k12*k19*k06*k14*k16+k01.*A*k04*k09*k12*k19*k06*k14*k17+k01.*A*k04*k09*k12*k19*k06*k17*k15+k01.*A*k04*k09*k12*k19*k07*k14*k16+k01.*A*k04*k09*k12*k19*k07*k14*k17+k01.*A*k04*k09*k12*k19*k07*k17*k15+k01.*A*k09*k12*k19*k07*k14*k05*k16+k01.*A*k09*k12*k19*k07*k14*k17*k05+k01.*A*k09*k12*k19*k07*k17*k05*k15+k01.*A*k04*k09*k12*k06*k14*k16*k18.*P2+k01.*A*k04*k09*k12*k07*k14*k16*k18.*P2+k01.*A*k09*k12*k07*k14*k05*k16*k18.*P2+k01.*A*k04*k09*k19*k06*k17*k15*k13.*D+k01.*A*k04*k09*k19*k07*k17*k15*k13.*D+k01.*A*k09*k19*k07*k17*k05*k15*k13.*D+k01.*A*k04*k12*k19*k06*k14*k08.*P1*k16+k01.*A*k04*k12*k19*k06*k14*k17*k08.*P1+k01.*A*k04*k12*k19*k06*k17*k08.*P1*k15+k01.*A*k04*k12*k06*k14*k08.*P1*k16*k18.*P2+k01.*A*k04*k19*k06*k17*k08.*P1*k15*k13.*D+k04*k06*k08.*P1*k10.*Q*k12*k19*k14*k16+k04*k06*k08.*P1*k10.*Q*k12*k19*k14*k17+k04*k06*k08.*P1*k10.*Q*k12*k19*k17*k15+k04*k06*k08.*P1*k10.*Q*k12*k14*k16*k18.*P2+k04*k06*k08.*P1*k10.*Q*k19*k17*k15*k13.*D;
E3 = k03.*B*k06*k01.*A*k09*k12*k19*k14*k16+k03.*B*k06*k01.*A*k09*k12*k19*k14*k17+k03.*B*k06*k01.*A*k09*k12*k19*k17*k15+k03.*B*k01.*A*k09*k12*k19*k07*k14*k16+k03.*B*k01.*A*k09*k12*k19*k07*k14*k17+k03.*B*k01.*A*k09*k12*k19*k07*k17*k15+k06*k08.*P1*k10.*Q*k02*k12*k19*k14*k16+k06*k08.*P1*k10.*Q*k02*k12*k19*k14*k17+k06*k08.*P1*k10.*Q*k02*k12*k19*k17*k15+k03.*B*k06*k01.*A*k09*k12*k14*k16*k18.*P2+k03.*B*k01.*A*k09*k12*k07*k14*k16*k18.*P2+k06*k08.*P1*k10.*Q*k02*k12*k14*k16*k18.*P2+k03.*B*k06*k01.*A*k09*k19*k17*k15*k13.*D+k03.*B*k01.*A*k09*k19*k07*k17*k15*k13.*D+k06*k08.*P1*k10.*Q*k02*k19*k17*k15*k13.*D+k03.*B*k06*k01.*A*k08.*P1*k12*k19*k14*k16+k03.*B*k06*k01.*A*k08.*P1*k12*k19*k14*k17+k03.*B*k06*k01.*A*k08.*P1*k12*k19*k17*k15+k03.*B*k06*k01.*A*k08.*P1*k12*k14*k16*k18.*P2+k03.*B*k06*k01.*A*k08.*P1*k19*k17*k15*k13.*D+k03.*B*k06*k08.*P1*k10.*Q*k12*k19*k14*k16+k03.*B*k06*k08.*P1*k10.*Q*k12*k19*k14*k17+k03.*B*k06*k08.*P1*k10.*Q*k12*k19*k17*k15+k03.*B*k06*k08.*P1*k10.*Q*k12*k14*k16*k18.*P2+k03.*B*k06*k08.*P1*k10.*Q*k19*k17*k15*k13.*D;
E4 = k05*k03.*B*k01.*A*k09*k12*k19*k14*k16+k05*k03.*B*k01.*A*k09*k12*k19*k14*k17+k05*k03.*B*k01.*A*k09*k12*k19*k17*k15+k08.*P1*k10.*Q*k02*k12*k19*k04*k14*k16+k08.*P1*k10.*Q*k02*k12*k19*k04*k14*k17+k08.*P1*k10.*Q*k02*k12*k19*k04*k17*k15+k05*k08.*P1*k10.*Q*k02*k12*k19*k14*k16+k05*k08.*P1*k10.*Q*k02*k12*k19*k14*k17+k05*k08.*P1*k10.*Q*k02*k12*k19*k17*k15+k05*k03.*B*k01.*A*k09*k12*k14*k16*k18.*P2+k08.*P1*k10.*Q*k02*k12*k04*k14*k16*k18.*P2+k05*k08.*P1*k10.*Q*k02*k12*k14*k16*k18.*P2+k05*k03.*B*k01.*A*k09*k19*k17*k15*k13.*D+k08.*P1*k10.*Q*k02*k19*k04*k17*k15*k13.*D+k05*k08.*P1*k10.*Q*k02*k19*k17*k15*k13.*D+k05*k08.*P1*k03.*B*k01.*A*k12*k19*k14*k16+k05*k08.*P1*k03.*B*k01.*A*k12*k19*k14*k17+k05*k08.*P1*k03.*B*k01.*A*k12*k19*k17*k15+k05*k08.*P1*k03.*B*k01.*A*k12*k14*k16*k18.*P2+k05*k08.*P1*k03.*B*k01.*A*k19*k17*k15*k13.*D+k05*k08.*P1*k03.*B*k10.*Q*k12*k19*k14*k16+k05*k08.*P1*k03.*B*k10.*Q*k12*k19*k14*k17+k05*k08.*P1*k03.*B*k10.*Q*k12*k19*k17*k15+k05*k08.*P1*k03.*B*k10.*Q*k12*k14*k16*k18.*P2+k05*k08.*P1*k03.*B*k10.*Q*k19*k17*k15*k13.*D;
E5 = k10.*Q*k02*k12*k19*k04*k14*k06*k16+k10.*Q*k02*k12*k19*k04*k14*k17*k06+k10.*Q*k02*k12*k19*k04*k17*k06*k15+k10.*Q*k07*k02*k12*k19*k04*k14*k16+k10.*Q*k07*k02*k12*k19*k04*k14*k17+k10.*Q*k07*k02*k12*k19*k04*k17*k15+k10.*Q*k07*k02*k12*k19*k05*k14*k16+k10.*Q*k07*k02*k12*k19*k05*k14*k17+k10.*Q*k07*k02*k12*k19*k05*k17*k15+k10.*Q*k02*k12*k04*k14*k06*k16*k18.*P2+k10.*Q*k07*k02*k12*k04*k14*k16*k18.*P2+k10.*Q*k07*k02*k12*k05*k14*k16*k18.*P2+k10.*Q*k02*k19*k04*k17*k06*k15*k13.*D+k10.*Q*k07*k02*k19*k04*k17*k15*k13.*D+k10.*Q*k07*k02*k19*k05*k17*k15*k13.*D+k07*k05*k03.*B*k01.*A*k12*k19*k14*k16+k07*k05*k03.*B*k01.*A*k12*k19*k14*k17+k07*k05*k03.*B*k01.*A*k12*k19*k17*k15+k07*k05*k03.*B*k01.*A*k12*k14*k16*k18.*P2+k07*k05*k03.*B*k01.*A*k19*k17*k15*k13.*D+k10.*Q*k07*k12*k19*k05*k14*k03.*B*k16+k10.*Q*k07*k12*k19*k05*k14*k17*k03.*B+k10.*Q*k07*k12*k19*k05*k17*k03.*B*k15+k10.*Q*k07*k12*k05*k14*k03.*B*k16*k18.*P2+k10.*Q*k07*k19*k05*k17*k03.*B*k15*k13.*D;
E6 = k11.*C*k14*k02*k09*k19*k16*k04*k06+k11.*C*k14*k02*k09*k19*k04*k17*k06+k11.*C*k02*k09*k19*k04*k17*k06*k15+k11.*C*k14*k02*k09*k19*k16*k04*k07+k11.*C*k14*k02*k09*k19*k04*k07*k17+k11.*C*k02*k09*k19*k04*k07*k17*k15+k11.*C*k14*k02*k09*k19*k16*k07*k05+k11.*C*k14*k02*k09*k19*k07*k17*k05+k11.*C*k02*k09*k19*k07*k17*k05*k15+k11.*C*k14*k02*k09*k16*k04*k18.*P2*k06+k11.*C*k14*k02*k09*k16*k04*k07*k18.*P2+k11.*C*k14*k02*k09*k16*k07*k18.*P2*k05+k14*k16*k18.*P2*k20.*R*k02*k09*k04*k06+k14*k16*k18.*P2*k20.*R*k02*k09*k04*k07+k14*k16*k18.*P2*k20.*R*k02*k09*k07*k05+k11.*C*k14*k02*k19*k16*k04*k06*k08.*P1+k11.*C*k14*k02*k19*k04*k17*k06*k08.*P1+k11.*C*k02*k19*k04*k17*k06*k15*k08.*P1+k11.*C*k14*k02*k16*k04*k18.*P2*k06*k08.*P1+k14*k16*k18.*P2*k20.*R*k02*k04*k06*k08.*P1+k11.*C*k14*k09*k19*k16*k07*k05*k03.*B+k11.*C*k14*k09*k19*k07*k17*k05*k03.*B+k11.*C*k09*k19*k07*k17*k05*k15*k03.*B+k11.*C*k14*k09*k16*k07*k18.*P2*k05*k03.*B+k14*k16*k18.*P2*k20.*R*k09*k07*k05*k03.*B;
E7 = k13.*D*k16*k11.*C*k02*k09*k19*k04*k06+k13.*D*k11.*C*k02*k09*k19*k04*k17*k06+k16*k18.*P2*k20.*R*k02*k09*k12*k04*k06+k13.*D*k16*k11.*C*k02*k09*k19*k04*k07+k13.*D*k11.*C*k02*k09*k19*k04*k07*k17+k16*k18.*P2*k20.*R*k02*k09*k12*k04*k07+k13.*D*k16*k11.*C*k02*k09*k19*k07*k05+k13.*D*k11.*C*k02*k09*k19*k07*k17*k05+k16*k18.*P2*k20.*R*k02*k09*k12*k07*k05+k13.*D*k16*k11.*C*k18.*P2*k02*k09*k04*k06+k13.*D*k16*k11.*C*k18.*P2*k02*k09*k04*k07+k13.*D*k16*k11.*C*k18.*P2*k02*k09*k07*k05+k13.*D*k16*k18.*P2*k20.*R*k02*k09*k04*k06+k13.*D*k16*k18.*P2*k20.*R*k02*k09*k04*k07+k13.*D*k16*k18.*P2*k20.*R*k02*k09*k07*k05+k13.*D*k16*k11.*C*k02*k19*k04*k06*k08.*P1+k13.*D*k11.*C*k02*k19*k04*k17*k06*k08.*P1+k16*k18.*P2*k20.*R*k02*k12*k04*k06*k08.*P1+k13.*D*k16*k11.*C*k18.*P2*k02*k04*k06*k08.*P1+k13.*D*k16*k18.*P2*k20.*R*k02*k04*k06*k08.*P1+k13.*D*k16*k11.*C*k09*k19*k07*k05*k03.*B+k13.*D*k11.*C*k09*k19*k07*k17*k05*k03.*B+k16*k18.*P2*k20.*R*k09*k12*k07*k05*k03.*B+k13.*D*k16*k11.*C*k18.*P2*k09*k07*k05*k03.*B+k13.*D*k16*k18.*P2*k20.*R*k09*k07*k05*k03.*B;
E8 = k15*k13.*D*k11.*C*k02*k09*k19*k04*k06+k18.*P2*k20.*R*k02*k09*k12*k04*k14*k06+k15*k18.*P2*k20.*R*k02*k09*k12*k04*k06+k15*k13.*D*k11.*C*k02*k09*k19*k04*k07+k18.*P2*k20.*R*k02*k09*k12*k04*k07*k14+k15*k18.*P2*k20.*R*k02*k09*k12*k04*k07+k15*k13.*D*k11.*C*k02*k09*k19*k07*k05+k18.*P2*k20.*R*k02*k09*k12*k07*k14*k05+k15*k18.*P2*k20.*R*k02*k09*k12*k07*k05+k15*k18.*P2*k13.*D*k11.*C*k02*k09*k04*k06+k15*k18.*P2*k13.*D*k11.*C*k02*k09*k04*k07+k15*k18.*P2*k13.*D*k11.*C*k02*k09*k07*k05+k15*k18.*P2*k13.*D*k20.*R*k02*k09*k04*k06+k15*k18.*P2*k13.*D*k20.*R*k02*k09*k04*k07+k15*k18.*P2*k13.*D*k20.*R*k02*k09*k07*k05+k15*k13.*D*k11.*C*k02*k19*k04*k06*k08.*P1+k18.*P2*k20.*R*k02*k12*k04*k14*k06*k08.*P1+k15*k18.*P2*k20.*R*k02*k12*k04*k06*k08.*P1+k15*k18.*P2*k13.*D*k11.*C*k02*k04*k06*k08.*P1+k15*k18.*P2*k13.*D*k20.*R*k02*k04*k06*k08.*P1+k15*k13.*D*k11.*C*k09*k19*k07*k05*k03.*B+k18.*P2*k20.*R*k09*k12*k07*k14*k05*k03.*B+k15*k18.*P2*k20.*R*k09*k12*k07*k05*k03.*B+k15*k18.*P2*k13.*D*k11.*C*k09*k07*k05*k03.*B+k15*k18.*P2*k13.*D*k20.*R*k09*k07*k05*k03.*B;
E9 = k20.*R*k02*k09*k12*k04*k14*k06*k16+k20.*R*k17*k02*k09*k12*k04*k14*k06+k20.*R*k17*k02*k09*k12*k15*k04*k06+k20.*R*k02*k09*k12*k04*k07*k14*k16+k20.*R*k17*k02*k09*k12*k04*k07*k14+k20.*R*k17*k02*k09*k12*k15*k04*k07+k20.*R*k02*k09*k12*k07*k14*k05*k16+k20.*R*k17*k02*k09*k12*k07*k14*k05+k20.*R*k17*k02*k09*k12*k15*k07*k05+k17*k15*k13.*D*k11.*C*k02*k09*k04*k06+k17*k15*k13.*D*k11.*C*k02*k09*k04*k07+k17*k15*k13.*D*k11.*C*k02*k09*k07*k05+k20.*R*k17*k02*k09*k15*k04*k13.*D*k06+k20.*R*k17*k02*k09*k15*k04*k07*k13.*D+k20.*R*k17*k02*k09*k15*k07*k13.*D*k05+k20.*R*k02*k12*k04*k14*k06*k16*k08.*P1+k20.*R*k17*k02*k12*k04*k14*k06*k08.*P1+k20.*R*k17*k02*k12*k15*k04*k06*k08.*P1+k17*k15*k13.*D*k11.*C*k02*k04*k06*k08.*P1+k20.*R*k17*k02*k15*k04*k13.*D*k06*k08.*P1+k20.*R*k09*k12*k07*k14*k05*k16*k03.*B+k20.*R*k17*k09*k12*k07*k14*k05*k03.*B+k20.*R*k17*k09*k12*k15*k07*k05*k03.*B+k17*k15*k13.*D*k11.*C*k09*k07*k05*k03.*B+k20.*R*k17*k09*k15*k07*k13.*D*k05*k03.*B;
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
v = +k17.*E8-k18.*P2.*E9;