function [v,E1,E2] = r_ERG131(X,K)
% Metabolites definition 
A = X(1,:);
B = X(2,:);
I1 = X(3,:);
I2 = X(4,:);
P = X(5,:);
Q = X(6,:);
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
%  Numerator terms
E1 = k02*k11*k14*k16*k04*k06.*P*k08+k02*k11*k14*k16*k04*k09*k06.*P+k02*k11*k14*k16*k04*k09*k07.*B+k02*k11*k14*k16*k09*k07.*B*k05+k02*k14*k16*k04*k06.*P*k08*k10+k11*k14*k16*k09*k07.*B*k05*k03;
E2 = k01.*A*k04*k11*k14*k16*k06.*P*k08+k01.*A*k04*k11*k14*k16*k06.*P*k09+k01.*A*k04*k11*k14*k16*k09*k07.*B+k01.*A*k11*k14*k16*k09*k07.*B*k05+k01.*A*k04*k14*k16*k06.*P*k08*k10+k04*k06.*P*k08*k10*k12.*Q*k14*k16;
E3 = k03*k06.*P*k01.*A*k08*k11*k14*k16+k03*k06.*P*k01.*A*k11*k14*k16*k09+k03*k01.*A*k11*k14*k16*k09*k07.*B+k06.*P*k08*k10*k12.*Q*k02*k14*k16+k03*k06.*P*k01.*A*k08*k14*k16*k10+k03*k06.*P*k08*k10*k12.*Q*k14*k16;
E4 = k05*k08*k03*k01.*A*k11*k14*k16+k05*k03*k01.*A*k11*k14*k16*k09+k08*k10*k12.*Q*k02*k14*k16*k04+k05*k08*k10*k12.*Q*k02*k14*k16+k05*k08*k03*k10*k01.*A*k14*k16+k05*k08*k03*k10*k12.*Q*k14*k16;
E5 = k07.*B*k05*k03*k01.*A*k11*k14*k16+k10*k12.*Q*k02*k14*k16*k04*k06.*P+k07.*B*k10*k12.*Q*k02*k14*k16*k04+k07.*B*k10*k05*k12.*Q*k02*k14*k16+k07.*B*k10*k05*k03*k01.*A*k14*k16+k07.*B*k10*k05*k12.*Q*k03*k14*k16;
E6 = k12.*Q*k02*k14*k16*k04*k06.*P*k08+k12.*Q*k09*k02*k14*k16*k04*k06.*P+k12.*Q*k09*k02*k14*k16*k07.*B*k04+k12.*Q*k09*k02*k14*k16*k07.*B*k05+k09*k07.*B*k05*k03*k01.*A*k14*k16+k12.*Q*k09*k14*k16*k07.*B*k05*k03;
E7 = k13.*I1*k02*k11*k16*k04*k06.*P*k08+k13.*I1*k02*k11*k16*k04*k09*k06.*P+k13.*I1*k02*k11*k16*k04*k09*k07.*B+k13.*I1*k02*k11*k16*k09*k07.*B*k05+k13.*I1*k02*k16*k04*k06.*P*k08*k10+k13.*I1*k11*k16*k09*k07.*B*k05*k03;
E8 = k15.*I2*k02*k11*k14*k04*k06.*P*k08+k15.*I2*k02*k11*k14*k04*k09*k06.*P+k15.*I2*k02*k11*k14*k04*k09*k07.*B+k15.*I2*k02*k11*k14*k09*k07.*B*k05+k15.*I2*k02*k14*k04*k06.*P*k08*k10+k15.*I2*k11*k14*k09*k07.*B*k05*k03;
% Denominator terms
Denominator = E1+E2+E3+E4+E5+E6+E7+E8;
% Enzyme abundances terms
E1 = E1./Denominator;
E2 = E2./Denominator;
E3 = E3./Denominator;
E4 = E4./Denominator;
E5 = E5./Denominator;
E6 = E6./Denominator;
E7 = E7./Denominator;
E8 = E8./Denominator;
% Reaction rate 
v = +k05.*E3-k06.*P.*E4;