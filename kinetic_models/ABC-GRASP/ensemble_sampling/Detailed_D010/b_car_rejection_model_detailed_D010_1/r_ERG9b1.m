function [v,E1,E2] = r_ERG9b1(X,K)
% Metabolites definition 
A = X(1,:);
B = X(2,:);
C = X(3,:);
D = X(4,:);
P1 = X(5,:);
P2 = X(6,:);
Q = X(7,:);
R = X(8,:);
S = X(9,:);
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
k21 = K(21);
k22 = K(22);
k23 = K(23);
k24 = K(24);
%  Numerator terms
E1 = k02*k09*k12*k23*k04*k14*k06*k16.*P2*k18*k20+k02*k09*k12*k23*k04*k14*k21*k06*k16.*P2*k18+k02*k09*k12*k23*k04*k14*k21*k06*k16.*P2*k19+k02*k09*k12*k23*k04*k14*k21*k06*k19*k17.*D+k02*k09*k12*k23*k04*k21*k06*k19*k17.*D*k15+k02*k09*k12*k23*k04*k07*k14*k16.*P2*k18*k20+k02*k09*k12*k23*k04*k07*k14*k21*k16.*P2*k18+k02*k09*k12*k23*k04*k07*k14*k21*k16.*P2*k19+k02*k09*k12*k23*k04*k07*k14*k21*k19*k17.*D+k02*k09*k12*k23*k04*k07*k21*k19*k17.*D*k15+k02*k09*k12*k23*k07*k14*k05*k16.*P2*k18*k20+k02*k09*k12*k23*k07*k14*k21*k05*k16.*P2*k18+k02*k09*k12*k23*k07*k14*k21*k05*k16.*P2*k19+k02*k09*k12*k23*k07*k14*k21*k05*k19*k17.*D+k02*k09*k12*k23*k07*k21*k05*k19*k17.*D*k15+k02*k09*k12*k04*k14*k06*k16.*P2*k18*k20*k22.*R+k02*k09*k12*k04*k07*k14*k16.*P2*k18*k20*k22.*R+k02*k09*k12*k07*k14*k05*k16.*P2*k18*k20*k22.*R+k02*k09*k23*k04*k21*k06*k19*k17.*D*k15*k13+k02*k09*k23*k04*k07*k21*k19*k17.*D*k15*k13+k02*k09*k23*k07*k21*k05*k19*k17.*D*k15*k13+k02*k12*k23*k04*k14*k06*k16.*P2*k08.*P1*k18*k20+k02*k12*k23*k04*k14*k21*k06*k16.*P2*k08.*P1*k18+k02*k12*k23*k04*k14*k21*k06*k16.*P2*k19*k08.*P1+k02*k12*k23*k04*k14*k21*k06*k19*k08.*P1*k17.*D+k02*k12*k23*k04*k21*k06*k19*k08.*P1*k17.*D*k15+k02*k12*k04*k14*k06*k16.*P2*k08.*P1*k18*k20*k22.*R+k02*k23*k04*k21*k06*k19*k08.*P1*k17.*D*k15*k13+k09*k12*k23*k07*k14*k05*k16.*P2*k03.*B*k18*k20+k09*k12*k23*k07*k14*k21*k05*k16.*P2*k03.*B*k18+k09*k12*k23*k07*k14*k21*k05*k16.*P2*k19*k03.*B+k09*k12*k23*k07*k14*k21*k05*k19*k03.*B*k17.*D+k09*k12*k23*k07*k21*k05*k19*k03.*B*k17.*D*k15+k09*k12*k07*k14*k05*k16.*P2*k03.*B*k18*k20*k22.*R+k09*k23*k07*k21*k05*k19*k03.*B*k17.*D*k15*k13;
E2 = k01.*A*k04*k09*k12*k23*k06*k14*k16.*P2*k18*k20+k01.*A*k04*k09*k12*k23*k06*k14*k21*k16.*P2*k18+k01.*A*k04*k09*k12*k23*k06*k14*k21*k16.*P2*k19+k01.*A*k04*k09*k12*k23*k06*k14*k21*k19*k17.*D+k01.*A*k04*k09*k12*k23*k06*k21*k19*k17.*D*k15+k01.*A*k04*k09*k12*k23*k07*k14*k16.*P2*k18*k20+k01.*A*k04*k09*k12*k23*k07*k14*k21*k16.*P2*k18+k01.*A*k04*k09*k12*k23*k07*k14*k21*k16.*P2*k19+k01.*A*k04*k09*k12*k23*k07*k14*k21*k19*k17.*D+k01.*A*k04*k09*k12*k23*k07*k21*k19*k17.*D*k15+k01.*A*k09*k12*k23*k07*k14*k05*k16.*P2*k18*k20+k01.*A*k09*k12*k23*k07*k14*k21*k05*k16.*P2*k18+k01.*A*k09*k12*k23*k07*k14*k21*k05*k16.*P2*k19+k01.*A*k09*k12*k23*k07*k14*k21*k05*k19*k17.*D+k01.*A*k09*k12*k23*k07*k21*k05*k19*k17.*D*k15+k01.*A*k04*k09*k12*k06*k14*k16.*P2*k18*k20*k22.*R+k01.*A*k04*k09*k12*k07*k14*k16.*P2*k18*k20*k22.*R+k01.*A*k09*k12*k07*k14*k05*k16.*P2*k18*k20*k22.*R+k01.*A*k04*k09*k23*k06*k21*k19*k17.*D*k15*k13+k01.*A*k04*k09*k23*k07*k21*k19*k17.*D*k15*k13+k01.*A*k09*k23*k07*k21*k05*k19*k17.*D*k15*k13+k01.*A*k04*k12*k23*k06*k14*k08.*P1*k16.*P2*k18*k20+k01.*A*k04*k12*k23*k06*k14*k21*k08.*P1*k16.*P2*k18+k01.*A*k04*k12*k23*k06*k14*k21*k08.*P1*k16.*P2*k19+k01.*A*k04*k12*k23*k06*k14*k21*k08.*P1*k19*k17.*D+k01.*A*k04*k12*k23*k06*k21*k08.*P1*k19*k17.*D*k15+k01.*A*k04*k12*k06*k14*k08.*P1*k16.*P2*k18*k20*k22.*R+k01.*A*k04*k23*k06*k21*k08.*P1*k19*k17.*D*k15*k13+k04*k06*k08.*P1*k10.*Q*k12*k23*k14*k16.*P2*k18*k20+k04*k06*k08.*P1*k10.*Q*k12*k23*k14*k21*k16.*P2*k18+k04*k06*k08.*P1*k10.*Q*k12*k23*k14*k21*k16.*P2*k19+k04*k06*k08.*P1*k10.*Q*k12*k23*k14*k21*k19*k17.*D+k04*k06*k08.*P1*k10.*Q*k12*k23*k21*k19*k17.*D*k15+k04*k06*k08.*P1*k10.*Q*k12*k14*k16.*P2*k18*k20*k22.*R+k04*k06*k08.*P1*k10.*Q*k23*k21*k19*k17.*D*k15*k13;
E3 = k03.*B*k06*k01.*A*k09*k12*k23*k14*k16.*P2*k18*k20+k03.*B*k06*k01.*A*k09*k12*k23*k14*k21*k16.*P2*k18+k03.*B*k06*k01.*A*k09*k12*k23*k14*k21*k16.*P2*k19+k03.*B*k06*k01.*A*k09*k12*k23*k14*k21*k19*k17.*D+k03.*B*k06*k01.*A*k09*k12*k23*k21*k19*k17.*D*k15+k03.*B*k01.*A*k09*k12*k23*k07*k14*k16.*P2*k18*k20+k03.*B*k01.*A*k09*k12*k23*k07*k14*k21*k16.*P2*k18+k03.*B*k01.*A*k09*k12*k23*k07*k14*k21*k16.*P2*k19+k03.*B*k01.*A*k09*k12*k23*k07*k14*k21*k19*k17.*D+k03.*B*k01.*A*k09*k12*k23*k07*k21*k19*k17.*D*k15+k06*k08.*P1*k10.*Q*k02*k12*k23*k14*k16.*P2*k18*k20+k06*k08.*P1*k10.*Q*k02*k12*k23*k14*k21*k16.*P2*k18+k06*k08.*P1*k10.*Q*k02*k12*k23*k14*k21*k16.*P2*k19+k06*k08.*P1*k10.*Q*k02*k12*k23*k14*k21*k19*k17.*D+k06*k08.*P1*k10.*Q*k02*k12*k23*k21*k19*k17.*D*k15+k03.*B*k06*k01.*A*k09*k12*k14*k16.*P2*k18*k20*k22.*R+k03.*B*k01.*A*k09*k12*k07*k14*k16.*P2*k18*k20*k22.*R+k06*k08.*P1*k10.*Q*k02*k12*k14*k16.*P2*k18*k20*k22.*R+k03.*B*k06*k01.*A*k09*k23*k21*k19*k17.*D*k15*k13+k03.*B*k01.*A*k09*k23*k07*k21*k19*k17.*D*k15*k13+k06*k08.*P1*k10.*Q*k02*k23*k21*k19*k17.*D*k15*k13+k03.*B*k06*k01.*A*k08.*P1*k12*k23*k14*k16.*P2*k18*k20+k03.*B*k06*k01.*A*k08.*P1*k12*k23*k14*k21*k16.*P2*k18+k03.*B*k06*k01.*A*k08.*P1*k12*k23*k14*k21*k16.*P2*k19+k03.*B*k06*k01.*A*k08.*P1*k12*k23*k14*k21*k19*k17.*D+k03.*B*k06*k01.*A*k08.*P1*k12*k23*k21*k19*k17.*D*k15+k03.*B*k06*k01.*A*k08.*P1*k12*k14*k16.*P2*k18*k20*k22.*R+k03.*B*k06*k01.*A*k08.*P1*k23*k21*k19*k17.*D*k15*k13+k03.*B*k06*k08.*P1*k10.*Q*k12*k23*k14*k16.*P2*k18*k20+k03.*B*k06*k08.*P1*k10.*Q*k12*k23*k14*k21*k16.*P2*k18+k03.*B*k06*k08.*P1*k10.*Q*k12*k23*k14*k21*k16.*P2*k19+k03.*B*k06*k08.*P1*k10.*Q*k12*k23*k14*k21*k19*k17.*D+k03.*B*k06*k08.*P1*k10.*Q*k12*k23*k21*k19*k17.*D*k15+k03.*B*k06*k08.*P1*k10.*Q*k12*k14*k16.*P2*k18*k20*k22.*R+k03.*B*k06*k08.*P1*k10.*Q*k23*k21*k19*k17.*D*k15*k13;
E4 = k05*k03.*B*k01.*A*k09*k12*k23*k14*k16.*P2*k18*k20+k05*k03.*B*k01.*A*k09*k12*k23*k14*k21*k16.*P2*k18+k05*k03.*B*k01.*A*k09*k12*k23*k14*k21*k16.*P2*k19+k05*k03.*B*k01.*A*k09*k12*k23*k14*k21*k19*k17.*D+k05*k03.*B*k01.*A*k09*k12*k23*k21*k19*k17.*D*k15+k08.*P1*k10.*Q*k02*k12*k23*k04*k14*k16.*P2*k18*k20+k08.*P1*k10.*Q*k02*k12*k23*k04*k14*k21*k16.*P2*k18+k08.*P1*k10.*Q*k02*k12*k23*k04*k14*k21*k16.*P2*k19+k08.*P1*k10.*Q*k02*k12*k23*k04*k14*k21*k19*k17.*D+k08.*P1*k10.*Q*k02*k12*k23*k04*k21*k19*k17.*D*k15+k05*k08.*P1*k10.*Q*k02*k12*k23*k14*k16.*P2*k18*k20+k05*k08.*P1*k10.*Q*k02*k12*k23*k14*k21*k16.*P2*k18+k05*k08.*P1*k10.*Q*k02*k12*k23*k14*k21*k16.*P2*k19+k05*k08.*P1*k10.*Q*k02*k12*k23*k14*k21*k19*k17.*D+k05*k08.*P1*k10.*Q*k02*k12*k23*k21*k19*k17.*D*k15+k05*k03.*B*k01.*A*k09*k12*k14*k16.*P2*k18*k20*k22.*R+k08.*P1*k10.*Q*k02*k12*k04*k14*k16.*P2*k18*k20*k22.*R+k05*k08.*P1*k10.*Q*k02*k12*k14*k16.*P2*k18*k20*k22.*R+k05*k03.*B*k01.*A*k09*k23*k21*k19*k17.*D*k15*k13+k08.*P1*k10.*Q*k02*k23*k04*k21*k19*k17.*D*k15*k13+k05*k08.*P1*k10.*Q*k02*k23*k21*k19*k17.*D*k15*k13+k05*k08.*P1*k03.*B*k01.*A*k12*k23*k14*k16.*P2*k18*k20+k05*k08.*P1*k03.*B*k01.*A*k12*k23*k14*k21*k16.*P2*k18+k05*k08.*P1*k03.*B*k01.*A*k12*k23*k14*k21*k16.*P2*k19+k05*k08.*P1*k03.*B*k01.*A*k12*k23*k14*k21*k19*k17.*D+k05*k08.*P1*k03.*B*k01.*A*k12*k23*k21*k19*k17.*D*k15+k05*k08.*P1*k03.*B*k01.*A*k12*k14*k16.*P2*k18*k20*k22.*R+k05*k08.*P1*k03.*B*k01.*A*k23*k21*k19*k17.*D*k15*k13+k05*k08.*P1*k03.*B*k10.*Q*k12*k23*k14*k16.*P2*k18*k20+k05*k08.*P1*k03.*B*k10.*Q*k12*k23*k14*k21*k16.*P2*k18+k05*k08.*P1*k03.*B*k10.*Q*k12*k23*k14*k21*k16.*P2*k19+k05*k08.*P1*k03.*B*k10.*Q*k12*k23*k14*k21*k19*k17.*D+k05*k08.*P1*k03.*B*k10.*Q*k12*k23*k21*k19*k17.*D*k15+k05*k08.*P1*k03.*B*k10.*Q*k12*k14*k16.*P2*k18*k20*k22.*R+k05*k08.*P1*k03.*B*k10.*Q*k23*k21*k19*k17.*D*k15*k13;
E5 = k10.*Q*k02*k12*k23*k04*k14*k06*k16.*P2*k18*k20+k10.*Q*k02*k12*k23*k04*k14*k21*k06*k16.*P2*k18+k10.*Q*k02*k12*k23*k04*k14*k21*k06*k16.*P2*k19+k10.*Q*k02*k12*k23*k04*k14*k21*k06*k19*k17.*D+k10.*Q*k02*k12*k23*k04*k21*k06*k19*k17.*D*k15+k10.*Q*k07*k02*k12*k23*k04*k14*k16.*P2*k18*k20+k10.*Q*k07*k02*k12*k23*k04*k14*k21*k16.*P2*k18+k10.*Q*k07*k02*k12*k23*k04*k14*k21*k16.*P2*k19+k10.*Q*k07*k02*k12*k23*k04*k14*k21*k19*k17.*D+k10.*Q*k07*k02*k12*k23*k04*k21*k19*k17.*D*k15+k10.*Q*k07*k02*k12*k23*k05*k14*k16.*P2*k18*k20+k10.*Q*k07*k02*k12*k23*k05*k14*k21*k16.*P2*k18+k10.*Q*k07*k02*k12*k23*k05*k14*k21*k16.*P2*k19+k10.*Q*k07*k02*k12*k23*k05*k14*k21*k19*k17.*D+k10.*Q*k07*k02*k12*k23*k05*k21*k19*k17.*D*k15+k10.*Q*k02*k12*k04*k14*k06*k16.*P2*k18*k20*k22.*R+k10.*Q*k07*k02*k12*k04*k14*k16.*P2*k18*k20*k22.*R+k10.*Q*k07*k02*k12*k05*k14*k16.*P2*k18*k20*k22.*R+k10.*Q*k02*k23*k04*k21*k06*k19*k17.*D*k15*k13+k10.*Q*k07*k02*k23*k04*k21*k19*k17.*D*k15*k13+k10.*Q*k07*k02*k23*k05*k21*k19*k17.*D*k15*k13+k07*k05*k03.*B*k01.*A*k12*k23*k14*k16.*P2*k18*k20+k07*k05*k03.*B*k01.*A*k12*k23*k14*k21*k16.*P2*k18+k07*k05*k03.*B*k01.*A*k12*k23*k14*k21*k16.*P2*k19+k07*k05*k03.*B*k01.*A*k12*k23*k14*k21*k19*k17.*D+k07*k05*k03.*B*k01.*A*k12*k23*k21*k19*k17.*D*k15+k07*k05*k03.*B*k01.*A*k12*k14*k16.*P2*k18*k20*k22.*R+k07*k05*k03.*B*k01.*A*k23*k21*k19*k17.*D*k15*k13+k10.*Q*k07*k12*k23*k05*k14*k03.*B*k16.*P2*k18*k20+k10.*Q*k07*k12*k23*k05*k14*k21*k03.*B*k16.*P2*k18+k10.*Q*k07*k12*k23*k05*k14*k21*k03.*B*k16.*P2*k19+k10.*Q*k07*k12*k23*k05*k14*k21*k03.*B*k19*k17.*D+k10.*Q*k07*k12*k23*k05*k21*k03.*B*k19*k17.*D*k15+k10.*Q*k07*k12*k05*k14*k03.*B*k16.*P2*k18*k20*k22.*R+k10.*Q*k07*k23*k05*k21*k03.*B*k19*k17.*D*k15*k13;
E6 = k11.*C*k14*k02*k09*k23*k16.*P2*k04*k18*k06*k20+k11.*C*k14*k02*k09*k23*k16.*P2*k04*k21*k18*k06+k11.*C*k14*k02*k09*k23*k16.*P2*k04*k21*k06*k19+k11.*C*k14*k02*k09*k23*k04*k21*k06*k19*k17.*D+k11.*C*k02*k09*k23*k04*k21*k06*k19*k17.*D*k15+k11.*C*k14*k02*k09*k23*k16.*P2*k04*k07*k18*k20+k11.*C*k14*k02*k09*k23*k16.*P2*k04*k07*k21*k18+k11.*C*k14*k02*k09*k23*k16.*P2*k04*k07*k21*k19+k11.*C*k14*k02*k09*k23*k04*k07*k21*k19*k17.*D+k11.*C*k02*k09*k23*k04*k07*k21*k19*k17.*D*k15+k11.*C*k14*k02*k09*k23*k16.*P2*k07*k18*k05*k20+k11.*C*k14*k02*k09*k23*k16.*P2*k07*k21*k18*k05+k11.*C*k14*k02*k09*k23*k16.*P2*k07*k21*k05*k19+k11.*C*k14*k02*k09*k23*k07*k21*k05*k19*k17.*D+k11.*C*k02*k09*k23*k07*k21*k05*k19*k17.*D*k15+k11.*C*k14*k02*k09*k16.*P2*k04*k18*k06*k20*k22.*R+k11.*C*k14*k02*k09*k16.*P2*k04*k07*k18*k20*k22.*R+k11.*C*k14*k02*k09*k16.*P2*k07*k18*k05*k20*k22.*R+k14*k16.*P2*k18*k20*k22.*R*k24.*S*k02*k09*k04*k06+k14*k16.*P2*k18*k20*k22.*R*k24.*S*k02*k09*k04*k07+k14*k16.*P2*k18*k20*k22.*R*k24.*S*k02*k09*k07*k05+k11.*C*k14*k02*k23*k16.*P2*k04*k18*k06*k20*k08.*P1+k11.*C*k14*k02*k23*k16.*P2*k04*k21*k18*k06*k08.*P1+k11.*C*k14*k02*k23*k16.*P2*k04*k21*k06*k19*k08.*P1+k11.*C*k14*k02*k23*k04*k21*k06*k19*k08.*P1*k17.*D+k11.*C*k02*k23*k04*k21*k06*k19*k08.*P1*k17.*D*k15+k11.*C*k14*k02*k16.*P2*k04*k18*k06*k20*k08.*P1*k22.*R+k14*k16.*P2*k18*k20*k22.*R*k24.*S*k02*k04*k06*k08.*P1+k11.*C*k14*k09*k23*k16.*P2*k07*k18*k05*k20*k03.*B+k11.*C*k14*k09*k23*k16.*P2*k07*k21*k18*k05*k03.*B+k11.*C*k14*k09*k23*k16.*P2*k07*k21*k05*k19*k03.*B+k11.*C*k14*k09*k23*k07*k21*k05*k19*k03.*B*k17.*D+k11.*C*k09*k23*k07*k21*k05*k19*k03.*B*k17.*D*k15+k11.*C*k14*k09*k16.*P2*k07*k18*k05*k20*k03.*B*k22.*R+k14*k16.*P2*k18*k20*k22.*R*k24.*S*k09*k07*k05*k03.*B;
E7 = k13*k16.*P2*k11.*C*k18*k02*k09*k23*k20*k04*k06+k13*k16.*P2*k11.*C*k18*k02*k09*k23*k04*k21*k06+k13*k16.*P2*k11.*C*k02*k09*k23*k04*k21*k06*k19+k13*k11.*C*k02*k09*k23*k04*k21*k06*k19*k17.*D+k16.*P2*k18*k20*k22.*R*k24.*S*k02*k09*k12*k04*k06+k13*k16.*P2*k11.*C*k18*k02*k09*k23*k20*k04*k07+k13*k16.*P2*k11.*C*k18*k02*k09*k23*k04*k07*k21+k13*k16.*P2*k11.*C*k02*k09*k23*k04*k07*k21*k19+k13*k11.*C*k02*k09*k23*k04*k07*k21*k19*k17.*D+k16.*P2*k18*k20*k22.*R*k24.*S*k02*k09*k12*k04*k07+k13*k16.*P2*k11.*C*k18*k02*k09*k23*k20*k07*k05+k13*k16.*P2*k11.*C*k18*k02*k09*k23*k07*k21*k05+k13*k16.*P2*k11.*C*k02*k09*k23*k07*k21*k05*k19+k13*k11.*C*k02*k09*k23*k07*k21*k05*k19*k17.*D+k16.*P2*k18*k20*k22.*R*k24.*S*k02*k09*k12*k07*k05+k13*k16.*P2*k11.*C*k18*k02*k09*k20*k04*k22.*R*k06+k13*k16.*P2*k11.*C*k18*k02*k09*k20*k04*k07*k22.*R+k13*k16.*P2*k11.*C*k18*k02*k09*k20*k07*k22.*R*k05+k13*k16.*P2*k18*k20*k22.*R*k24.*S*k02*k09*k04*k06+k13*k16.*P2*k18*k20*k22.*R*k24.*S*k02*k09*k04*k07+k13*k16.*P2*k18*k20*k22.*R*k24.*S*k02*k09*k07*k05+k13*k16.*P2*k11.*C*k18*k02*k23*k20*k04*k06*k08.*P1+k13*k16.*P2*k11.*C*k18*k02*k23*k04*k21*k06*k08.*P1+k13*k16.*P2*k11.*C*k02*k23*k04*k21*k06*k19*k08.*P1+k13*k11.*C*k02*k23*k04*k21*k06*k19*k08.*P1*k17.*D+k16.*P2*k18*k20*k22.*R*k24.*S*k02*k12*k04*k06*k08.*P1+k13*k16.*P2*k11.*C*k18*k02*k20*k04*k22.*R*k06*k08.*P1+k13*k16.*P2*k18*k20*k22.*R*k24.*S*k02*k04*k06*k08.*P1+k13*k16.*P2*k11.*C*k18*k09*k23*k20*k07*k05*k03.*B+k13*k16.*P2*k11.*C*k18*k09*k23*k07*k21*k05*k03.*B+k13*k16.*P2*k11.*C*k09*k23*k07*k21*k05*k19*k03.*B+k13*k11.*C*k09*k23*k07*k21*k05*k19*k03.*B*k17.*D+k16.*P2*k18*k20*k22.*R*k24.*S*k09*k12*k07*k05*k03.*B+k13*k16.*P2*k11.*C*k18*k09*k20*k07*k22.*R*k05*k03.*B+k13*k16.*P2*k18*k20*k22.*R*k24.*S*k09*k07*k05*k03.*B;
E8 = k15*k18*k13*k20*k11.*C*k02*k09*k23*k04*k06+k15*k18*k13*k11.*C*k02*k09*k23*k04*k21*k06+k15*k13*k11.*C*k02*k09*k23*k04*k21*k06*k19+k18*k20*k22.*R*k24.*S*k02*k09*k12*k04*k14*k06+k15*k18*k20*k22.*R*k24.*S*k02*k09*k12*k04*k06+k15*k18*k13*k20*k11.*C*k02*k09*k23*k04*k07+k15*k18*k13*k11.*C*k02*k09*k23*k04*k07*k21+k15*k13*k11.*C*k02*k09*k23*k04*k07*k21*k19+k18*k20*k22.*R*k24.*S*k02*k09*k12*k04*k07*k14+k15*k18*k20*k22.*R*k24.*S*k02*k09*k12*k04*k07+k15*k18*k13*k20*k11.*C*k02*k09*k23*k07*k05+k15*k18*k13*k11.*C*k02*k09*k23*k07*k21*k05+k15*k13*k11.*C*k02*k09*k23*k07*k21*k05*k19+k18*k20*k22.*R*k24.*S*k02*k09*k12*k07*k14*k05+k15*k18*k20*k22.*R*k24.*S*k02*k09*k12*k07*k05+k15*k18*k13*k20*k11.*C*k22.*R*k02*k09*k04*k06+k15*k18*k13*k20*k11.*C*k22.*R*k02*k09*k04*k07+k15*k18*k13*k20*k11.*C*k22.*R*k02*k09*k07*k05+k15*k18*k13*k20*k22.*R*k24.*S*k02*k09*k04*k06+k15*k18*k13*k20*k22.*R*k24.*S*k02*k09*k04*k07+k15*k18*k13*k20*k22.*R*k24.*S*k02*k09*k07*k05+k15*k18*k13*k20*k11.*C*k02*k23*k04*k06*k08.*P1+k15*k18*k13*k11.*C*k02*k23*k04*k21*k06*k08.*P1+k15*k13*k11.*C*k02*k23*k04*k21*k06*k19*k08.*P1+k18*k20*k22.*R*k24.*S*k02*k12*k04*k14*k06*k08.*P1+k15*k18*k20*k22.*R*k24.*S*k02*k12*k04*k06*k08.*P1+k15*k18*k13*k20*k11.*C*k22.*R*k02*k04*k06*k08.*P1+k15*k18*k13*k20*k22.*R*k24.*S*k02*k04*k06*k08.*P1+k15*k18*k13*k20*k11.*C*k09*k23*k07*k05*k03.*B+k15*k18*k13*k11.*C*k09*k23*k07*k21*k05*k03.*B+k15*k13*k11.*C*k09*k23*k07*k21*k05*k19*k03.*B+k18*k20*k22.*R*k24.*S*k09*k12*k07*k14*k05*k03.*B+k15*k18*k20*k22.*R*k24.*S*k09*k12*k07*k05*k03.*B+k15*k18*k13*k20*k11.*C*k22.*R*k09*k07*k05*k03.*B+k15*k18*k13*k20*k22.*R*k24.*S*k09*k07*k05*k03.*B;
E9 = k17.*D*k20*k15*k13*k11.*C*k02*k09*k23*k04*k06+k17.*D*k15*k13*k11.*C*k02*k09*k23*k04*k21*k06+k20*k22.*R*k24.*S*k02*k09*k12*k04*k14*k06*k16.*P2+k17.*D*k20*k22.*R*k24.*S*k02*k09*k12*k04*k14*k06+k17.*D*k20*k15*k22.*R*k24.*S*k02*k09*k12*k04*k06+k17.*D*k20*k15*k13*k11.*C*k02*k09*k23*k04*k07+k17.*D*k15*k13*k11.*C*k02*k09*k23*k04*k07*k21+k20*k22.*R*k24.*S*k02*k09*k12*k04*k07*k14*k16.*P2+k17.*D*k20*k22.*R*k24.*S*k02*k09*k12*k04*k07*k14+k17.*D*k20*k15*k22.*R*k24.*S*k02*k09*k12*k04*k07+k17.*D*k20*k15*k13*k11.*C*k02*k09*k23*k07*k05+k17.*D*k15*k13*k11.*C*k02*k09*k23*k07*k21*k05+k20*k22.*R*k24.*S*k02*k09*k12*k07*k14*k05*k16.*P2+k17.*D*k20*k22.*R*k24.*S*k02*k09*k12*k07*k14*k05+k17.*D*k20*k15*k22.*R*k24.*S*k02*k09*k12*k07*k05+k17.*D*k20*k15*k22.*R*k13*k11.*C*k02*k09*k04*k06+k17.*D*k20*k15*k22.*R*k13*k11.*C*k02*k09*k04*k07+k17.*D*k20*k15*k22.*R*k13*k11.*C*k02*k09*k07*k05+k17.*D*k20*k15*k22.*R*k13*k24.*S*k02*k09*k04*k06+k17.*D*k20*k15*k22.*R*k13*k24.*S*k02*k09*k04*k07+k17.*D*k20*k15*k22.*R*k13*k24.*S*k02*k09*k07*k05+k17.*D*k20*k15*k13*k11.*C*k02*k23*k04*k06*k08.*P1+k17.*D*k15*k13*k11.*C*k02*k23*k04*k21*k06*k08.*P1+k20*k22.*R*k24.*S*k02*k12*k04*k14*k06*k16.*P2*k08.*P1+k17.*D*k20*k22.*R*k24.*S*k02*k12*k04*k14*k06*k08.*P1+k17.*D*k20*k15*k22.*R*k24.*S*k02*k12*k04*k06*k08.*P1+k17.*D*k20*k15*k22.*R*k13*k11.*C*k02*k04*k06*k08.*P1+k17.*D*k20*k15*k22.*R*k13*k24.*S*k02*k04*k06*k08.*P1+k17.*D*k20*k15*k13*k11.*C*k09*k23*k07*k05*k03.*B+k17.*D*k15*k13*k11.*C*k09*k23*k07*k21*k05*k03.*B+k20*k22.*R*k24.*S*k09*k12*k07*k14*k05*k16.*P2*k03.*B+k17.*D*k20*k22.*R*k24.*S*k09*k12*k07*k14*k05*k03.*B+k17.*D*k20*k15*k22.*R*k24.*S*k09*k12*k07*k05*k03.*B+k17.*D*k20*k15*k22.*R*k13*k11.*C*k09*k07*k05*k03.*B+k17.*D*k20*k15*k22.*R*k13*k24.*S*k09*k07*k05*k03.*B;
E10 = k19*k17.*D*k15*k13*k11.*C*k02*k09*k23*k04*k06+k22.*R*k24.*S*k02*k09*k12*k04*k14*k06*k16.*P2*k18+k19*k22.*R*k24.*S*k02*k09*k12*k04*k14*k06*k16.*P2+k19*k22.*R*k17.*D*k24.*S*k02*k09*k12*k04*k14*k06+k19*k22.*R*k17.*D*k24.*S*k15*k02*k09*k12*k04*k06+k19*k17.*D*k15*k13*k11.*C*k02*k09*k23*k04*k07+k22.*R*k24.*S*k02*k09*k12*k04*k07*k14*k16.*P2*k18+k19*k22.*R*k24.*S*k02*k09*k12*k04*k07*k14*k16.*P2+k19*k22.*R*k17.*D*k24.*S*k02*k09*k12*k04*k07*k14+k19*k22.*R*k17.*D*k24.*S*k15*k02*k09*k12*k04*k07+k19*k17.*D*k15*k13*k11.*C*k02*k09*k23*k07*k05+k22.*R*k24.*S*k02*k09*k12*k07*k14*k05*k16.*P2*k18+k19*k22.*R*k24.*S*k02*k09*k12*k07*k14*k05*k16.*P2+k19*k22.*R*k17.*D*k24.*S*k02*k09*k12*k07*k14*k05+k19*k22.*R*k17.*D*k24.*S*k15*k02*k09*k12*k07*k05+k19*k22.*R*k17.*D*k15*k13*k11.*C*k02*k09*k04*k06+k19*k22.*R*k17.*D*k15*k13*k11.*C*k02*k09*k04*k07+k19*k22.*R*k17.*D*k15*k13*k11.*C*k02*k09*k07*k05+k19*k22.*R*k17.*D*k24.*S*k15*k02*k09*k13*k04*k06+k19*k22.*R*k17.*D*k24.*S*k15*k02*k09*k13*k04*k07+k19*k22.*R*k17.*D*k24.*S*k15*k02*k09*k13*k07*k05+k19*k17.*D*k15*k13*k11.*C*k02*k23*k04*k06*k08.*P1+k22.*R*k24.*S*k02*k12*k04*k14*k06*k16.*P2*k08.*P1*k18+k19*k22.*R*k24.*S*k02*k12*k04*k14*k06*k16.*P2*k08.*P1+k19*k22.*R*k17.*D*k24.*S*k02*k12*k04*k14*k06*k08.*P1+k19*k22.*R*k17.*D*k24.*S*k15*k02*k12*k04*k06*k08.*P1+k19*k22.*R*k17.*D*k15*k13*k11.*C*k02*k04*k06*k08.*P1+k19*k22.*R*k17.*D*k24.*S*k15*k02*k13*k04*k06*k08.*P1+k19*k17.*D*k15*k13*k11.*C*k09*k23*k07*k05*k03.*B+k22.*R*k24.*S*k09*k12*k07*k14*k05*k16.*P2*k03.*B*k18+k19*k22.*R*k24.*S*k09*k12*k07*k14*k05*k16.*P2*k03.*B+k19*k22.*R*k17.*D*k24.*S*k09*k12*k07*k14*k05*k03.*B+k19*k22.*R*k17.*D*k24.*S*k15*k09*k12*k07*k05*k03.*B+k19*k22.*R*k17.*D*k15*k13*k11.*C*k09*k07*k05*k03.*B+k19*k22.*R*k17.*D*k24.*S*k15*k09*k13*k07*k05*k03.*B;
E11 = k24.*S*k02*k09*k12*k04*k14*k06*k16.*P2*k18*k20+k24.*S*k21*k02*k09*k12*k04*k14*k06*k16.*P2*k18+k24.*S*k21*k02*k09*k12*k19*k04*k14*k06*k16.*P2+k24.*S*k21*k02*k09*k12*k19*k04*k14*k17.*D*k06+k24.*S*k21*k02*k09*k12*k19*k04*k17.*D*k06*k15+k24.*S*k02*k09*k12*k04*k07*k14*k16.*P2*k18*k20+k24.*S*k21*k02*k09*k12*k04*k07*k14*k16.*P2*k18+k24.*S*k21*k02*k09*k12*k19*k04*k07*k14*k16.*P2+k24.*S*k21*k02*k09*k12*k19*k04*k07*k14*k17.*D+k24.*S*k21*k02*k09*k12*k19*k04*k07*k17.*D*k15+k24.*S*k02*k09*k12*k07*k14*k05*k16.*P2*k18*k20+k24.*S*k21*k02*k09*k12*k07*k14*k05*k16.*P2*k18+k24.*S*k21*k02*k09*k12*k19*k07*k14*k05*k16.*P2+k24.*S*k21*k02*k09*k12*k19*k07*k14*k17.*D*k05+k24.*S*k21*k02*k09*k12*k19*k07*k17.*D*k05*k15+k21*k19*k17.*D*k15*k13*k11.*C*k02*k09*k04*k06+k21*k19*k17.*D*k15*k13*k11.*C*k02*k09*k04*k07+k21*k19*k17.*D*k15*k13*k11.*C*k02*k09*k07*k05+k24.*S*k21*k02*k09*k19*k04*k17.*D*k06*k15*k13+k24.*S*k21*k02*k09*k19*k04*k07*k17.*D*k15*k13+k24.*S*k21*k02*k09*k19*k07*k17.*D*k05*k15*k13+k24.*S*k02*k12*k04*k14*k06*k16.*P2*k08.*P1*k18*k20+k24.*S*k21*k02*k12*k04*k14*k06*k16.*P2*k08.*P1*k18+k24.*S*k21*k02*k12*k19*k04*k14*k06*k16.*P2*k08.*P1+k24.*S*k21*k02*k12*k19*k04*k14*k17.*D*k06*k08.*P1+k24.*S*k21*k02*k12*k19*k04*k17.*D*k06*k15*k08.*P1+k21*k19*k17.*D*k15*k13*k11.*C*k02*k04*k06*k08.*P1+k24.*S*k21*k02*k19*k04*k17.*D*k06*k15*k08.*P1*k13+k24.*S*k09*k12*k07*k14*k05*k16.*P2*k03.*B*k18*k20+k24.*S*k21*k09*k12*k07*k14*k05*k16.*P2*k03.*B*k18+k24.*S*k21*k09*k12*k19*k07*k14*k05*k16.*P2*k03.*B+k24.*S*k21*k09*k12*k19*k07*k14*k17.*D*k05*k03.*B+k24.*S*k21*k09*k12*k19*k07*k17.*D*k05*k15*k03.*B+k21*k19*k17.*D*k15*k13*k11.*C*k09*k07*k05*k03.*B+k24.*S*k21*k09*k19*k07*k17.*D*k05*k15*k03.*B*k13;
% Denominator terms
Denominator = E1+E2+E3+E4+E5+E6+E7+E8+E9+E10+E11;
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
E10 = E10./Denominator;
E11 = E11./Denominator;
% Reaction rate 
v = +k15.*E7-k16.*P2.*E8;