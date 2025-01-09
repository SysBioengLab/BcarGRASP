function [f,grad] = b_car_rejection_model_simple_D025_Kinetics1(x,xconst,model,fixedExch,Sred,kinInactRxns,subunits,flag)
% Pre-allocation of memory
h = 1e-8;
% Defining metabolite and enzyme species
if flag==1
x = x(:);
xconst = xconst(:);
v = zeros(22,52);
E = zeros(22,52);
x = [x,x(:,ones(1,51)) + diag(h*1i*ones(51,1))];
xconst = [xconst,xconst(:,ones(1,51))];
else
v = zeros(22,size(x,2));
E = zeros(22,size(x,2));
end
% Defining metabolite and enzyme species
m_aacoa = x(1,:);
m_hmgcoa = x(2,:);
m_mev_R = x(3,:);
m_5pmev = x(4,:);
m_5dpmev = x(5,:);
m_ipdp = x(6,:);
m_dmpp = x(7,:);
m_grdp = x(8,:);
m_frdp = x(9,:);
m_ggdp = x(10,:);
m_psql = x(11,:);
m_pphy = x(12,:);
m_phy = x(13,:);
m_z_car = x(14,:);
m_lyc = x(15,:);
m_g_car = x(16,:);
m_b_car = x(17,:);
m_accoa = x(18,:);
m_coa = x(19,:);
m_nadph = x(20,:);
m_nadp = x(21,:);
m_atp = x(22,:);
m_adp = x(23,:);
m_pi = x(24,:);
m_co2 = x(25,:);
m_ppi = x(26,:);
m_sql = x(27,:);
m_fad = x(28,:);
m_fadh2 = x(29,:);
E(1,:) = x(30,:);
E(2,:) = x(31,:);
E(3,:) = x(32,:);
E(4,:) = x(33,:);
E(5,:) = x(34,:);
E(6,:) = x(35,:);
E(7,:) = x(36,:);
E(8,:) = x(37,:);
E(9,:) = x(38,:);
E(10,:) = x(39,:);
E(11,:) = x(40,:);
E(12,:) = x(41,:);
E(13,:) = x(42,:);
E(14,:) = x(43,:);
E(15,:) = x(44,:);
E(16,:) = x(45,:);
E(17,:) = x(46,:);
E(18,:) = x(47,:);
E(19,:) = x(48,:);
E(20,:) = x(49,:);
E(21,:) = x(50,:);
E(22,:) = x(51,:);
% Reaction rates
v(1,:) = r_ERG101([m_accoa;m_accoa;m_coa;m_aacoa],model.rxnParams(1).kineticParams);
v(2,:) = r_ERG131([m_accoa;m_aacoa;m_coa;m_hmgcoa],model.rxnParams(2).kineticParams);
v(3,:) = r_HMG11([m_hmgcoa;m_nadph;m_nadph;m_nadp;m_coa;m_nadp;m_mev_R],model.rxnParams(3).kineticParams);
v(4,:) = r_HMG21([m_hmgcoa;m_nadph;m_nadph;m_nadp;m_coa;m_nadp;m_mev_R],model.rxnParams(4).kineticParams);
v(5,:) = r_ERG121([m_mev_R;m_atp;m_5pmev;m_adp],model.rxnParams(5).kineticParams);
v(6,:) = r_ERG81([m_5pmev;m_atp;m_5dpmev;m_adp],model.rxnParams(6).kineticParams);
v(7,:) = r_MVD11([m_5dpmev;m_atp;m_adp;m_pi;m_co2;m_ipdp],model.rxnParams(7).kineticParams);
v(8,:) = r_IDI11([m_ipdp;m_dmpp],model.rxnParams(8).kineticParams);
v(9,:) = r_ERG20a1([m_dmpp;m_ipdp;m_ppi;m_grdp],model.rxnParams(9).kineticParams);
v(10,:) = r_ERG20b1([m_grdp;m_ipdp;m_ppi;m_frdp],model.rxnParams(10).kineticParams);
v(11,:) = r_BTS11([m_frdp;m_ipdp;m_ppi;m_ggdp],model.rxnParams(11).kineticParams);
v(12,:) = r_CrtE1([m_frdp;m_ipdp;m_ppi;m_ggdp],model.rxnParams(12).kineticParams);
v(13,:) = r_ERG9a1([m_frdp;m_frdp;m_ppi;m_psql],model.rxnParams(13).kineticParams);
v(14,:) = r_ERG9b1([m_psql;m_nadph;m_ppi;m_nadp;m_sql],model.rxnParams(14).kineticParams);
v(15,:) = r_CrtIa1([m_fad;m_fad;m_phy;m_fadh2;m_fadh2;m_z_car],model.rxnParams(15).kineticParams);
v(16,:) = r_CrtIb1([m_fad;m_fad;m_z_car;m_fadh2;m_fadh2;m_lyc],model.rxnParams(16).kineticParams);
v(17,:) = r_CrtYa1([m_lyc;m_g_car],model.rxnParams(17).kineticParams);
v(18,:) = r_CrtYb1([m_g_car;m_b_car],model.rxnParams(18).kineticParams);
v(19,:) = r_CrtBb1([m_pphy;m_ppi;m_phy],model.rxnParams(19).kineticParams);
v(20,:) = r_CrtBa1([m_ggdp;m_ggdp;m_ppi;m_pphy],model.rxnParams(20).kineticParams);
v(21,:) = r_EX_lyc1([m_lyc],model.rxnParams(21).kineticParams);
v(22,:) = r_EX_b_car1([m_b_car],model.rxnParams(22).kineticParams);
if flag==1
% Final rates
y = sum((Sred*(E.*v)).^2);
f = real(y(1));
if (nargout>1) % gradient is required
grad = imag(y(2:end))/h;
end
else
f = E.*v;
grad = [];
end