function y = b_car_rejection_model_simple_D010_Kinetics1_ode(y,Eref,metsRefConc,xconst,model,fixedExch,Sred,kinInactRxns,subunits,flag)
v = zeros(22,size(Eref,2));
E = zeros(22,size(Eref,2));
m_aacoa = y(1,:);
m_hmgcoa = y(2,:);
m_mev_R = y(3,:);
m_5pmev = y(4,:);
m_5dpmev = y(5,:);
m_ipdp = y(6,:);
m_dmpp = y(7,:);
m_grdp = y(8,:);
m_frdp = y(9,:);
m_ggdp = y(10,:);
m_psql = y(11,:);
m_pphy = y(12,:);
m_phy = y(13,:);
m_z_car = y(14,:);
m_lyc = y(15,:);
m_g_car = y(16,:);
m_b_car = y(17,:);
m_accoa = y(18,:);
m_coa = y(19,:);
m_nadph = y(20,:);
m_nadp = y(21,:);
m_atp = y(22,:);
m_adp = y(23,:);
m_pi = y(24,:);
m_co2 = y(25,:);
m_ppi = y(26,:);
m_sql = y(27,:);
m_fad = y(28,:);
m_fadh2 = y(29,:);
E(1,:) = Eref(1,:);
E(2,:) = Eref(2,:);
E(3,:) = Eref(3,:);
E(4,:) = Eref(4,:);
E(5,:) = Eref(5,:);
E(6,:) = Eref(6,:);
E(7,:) = Eref(7,:);
E(8,:) = Eref(8,:);
E(9,:) = Eref(9,:);
E(10,:) = Eref(10,:);
E(11,:) = Eref(11,:);
E(12,:) = Eref(12,:);
E(13,:) = Eref(13,:);
E(14,:) = Eref(14,:);
E(15,:) = Eref(15,:);
E(16,:) = Eref(16,:);
E(17,:) = Eref(17,:);
E(18,:) = Eref(18,:);
E(19,:) = Eref(19,:);
E(20,:) = Eref(20,:);
E(21,:) = Eref(21,:);
E(22,:) = Eref(22,:);
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

y = (1./(metsRefConc.*10^3)) .* (Sred*(E.*v));