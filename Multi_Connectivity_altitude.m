clc
clear
close all

% system parameters

fc_b = 2;          % carrier frequency of BS (GHz)
fc_u = 2;          % carrier frequency of UAV (GHz)
fc_h = 2;          % carrier frequency of HAPS (GHz)  
fc_s_S = 2;        % carrier frequency of satellite in S-band (GHz)
fc_s_Ka = 30;      % carrier frequency of satellite in Ka-band (GHz)

BW_RB = 0.2e6;       % bandwidth of each PRB (Hz)
Nc = 4;              % number of PRB allocated to each user (allocated bandwidth): {1,...,6}
Nf = 1;              % frequency diversity factor for BU links: {1,2}

BW_bu_tot = 10e6;    % total bandwidth of BU (Hz) ---> 50 PRB
BW_bu = Nc*BW_RB;    % bandwidth of each BU link (Hz): maximum bandwidth of URLLC user = coherence BW (1.2 MHz)

BW_uu_tot = 10e6;    % total bandwidth of UU (Hz)
BW_uu = Nc*BW_RB;    % bandwidth of each UU link (Hz)

BW_hu_tot = 10e6;    % total bandwidth of HU (Hz)
BW_hu = Nc*BW_RB;    % bandwidth of each HU link (Hz)

BW_su_tot = 10e6;    % total bandwidth of SU (Hz)
BW_su = Nc*BW_RB;    % bandwidth of each SU link (Hz)

BW_bh = 1e6;       % bandwidth of BH link (Hz)
BW_gs = 1e6;       % bandwidth of GS link (Hz)

R = 500e3;        % data rate (bps)

Pb = 10^(46/10);   % BS Tx power (mW)
Pu = 10^(23/10);   % UAV Tx power (mW)
Ph = 10^(46/10);   % HAPS Tx power (mW)
Pg = 10^(50/10);   % gateway Tx power (mW)
Ps = 10^(50/10);   % sat Tx power (mW)

thetaT = 102*pi/180;         % BS downtilt angle
thetaT_gat = 15*pi/180;      % HAP/Sat gateway uptilt angle
ge_max = 8;        % maximum element gain (dBi)
Ne = 8;            % number of elements of BS antenna
Gt_hmax = 32;      % HAPS max Tx antenna gain
Gr_h = 32;         % HAPS Rx antenna gain
Gt_u = 0;          % UAV Tx antenna gain
Gr_u = 0;          % UAV Rx antenna gain
Gt_smax = 38;      % Sat max Tx antenna gain
Gr_s = 38;         % Sat Rx antenna gain
Gt_g = 41;         % gateway Tx antenna gain

NF_b = 5;          % BS Rx noise figure (dB)
NF_u = 9;          % UAV Rx noise figure
NF_h = 5;          % HAPS Rx noise figure
NF_s = 5;          % Sat Rx noise figure

N0 = 10^(-174/10);  % single-side noise (mW/Hz)
packet_size = 256;  % bits (32 bytes)
Re = 6371e3;        % radius of Earth (m)
c = 3e8;            % speed of light (m/s)

lambda = 2e3;       % average arrival rate of packets (packet/s)
eff_bw = 4*lambda;  % effective bandwidth (packet/s)
Dq_max = 0.7e-3;    % qeueuing delay bound (s) (or 0.8ms if it's necessary)

hb = 25;            % height of BS (meter)
hu = 300:100:1000;           % altitude of UAV
hh = 20e3;          % altitude of HAPS

hs = 1110e3;          % altitude of satellite (Starlink-I constellation)
ns = 4425;            % number of satellites
theta_0 = 15*pi/180;  % minimum elevation angle of receiver antenna (radian)

nb = 37;              % number of cells (BSs)
nb_JT = 3;            % size of CoMP cluster
Pr_intrf = 0.05;       % Interference probability
ISD = 500;            % inter-site distance in UMa (m)
Rnet = 1.6e3;         % area for UAV distribution = sqrt(2*sqrt(3)/pi)*3*ISD (approximately)
nu = 10;              % number of UAVs
nu_r = 3;             % size of swarm of coordinated UAVs
r_bh = 5e3;           % 2D distance of the HAP from its ground Station
d_2D = 150;           % 2D distance of the UAV from the serving BS (cell radius = sqrt(2*sqrt(3)/pi)*ISD/2 = 260 m)

xyz_HAP = [0,0,hh];
xyz_sat = [0,0,hs];
xyz_gat = [300e3,0,hb];  % suburban

E_th = 1e-5;       % error threshold
D_th = 20e-3;      % delay threshold
n_smp = 1e3;       % number of channel generation
n_itr = 2e3;       % number of iterations
Eb = 1e-6;         % backhual failure probability
Db = 1e-3;         % backhual delay
Db_c = 0.1e-3;     % backhual delay of CoMP
Df = 5e-3;         % feedback delay of CoMP

[Eq_b, Dq_b] = queue(eff_bw, Dq_max, lambda); % queueing delay violation
[Eq_u, Dq_u] = queue(eff_bw, Dq_max, lambda);
[Eq_h, Dq_h] = queue(eff_bw, Dq_max, lambda);
[Eq_g, Dq_g] = queue(eff_bw, Dq_max, lambda);
[Eq_s, Dq_s] = queue(eff_bw, Dq_max, lambda);

Dc = Df + 2*Db_c;              % CoMP delay
Ec = Eb + (1-Eb).* Eb^nb_JT;   % CoMP error

Gr_u = Gr_u - NF_u;
Gr_u = 10^(Gr_u/10);
Gt_u = 10^(Gt_u/10);  % omnidirectional
Gr_h = Gr_h - NF_h;
Gr_h = 10^(Gr_h/10);
Gr_s = Gr_s - NF_s;
Gr_s = 10^(Gr_s/10);
Gt_g = 10^(Gt_g/10);

Rc = sqrt(2*sqrt(3)/pi)*ISD/2;  % cell radius

ii = 1;
for hu = 300:100:1000
    
Dtot_bu = [];
Etot_bu = [];
Dtot_Kuu = [];
Etot_Kuu = [];

Dtot_bu_JT = [];
Etot_bu_JT = [];
Dtot_Kuu_JT = [];
Etot_Kuu_JT = [];

Dtot_hu = [];
Etot_hu = [];

Dtot_su_S = [];
Etot_su_S = [];

Dtot_su_Ka = [];
Etot_su_Ka = [];

for jj = 1:n_itr   % average over position
    
% random distribution of UAVs
[ xyz_BS, xyz_UAV ] = position( nb, hb, ISD, Rnet, nu, hu, d_2D );

d_uu = sqrt(sum((xyz_UAV(2:end,:) - repmat(xyz_UAV(1,:),[nu-1,1]))'.^2));    % distance of UAVs
[d_uu,index] = sort(d_uu,'ascend');
xyz_relay = xyz_UAV(2:end,:);
xyz_UAV = [xyz_UAV(1,:); xyz_relay(index,:)];

% DA2G and A2A communications
[dtot_bu, etot_bu, dtot_Kuu, etot_Kuu] = DA2G_A2A_com(xyz_BS, xyz_UAV, fc_b, fc_u, Nf, BW_bu, BW_uu, R, Pb, Pu, thetaT, ...
    ge_max, Ne, Gt_u, Gr_u, N0, packet_size, hb, hu, nb, nu, nu_r, d_uu, Pr_intrf, n_smp, Eb, Db, Eq_b, Dq_b, Eq_u, Dq_u);

Dtot_bu = [Dtot_bu, dtot_bu];
Etot_bu = [Etot_bu, etot_bu];
Dtot_Kuu = [Dtot_Kuu, dtot_Kuu];
Etot_Kuu = [Etot_Kuu, etot_Kuu];

% Joint transmission CoMP
[dtot_bu_JT, etot_bu_JT, dtot_Kuu_JT, etot_Kuu_JT] = JT_CoMP(xyz_BS, xyz_UAV, fc_b, fc_u, Nf, BW_bu, BW_uu, R, Pb, Pu, thetaT, ge_max, ...
    Ne, Gt_u, Gr_u, N0, packet_size, hb, hu, nb, nb_JT, nu, nu_r, d_uu, Rc, Pr_intrf, n_smp, Eb, Db, Ec, Dc, Eq_b, Dq_b, Eq_u, Dq_u);

Dtot_bu_JT = [Dtot_bu_JT, dtot_bu_JT];
Etot_bu_JT = [Etot_bu_JT, etot_bu_JT];
Dtot_Kuu_JT = [Dtot_Kuu_JT, dtot_Kuu_JT];
Etot_Kuu_JT = [Etot_Kuu_JT, etot_Kuu_JT];

% HAP communication
[dtot_hu, etot_hu] = HAP_com(xyz_BS, xyz_UAV, xyz_HAP, fc_b, fc_h, BW_bh, BW_hu, R, Ph, Pb, thetaT_gat, ge_max, Ne, Gt_hmax, ... 
    Gr_h, Gr_u, N0, packet_size, hb, hu, hh, nb, r_bh, Pr_intrf, n_smp, Eb, Db, Eq_b, Dq_b, Eq_h, Dq_h, c);

Dtot_hu = [Dtot_hu, dtot_hu];
Etot_hu = [Etot_hu, etot_hu];

% Satellite communication
[dtot_su_S, etot_su_S] = Sat_com_Sband(xyz_BS, xyz_gat, xyz_sat, xyz_UAV, fc_s_S, BW_su, BW_gs, R, Pg, Ps, thetaT_gat, Gt_smax, Ne, ...
    Gr_s, Gt_g, Gr_u, Re, theta_0, ns, N0, packet_size, hb, hu, hs, nb, Pr_intrf, n_smp, Eb, Db, Eq_g, Dq_g, Eq_s, Dq_s, c);

Dtot_su_S = [Dtot_su_S, dtot_su_S];
Etot_su_S = [Etot_su_S, etot_su_S];

[dtot_su_Ka, etot_su_Ka] = Sat_com_Kaband(xyz_BS, xyz_gat, xyz_sat, xyz_UAV, fc_s_Ka, BW_su, BW_gs, R, Pg, Ps, thetaT_gat, Gt_smax, Ne, ...
    Gr_s, Gt_g, Gr_u, Re, theta_0, ns, N0, packet_size, hb, hu, hs, nb, Pr_intrf, n_smp, Eb, Db, Eq_g, Dq_g, Eq_s, Dq_s, c);

Dtot_su_Ka = [Dtot_su_Ka, dtot_su_Ka];
Etot_su_Ka = [Etot_su_Ka, etot_su_Ka];

end

%===========================combining schemes==============================
% combining DA2G (ST) and A2A links
Dtot_BU1 = min([Dtot_bu; Dtot_Kuu(1,:)]);  
Etot_BU1 = Etot_bu .* Etot_Kuu(1,:);

Dtot_BU2 = min([Dtot_bu; Dtot_Kuu(2,:)]);
Etot_BU2 = Etot_bu .* Etot_Kuu(2,:);

Dtot_BU3 = min([Dtot_bu; Dtot_Kuu(3,:)]);
Etot_BU3 = Etot_bu .* Etot_Kuu(3,:);

% combining JT and A2A links
Dtot_BU1_JT = min([Dtot_bu_JT; Dtot_Kuu_JT(1,:)]);
Etot_BU1_JT = Etot_bu_JT .* Etot_Kuu_JT(1,:);

Dtot_BU2_JT = min([Dtot_bu_JT; Dtot_Kuu_JT(2,:)]);
Etot_BU2_JT = Etot_bu_JT .* Etot_Kuu_JT(2,:);

Dtot_BU3_JT = min([Dtot_bu_JT; Dtot_Kuu_JT(3,:)]);
Etot_BU3_JT = Etot_bu_JT .* Etot_Kuu_JT(3,:);

% combining HAP and DA2G and A2A links
Dtot_BH = min([Dtot_bu; Dtot_hu]);
Etot_BH = Etot_bu .* Etot_hu;

Dtot_BHU1 = min([Dtot_bu; Dtot_hu; Dtot_Kuu(1,:)]);
Etot_BHU1 = Etot_bu .* Etot_hu .* Etot_Kuu(1,:);

Dtot_BHU2 = min([Dtot_bu; Dtot_hu; Dtot_Kuu(2,:)]);
Etot_BHU2 = Etot_bu .* Etot_hu .* Etot_Kuu(2,:);

Dtot_BHU3 = min([Dtot_bu; Dtot_hu; Dtot_Kuu(3,:)]);
Etot_BHU3 = Etot_bu .* Etot_hu .* Etot_Kuu(3,:);

% combining HAP and JT and A2A links
Dtot_BH_JT = min([Dtot_bu_JT; Dtot_hu]);
Etot_BH_JT = Etot_bu_JT .* Etot_hu;

Dtot_BHU1_JT = min([Dtot_bu_JT; Dtot_hu; Dtot_Kuu_JT(1,:)]);
Etot_BHU1_JT = Etot_bu_JT .* Etot_hu .* Etot_Kuu_JT(1,:);

Dtot_BHU2_JT = min([Dtot_bu_JT; Dtot_hu; Dtot_Kuu_JT(2,:)]);
Etot_BHU2_JT = Etot_bu_JT .* Etot_hu .* Etot_Kuu_JT(2,:);

Dtot_BHU3_JT = min([Dtot_bu_JT; Dtot_hu; Dtot_Kuu_JT(3,:)]);
Etot_BHU3_JT = Etot_bu_JT .* Etot_hu .* Etot_Kuu_JT(3,:);

% combining Sat-Sband and DA2G and A2A links
Dtot_BS_S = min([Dtot_bu; Dtot_su_S]);
Etot_BS_S = Etot_bu .* Etot_su_S;

Dtot_BSU1_S = min([Dtot_bu; Dtot_su_S; Dtot_Kuu(1,:)]);
Etot_BSU1_S = Etot_bu .* Etot_su_S .* Etot_Kuu(1,:);

Dtot_BSU2_S = min([Dtot_bu; Dtot_su_S; Dtot_Kuu(2,:)]);
Etot_BSU2_S = Etot_bu .* Etot_su_S .* Etot_Kuu(2,:);

Dtot_BSU3_S = min([Dtot_bu; Dtot_su_S; Dtot_Kuu(3,:)]);
Etot_BSU3_S = Etot_bu .* Etot_su_S .* Etot_Kuu(3,:);

% combining Sat-Sband and JT and A2A links
Dtot_BS_S_JT = min([Dtot_bu_JT; Dtot_su_S]);
Etot_BS_S_JT = Etot_bu_JT .* Etot_su_S;

Dtot_BSU1_S_JT = min([Dtot_bu_JT; Dtot_su_S; Dtot_Kuu_JT(1,:)]);
Etot_BSU1_S_JT = Etot_bu_JT .* Etot_su_S .* Etot_Kuu_JT(1,:);

Dtot_BSU2_S_JT = min([Dtot_bu_JT; Dtot_su_S; Dtot_Kuu_JT(2,:)]);
Etot_BSU2_S_JT = Etot_bu_JT .* Etot_su_S .* Etot_Kuu_JT(2,:);

Dtot_BSU3_S_JT = min([Dtot_bu_JT; Dtot_su_S; Dtot_Kuu_JT(3,:)]);
Etot_BSU3_S_JT = Etot_bu_JT .* Etot_su_S .* Etot_Kuu_JT(3,:);

% combining Sat-Kaband and DA2G and A2A links
Dtot_BS_Ka = min([Dtot_bu; Dtot_su_Ka]);
Etot_BS_Ka = Etot_bu .* Etot_su_Ka;

Dtot_BSU1_Ka = min([Dtot_bu; Dtot_su_Ka; Dtot_Kuu(1,:)]);
Etot_BSU1_Ka = Etot_bu .* Etot_su_Ka .* Etot_Kuu(1,:);

Dtot_BSU2_Ka = min([Dtot_bu; Dtot_su_Ka; Dtot_Kuu(2,:)]);
Etot_BSU2_Ka = Etot_bu .* Etot_su_Ka .* Etot_Kuu(2,:);

Dtot_BSU3_Ka = min([Dtot_bu; Dtot_su_Ka; Dtot_Kuu(3,:)]);
Etot_BSU3_Ka = Etot_bu .* Etot_su_Ka .* Etot_Kuu(3,:);

% combining Sat-Kaband and JT and A2A links
Dtot_BS_Ka_JT = min([Dtot_bu_JT; Dtot_su_Ka]);
Etot_BS_Ka_JT = Etot_bu_JT .* Etot_su_Ka;

Dtot_BSU1_Ka_JT = min([Dtot_bu_JT; Dtot_su_Ka; Dtot_Kuu_JT(1,:)]);
Etot_BSU1_Ka_JT = Etot_bu_JT .* Etot_su_Ka .* Etot_Kuu_JT(1,:);

Dtot_BSU2_Ka_JT = min([Dtot_bu_JT; Dtot_su_Ka; Dtot_Kuu_JT(2,:)]);
Etot_BSU2_Ka_JT = Etot_bu_JT .* Etot_su_Ka .* Etot_Kuu_JT(2,:);

Dtot_BSU3_Ka_JT = min([Dtot_bu_JT; Dtot_su_Ka; Dtot_Kuu_JT(3,:)]);
Etot_BSU3_Ka_JT = Etot_bu_JT .* Etot_su_Ka .* Etot_Kuu_JT(3,:);

%=================== Network availability calculations ====================
Net_avail_bu(ii) = sum((Etot_bu <= E_th).*(Dtot_bu <= D_th))/(n_smp*n_itr);
Net_avail_Kuu(:,ii) = sum((Etot_Kuu <= E_th).*(Dtot_Kuu <= D_th),2)/(n_smp*n_itr);
Net_avail_BU1(ii) = sum((Etot_BU1 <= E_th).*(Dtot_BU1 <= D_th))/(n_smp*n_itr);
Net_avail_BU2(ii) = sum((Etot_BU2 <= E_th).*(Dtot_BU2 <= D_th))/(n_smp*n_itr);
Net_avail_BU3(ii) = sum((Etot_BU3 <= E_th).*(Dtot_BU3 <= D_th))/(n_smp*n_itr);

Net_avail_bu_JT(ii) = sum((Etot_bu_JT <= E_th).*(Dtot_bu_JT <= D_th))/(n_smp*n_itr);
Net_avail_Kuu_JT(:,ii) = sum((Etot_Kuu_JT <= E_th).*(Dtot_Kuu_JT <= D_th),2)/(n_smp*n_itr);
Net_avail_BU1_JT(ii) = sum((Etot_BU1_JT <= E_th).*(Dtot_BU1_JT <= D_th))/(n_smp*n_itr);
Net_avail_BU2_JT(ii) = sum((Etot_BU2_JT <= E_th).*(Dtot_BU2_JT <= D_th))/(n_smp*n_itr);
Net_avail_BU3_JT(ii) = sum((Etot_BU3_JT <= E_th).*(Dtot_BU3_JT <= D_th))/(n_smp*n_itr);

Net_avail_hu(ii) = sum((Etot_hu <= E_th).*(Dtot_hu <= D_th))/(n_smp*n_itr);
Net_avail_BH(ii) = sum((Etot_BH <= E_th).*(Dtot_BH <= D_th))/(n_smp*n_itr);
Net_avail_BHU1(ii) = sum((Etot_BHU1 <= E_th).*(Dtot_BHU1 <= D_th))/(n_smp*n_itr);
Net_avail_BHU2(ii) = sum((Etot_BHU2 <= E_th).*(Dtot_BHU2 <= D_th))/(n_smp*n_itr);
Net_avail_BHU3(ii) = sum((Etot_BHU3 <= E_th).*(Dtot_BHU3 <= D_th))/(n_smp*n_itr);

Net_avail_BH_JT(ii) = sum((Etot_BH_JT <= E_th).*(Dtot_BH_JT <= D_th))/(n_smp*n_itr);
Net_avail_BHU1_JT(ii) = sum((Etot_BHU1_JT <= E_th).*(Dtot_BHU1_JT <= D_th))/(n_smp*n_itr);
Net_avail_BHU2_JT(ii) = sum((Etot_BHU2_JT <= E_th).*(Dtot_BHU2_JT <= D_th))/(n_smp*n_itr);
Net_avail_BHU3_JT(ii) = sum((Etot_BHU3_JT <= E_th).*(Dtot_BHU3_JT <= D_th))/(n_smp*n_itr);

Net_avail_su_S(ii) = sum((Etot_su_S <= E_th).*(Dtot_su_S <= D_th))/(n_smp*n_itr);
Net_avail_BS_S(ii) = sum((Etot_BS_S <= E_th).*(Dtot_BS_S <= D_th))/(n_smp*n_itr);
Net_avail_BSU1_S(ii) = sum((Etot_BSU1_S <= E_th).*(Dtot_BSU1_S <= D_th))/(n_smp*n_itr);
Net_avail_BSU2_S(ii) = sum((Etot_BSU2_S <= E_th).*(Dtot_BSU2_S <= D_th))/(n_smp*n_itr);
Net_avail_BSU3_S(ii) = sum((Etot_BSU3_S <= E_th).*(Dtot_BSU3_S <= D_th))/(n_smp*n_itr);

Net_avail_BS_S_JT(ii) = sum((Etot_BS_S_JT <= E_th).*(Dtot_BS_S_JT <= D_th))/(n_smp*n_itr);
Net_avail_BSU1_S_JT(ii) = sum((Etot_BSU1_S_JT <= E_th).*(Dtot_BSU1_S_JT <= D_th))/(n_smp*n_itr);
Net_avail_BSU2_S_JT(ii) = sum((Etot_BSU2_S_JT <= E_th).*(Dtot_BSU2_S_JT <= D_th))/(n_smp*n_itr);
Net_avail_BSU3_S_JT(ii) = sum((Etot_BSU3_S_JT <= E_th).*(Dtot_BSU3_S_JT <= D_th))/(n_smp*n_itr);

Net_avail_su_Ka(ii) = sum((Etot_su_Ka <= E_th).*(Dtot_su_Ka <= D_th))/(n_smp*n_itr);
Net_avail_BS_Ka(ii) = sum((Etot_BS_Ka <= E_th).*(Dtot_BS_Ka <= D_th))/(n_smp*n_itr);
Net_avail_BSU1_Ka(ii) = sum((Etot_BSU1_Ka <= E_th).*(Dtot_BSU1_Ka <= D_th))/(n_smp*n_itr);
Net_avail_BSU2_Ka(ii) = sum((Etot_BSU2_Ka <= E_th).*(Dtot_BSU2_Ka <= D_th))/(n_smp*n_itr);
Net_avail_BSU3_Ka(ii) = sum((Etot_BSU3_Ka <= E_th).*(Dtot_BSU3_Ka <= D_th))/(n_smp*n_itr);

Net_avail_BS_Ka_JT(ii) = sum((Etot_BS_Ka_JT <= E_th).*(Dtot_BS_Ka_JT <= D_th))/(n_smp*n_itr);
Net_avail_BSU1_Ka_JT(ii) = sum((Etot_BSU1_Ka_JT <= E_th).*(Dtot_BSU1_Ka_JT <= D_th))/(n_smp*n_itr);
Net_avail_BSU2_Ka_JT(ii) = sum((Etot_BSU2_Ka_JT <= E_th).*(Dtot_BSU2_Ka_JT <= D_th))/(n_smp*n_itr);
Net_avail_BSU3_Ka_JT(ii) = sum((Etot_BSU3_Ka_JT <= E_th).*(Dtot_BSU3_Ka_JT <= D_th))/(n_smp*n_itr);

%========================= mean error analysis ============================
mean_Etot_bu(ii) = mean(Etot_bu);
mean_Etot_uu(ii) = mean(Etot_Kuu(1,:));
mean_Etot_BU1(ii) = mean(Etot_BU1);
mean_Etot_BU2(ii) = mean(Etot_BU2);
mean_Etot_BU3(ii) = mean(Etot_BU3);

mean_Etot_bu_JT(ii) = mean(Etot_bu_JT);
mean_Etot_uu_JT(ii) = mean(Etot_Kuu_JT(1,:));
mean_Etot_BU1_JT(ii) = mean(Etot_BU1_JT);
mean_Etot_BU2_JT(ii) = mean(Etot_BU2_JT);
mean_Etot_BU3_JT(ii) = mean(Etot_BU3_JT);

mean_Etot_hu(ii) = mean(Etot_hu);
mean_Etot_BH(ii) = mean(Etot_BH);
mean_Etot_BHU1(ii) = mean(Etot_BHU1);
mean_Etot_BHU2(ii) = mean(Etot_BHU2);
mean_Etot_BHU3(ii) = mean(Etot_BHU3);

mean_Etot_BH_JT(ii) = mean(Etot_BH_JT);
mean_Etot_BHU1_JT(ii) = mean(Etot_BHU1_JT);
mean_Etot_BHU2_JT(ii) = mean(Etot_BHU2_JT);
mean_Etot_BHU3_JT(ii) = mean(Etot_BHU3_JT);

mean_Etot_su_S(ii) = mean(Etot_su_S);
mean_Etot_BS_S(ii) = mean(Etot_BS_S);
mean_Etot_BSU1_S(ii) = mean(Etot_BSU1_S);
mean_Etot_BSU2_S(ii) = mean(Etot_BSU2_S);
mean_Etot_BSU3_S(ii) = mean(Etot_BSU3_S);

mean_Etot_BS_S_JT(ii) = mean(Etot_BS_S_JT);
mean_Etot_BSU1_S_JT(ii) = mean(Etot_BSU1_S_JT);
mean_Etot_BSU2_S_JT(ii) = mean(Etot_BSU2_S_JT);
mean_Etot_BSU3_S_JT(ii) = mean(Etot_BSU3_S_JT);

mean_Etot_su_Ka(ii) = mean(Etot_su_Ka);
mean_Etot_BS_Ka(ii) = mean(Etot_BS_Ka);
mean_Etot_BSU1_Ka(ii) = mean(Etot_BSU1_Ka);
mean_Etot_BSU2_Ka(ii) = mean(Etot_BSU2_Ka);
mean_Etot_BSU3_Ka(ii) = mean(Etot_BSU3_Ka);

mean_Etot_BS_Ka_JT(ii) = mean(Etot_BS_Ka_JT);
mean_Etot_BSU1_Ka_JT(ii) = mean(Etot_BSU1_Ka_JT);
mean_Etot_BSU2_Ka_JT(ii) = mean(Etot_BSU2_Ka_JT);
mean_Etot_BSU3_Ka_JT(ii) = mean(Etot_BSU3_Ka_JT);

%========================== Error calculation =============================
Err_bu(ii) = sum((Dtot_bu > D_th))/(n_smp*n_itr);
Err_uu(ii) = sum((Dtot_Kuu(1,:) > D_th))/(n_smp*n_itr);
Err_BU1(ii) = sum((Dtot_BU1 > D_th))/(n_smp*n_itr);
Err_BU2(ii) = sum((Dtot_BU2 > D_th))/(n_smp*n_itr);
Err_BU3(ii) = sum((Dtot_BU3 > D_th))/(n_smp*n_itr);

Err_bu_JT(ii) = sum((Dtot_bu_JT > D_th))/(n_smp*n_itr);
Err_uu_JT(ii) = sum((Dtot_Kuu_JT(1,:) > D_th))/(n_smp*n_itr);
Err_BU1_JT(ii) = sum((Dtot_BU1_JT > D_th))/(n_smp*n_itr);
Err_BU2_JT(ii) = sum((Dtot_BU2_JT > D_th))/(n_smp*n_itr);
Err_BU3_JT(ii) = sum((Dtot_BU3_JT > D_th))/(n_smp*n_itr);

Err_hu(ii) = sum((Dtot_hu > D_th))/(n_smp*n_itr);
Err_BH(ii) = sum((Dtot_BH > D_th))/(n_smp*n_itr);
Err_BHU1(ii) = sum((Dtot_BHU1 > D_th))/(n_smp*n_itr);
Err_BHU2(ii) = sum((Dtot_BHU2 > D_th))/(n_smp*n_itr);
Err_BHU3(ii) = sum((Dtot_BHU3 > D_th))/(n_smp*n_itr);

Err_su_S(ii) = sum((Dtot_su_S > D_th))/(n_smp*n_itr);
Err_BS_S(ii) = sum((Dtot_BS_S > D_th))/(n_smp*n_itr);
Err_BSU1_S(ii) = sum((Dtot_BSU1_S > D_th))/(n_smp*n_itr);
Err_BSU2_S(ii) = sum((Dtot_BSU2_S > D_th))/(n_smp*n_itr);
Err_BSU3_S(ii) = sum((Dtot_BSU3_S > D_th))/(n_smp*n_itr);

Err_su_Ka(ii) = sum((Dtot_su_Ka > D_th))/(n_smp*n_itr);
Err_BS_Ka(ii) = sum((Dtot_BS_Ka > D_th))/(n_smp*n_itr);
Err_BSU1_Ka(ii) = sum((Dtot_BSU1_Ka > D_th))/(n_smp*n_itr);
Err_BSU2_Ka(ii) = sum((Dtot_BSU2_Ka > D_th))/(n_smp*n_itr);
Err_BSU3_Ka(ii) = sum((Dtot_BSU3_Ka > D_th))/(n_smp*n_itr);

ii = ii+1;
end
