function [ dtot_bu, etot_bu, dtot_Kuu, etot_Kuu ] = DA2G_A2A_com( xyz_BS, xyz_UAV, fc_b, fc_u, Nf, BW_bu, BW_uu, R, Pb, Pu, thetaT, ...
    ge_max, Ne, Gt_u, Gr_u, N0, packet_size, hb, hu, nb, nu, nu_r, d_uu, Pr_intrf, n_smp, Eb, Db, Eq_b, Dq_b, Eq_u, Dq_u )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

R_bu = R; R_uu = R;

% BU and UU links
for i = 1:nu_r+1
    d_bu(i,:) = sqrt(sum((xyz_BS - repmat(xyz_UAV(i,:),[nb,1]))'.^2));    % 3D distance of BSs and UAVs
    r_bu(i,:) = sqrt(sum((xyz_BS(:,[1 2]) - repmat(xyz_UAV(i,[1 2]),[nb,1]))'.^2));    % 2D horizontal distance of BSs and UAVs
end

d_bu = sort(d_bu,2,'ascend');
r_bu = sort(r_bu,2,'ascend');
theta_bu = atan((hu-hb)./r_bu);

PL_bu_los = 28 + 20*log10(fc_b) + 22*log10(d_bu);
PL_bu_nlos = -17.5 + 20*log10(40*pi*fc_b/3) + (46 - 7*log10(hu))*log10(d_bu);

a1 = 0.3; a2 = 500; a3 = 20;   % for an urban scenario
K = floor(r_bu*sqrt(a1*a2)/1000-1);
for i = 1:nu_r+1
    for j = 1:nb
        Pr_los = 1;
        for k = 0:K(i,j)
            Pr_elem = 1-exp(-(hb-(k+0.5)*(hb-hu)/(K(i,j)+1))^2/(2*a3^2));
            Pr_los = Pr_elem*Pr_los;
        end
        Pr_los_bu(i,j) = Pr_los;
    end
end

PL_bu_dB = Pr_los_bu.*PL_bu_los + (1-Pr_los_bu).*PL_bu_nlos;
PL_uu_dB = 32.45 + 20*log10(fc_u) + 20*log10(d_uu);

PL_bu = 10.^(PL_bu_dB/10);
PL_uu = 10.^(PL_uu_dB/10);

PL_bu = repmat(PL_bu,[Nf,1]);

% rician channel generation
% rice factor of BU channel
for i = 1:nu_r+1
    for j = 1:nb
        if theta_bu(i,j) < 10*pi/180
            K_bu(i,j) = 10^(4/10);
        elseif theta_bu(i,j) >= 10*pi/180 && theta_bu(i,j) < 20*pi/180
            K_bu(i,j) = 10^(5/10);
        elseif theta_bu(i,j) >= 20*pi/180 && theta_bu(i,j) < 30*pi/180
            K_bu(i,j) = 10^(6/10);
        elseif theta_bu(i,j) >= 30*pi/180 && theta_bu(i,j) < 40*pi/180
            K_bu(i,j) = 10^(7/10);
        elseif theta_bu(i,j) >= 40*pi/180 && theta_bu(i,j) < 50*pi/180
            K_bu(i,j) = 10^(8/10);
        elseif theta_bu(i,j) >= 50*pi/180 && theta_bu(i,j) < 60*pi/180
            K_bu(i,j) = 10^(9/10);
        elseif theta_bu(i,j) >= 60*pi/180 && theta_bu(i,j) < 70*pi/180
            K_bu(i,j) = 10^(10/10);
        elseif theta_bu(i,j) >= 70*pi/180 && theta_bu(i,j) < 80*pi/180
            K_bu(i,j) = 10^(11/10);
        else
            K_bu(i,j) = 10^(12/10);
        end
    end
end
K_uu = 10^(12/10);    % rice factor of UU channel

sigma_bu = 1./sqrt(2*(K_bu+1));  % NLoS
sigma_uu = 1./sqrt(2*(K_uu+1));  % NLoS

rho_bu = sqrt(K_bu./(K_bu+1));   % LoS
rho_uu = sqrt(K_uu./(K_uu+1));   % LoS

sigma_bu = repmat(sigma_bu,[Nf,1]);
rho_bu = repmat(rho_bu,[Nf,1]);

for i = 1:Nf*(nu_r+1)
    for j = 1:nb
        h_bu(i,j,:) = rho_bu(i,j) + sigma_bu(i,j)*(randn(1,n_smp) + 1i*randn(1,n_smp));
    end
end
h_uu = rho_uu + sigma_uu*(randn(nu-1,n_smp) + 1i*randn(nu-1,n_smp));

g_bu = h_bu.*conj(h_bu);   % small-scale channel gain
g_uu = h_uu.*conj(h_uu);

% BS directional antenna
thetaZ_bu = pi/2 - theta_bu;    % zenith angle
ge = 10^(ge_max/10)*sin(thetaZ_bu).^2;     % element radiation pattern
gA = sin(Ne*pi*(cos(thetaZ_bu)-cos(thetaT))/2).^2./(Ne*sin(pi*(cos(thetaZ_bu)-cos(thetaT))/2).^2);     % array radiation pattern
Gt_b = ge.*gA;    % linear
Gt_b = repmat(Gt_b,[Nf,1]);

ICI = rand(1,nb);
I = ICI <= Pr_intrf; 
for i = 1:Nf*(nu_r+1)
    for j = 1:nb
        rx_pow_bu(i,j,:) = Pb*Gt_b(i,j)*Gr_u*g_bu(i,j,:)/PL_bu(i,j);
        BS_interf(i,j,:) = I(j)*rx_pow_bu(i,j,:);
    end
end

for i = 1:Nf*(nu_r+1)
    sinr_bu(i,:) = rx_pow_bu(i,1,:)./(BW_bu*N0 + sum(BS_interf(i,2:end,:),2));    % ST scheme
end

IUI = rand(1,nu-1);
II = IUI <= Pr_intrf; 
for j = 1:nu-1
    rx_pow_uu(j,:) = Pu*Gt_u*Gr_u*g_uu(j,:)/PL_uu(j);
    UAV_interf(j,:) = II(j)*rx_pow_uu(j,:);
end

for i = 1:nu_r
    sinr_uu(i,:) = rx_pow_uu(i,:)./(BW_uu*N0 + sum(UAV_interf(nu_r+2:end,:),1));    % nu_r parallel UAVs with JD
end
% sinr_uu2 = sum(rx_pow_uu(1:2,:))./(BW_uu*N0 + sum(rx_pow_uu(3:end,:),1));   % 2 parallel UAVs with MRC

V_bu = 1-(1+sinr_bu).^(-2);
V_uu = 1-(1+sinr_uu).^(-2);
% V_uu2 = 1-(1+sinr_uu2).^(-2);

Dt_bu = packet_size/R_bu;  % transmission delay
Dt_uu = packet_size/R_uu;

et_bu = qfunc((BW_bu*log2(1+sinr_bu)-R_bu)*log(2)./sqrt(BW_bu*V_bu/Dt_bu));  % decoding error
et_uu = qfunc((BW_uu*log2(1+sinr_uu)-R_uu)*log(2)./sqrt(BW_uu*V_uu/Dt_uu));
% et_uu2 = qfunc((BW_uu*log2(1+sinr_uu2)-R_uu)*log(2)./sqrt(BW_uu*V_uu2/Dt_uu));

if Nf == 1
    dtot_bu = Db + Dq_b + Dt_bu./(1-et_bu(1,:));
    etot_bu = 1-(1-Eb).*(1-Eq_b).*(1-et_bu(1,:));
    
    dkol_uu = Db + Dq_b + Dt_bu./(1-et_bu(2:end,:)) + Dq_u + Dt_uu./(1-et_uu);
    ekol_uu = 1-(1-Eb).*(1-Eq_b).*(1-et_bu(2:end,:)).*(1-Eq_u).*(1-et_uu);
else
    dtot_bu = Db + Dq_b + min(Dt_bu./(1-et_bu([1,nu_r+2],:)));
    etot_bu = 1-(1-Eb).*(1-Eq_b).*(1-et_bu(1,:).*et_bu(nu_r+2,:));
    
    dkol_uu = Db + Dq_b + Dt_bu./(1-[min(et_bu([2,6],:));min(et_bu([3,7],:));min(et_bu([4,8],:))]) + Dq_u + Dt_uu./(1-et_uu);
    ekol_uu = 1-(1-Eb).*(1-Eq_b).*(1-et_bu(2:nu_r+1,:).*et_bu(nu_r+3:end,:)).*(1-Eq_u).*(1-et_uu);
end 

% combining A2A links
for k = 1:nu_r  
    dtot_Kuu(k,:) = min(dkol_uu(1:k,:),[],1);
    etot_Kuu(k,:) = prod(ekol_uu(1:k,:),1);  % since each UAV has a different serving BS
end

end

