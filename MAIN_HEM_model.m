%% Hybrid network model of excitotoxicity
function [deda,dDA,kid,simtime,srnd,VtrajectorySTN,VtrajectoryGPE,VtrajectorySNC,Ca_trajectorySNC,mCatrajectorySNC,Ttime,Nstn,erCatrajectorySNC,Ca_er,Ca_mt,er_catrajectorySNC,mt_catrajectorySNC]=MAIN_HEM_model(dt,durr,peren,wstsn,scfa,apopthr,camtthr,cl,gion,gi_dose,dron,dr_dose,ccbon,ccb_dose,cbdon,cbd_dose,asbon,asb_dose,gpuon)
%% CREDITS
% Created by
% Vignayanandam R. Muddapu (Ph.D. scholar)
% C/o Prof. V. Srinivasa Chakravarthy
% Indian Institute of Technology Madras
% India

% SNc with ATP dynamics (Francis et.al., 2013)
% Dopamine synthesis, storage, release, metabolism and terminal autoreceptors (Bravo, 2012)
% Ca2+ induced apoptosis (Hong et.al., 2012)
% Calcium-induced calcium release (Marhl et.al., 2000)
% Energy Metabolism (Cloutier & Wellstead, 2010)
% STN, GPe - (Alekhya et.al., 2015)

%%

% clear;clc;
% time=clock;curdate=time(3);curmonth=time(2);
% peren=[20];
% stsn=[5.5];
% scfa=0.000005;
% apopthr=[0.8];
% durr=50000;%32500
% % filename='test2';
% wstsn=deci2str(stsn);
% peren1=deci2str(peren);
% apopthr1=deci2str(apopthr);
% scfa1=deci2str(scfa);
% %     filename=strcat('HybMod_SSG_GDA_stnlatRS_DAr0-5--15e-6__nRRP10_sncstn1_Rstnsnc5-5_',num2str(durr),'msec_',num2str(curdate),'_',num2str(curmonth));
% filename=strcat('H_MAp-1_V0_0-4_t',num2str(apopthr1),'_PD',num2str(peren1),'_allr_rA0-1_lamR10e-6-5e-6_Fexp6_DA1-5_3-5e-5_Rstnsnc',num2str(wstsn),'_std0-5_scfa',num2str(scfa1),'_1_',num2str(durr/1000),'sec_',num2str(curdate),'_',num2str(curmonth));

tic;

STNcells2track = [2, 4; 1, 8; 4, 3;];
%GPEcells2track = ...
%SCNcells2track = ...

% time1=clock;
% Time parameters & Random seeding
taustn=dt;
taugpe=dt;
tspan=dt:dt:durr;
Ttime=numel(tspan);
srnd = rng;
sfcaiapop=1;
k11f=0.1; % (muM*sec)-1

%%
% Energy deficit conditions
% nATP=ones(8,8);
% xmins=0.5;xmaxs=1;
% nATP=xmins+rand(8,8)*(xmaxs-xmins);

% ATP1=ATP0.*ones(8,8);
% tatp=numel(ATP1);
% peratp=round(peren*tatp/100);
% idx = randsample(1:tatp,tatp);
% pratp=idx(1:peratp);
% endf=10000/dt;
% rATP1=0.1;

%

pd=1;
% stsn=5.5; %5.5
% datp1=deci2str(datp);stsn1=deci2str(stsn);

% filename=strcat('HybMod_SSG_nRRP10_sncstn',num2str(pd),'_Rstnsnc',num2str(stsn1),'-sd-0-1_rATP1-datp',num2str(datp1),'_',num2str(durr),'msec_',num2str(curdate),'_',num2str(curmonth));

%%
%STN
nSTN=32; % (nSTNxnSTN network size)
Mstn=nSTN;
Nstn=nSTN;
Pstn=Mstn*Nstn;

%SNc
nSNc=8; % (nSNcxnSNc network size)
Msnc=nSNc;
Nsnc=nSNc;
Psnc=Msnc*Nsnc;

%GPe
nGPe=32; % (nGPexnGPe network size)
Mgpe=nGPe;
Ngpe=nGPe;
Pgpe=Mgpe*Ngpe;

% Neuron properties
%STN
astn=0.005; % (1/ms)
bstn=0.265; % (1/mV)
cstn=-65; % (mV)
dstn=1.5; % 2-Thibeault2013
vpeak_stn=30;
% taustn = 0.1;

%GPe
agpe=0.1; % (1/ms)
bgpe=0.2; % (1/mV)
cgpe=-65; % (mV)
dgpe=2;
vpeak_gpe=30;
% taugpe=0.1;

% Membrane capacitances
Cstn = 1; %(microF)
Cgpe = 1; %(microF)

% V, U initialization
%STN
Vstn = -62.5.*(rand(Mstn,Nstn)-0.5.*ones(Mstn,Nstn));
Ustn = ((-15)-(-5)).*rand(Mstn,Nstn) + (-5);
VtrajectorySTN = -62.5.*ones(Ttime,Nstn); 
mCatrajectorySNC = zeros(Ttime, 5);
erCatrajectorySNC = zeros(Ttime, 5);
er_catrajectorySNC = zeros(Ttime, 8);
mt_catrajectorySNC = zeros(Ttime, 8);
% Vstn = -62.5.*ones(Mstn,Nstn);
% Ustn = zeros(size(Vstn));
% Vstns=Vstn;Ustns=Ustn;

%GPe
Vgpe = -53.67.*(rand(Mgpe,Ngpe)-0.5.*ones(Mgpe,Ngpe));
Ugpe = ((-15)-(-5)).*rand(Mgpe,Ngpe) + (-5);
% Vgpe = -53.67.*ones(Mgpe,Ngpe);
% Ugpe = zeros(size(Vgpe));
% Vgpes=Vgpe;Ugpes=Ugpe;

% Currents initilization
%STN
Istn=3*ones(size(Vstn)); % 10 Hz->1.9
stn_zeros = zeros(size(Vstn));
stncurr_spk = stn_zeros;
stn_firings=[];stn_firings2=[];
% spkhststn=[];

%GPe
Igpe=4.25*ones(size(Vgpe));
gpe_zeros = zeros(size(Vgpe));
gpecurr_spk = gpe_zeros;
gpe_firings=[];gpe_firings2=[];

% psp variable initilization
%STN
h_nmdastn=stn_zeros;
h_ampastn=stn_zeros;
h_gs = stn_zeros;

%GPe
h_nmdagpe=gpe_zeros;
h_ampagpe=gpe_zeros;
hlat_gaba_gpe=gpe_zeros;

% decay constants(ms)
taunmda=160;
tauampa=6;
taugaba=4;

% dt/T in PSP
lam_nmda = dt/taunmda;
lam_ampa = dt/tauampa;
lam_gaba= dt/taugaba;

mg0=1; % magnesium conc.

% RMP of receptors
Enmda = 0;
Eampa = 0;
Egaba = -60;

xmin=-55;xmax=-45;
V_sncinit = xmin+rand(Msnc,Nsnc)*(xmax-xmin);
% V_sncinit = -49.42.*ones(Msnc,Nsnc);
Ca_iinit = 0.000188.*ones(Msnc,Nsnc);

Ca_trajectorySNC = 0.000188.*ones(1,Nsnc);

Na_iinit = 4.6876.*ones(Msnc,Nsnc);
K_iinit = 126.05893.*ones(Msnc,Nsnc);
Calbinit = 0.0026.*ones(Msnc,Nsnc);
Caminit = 0.0222.*ones(Msnc,Nsnc);
m_calinit = 0.006271.*ones(Msnc,Nsnc);
m_nainit = 0.0952.*ones(Msnc,Nsnc);
h_nainit = 0.1848.*ones(Msnc,Nsnc);
O_hcninit = 0.003.*ones(Msnc,Nsnc);
m_kdrinit = 0.0932.*ones(Msnc,Nsnc);
y_pcinit = 0.483.*ones(Msnc,Nsnc);
y_nkinit = 0.6213.*ones(Msnc,Nsnc);
Ca_erinit = 1.0*0.001.*ones(Msnc,Nsnc); %mM
Ca_mtinit = 0.4*0.001.*ones(Msnc,Nsnc); % mM % 0.4e-3
cdainit = 1e-4.*ones(Msnc,Nsnc); %mM%1e-4
vdainit = 500.*ones(Msnc,Nsnc); %mM 500
edainit = 4e-6.*ones(Msnc,Nsnc); %mM
Iextinit=0.*ones(Msnc,Nsnc); %Iext
ATPusedinit=0.*ones(Msnc,Nsnc);
calinit=1.*ones(Msnc,Nsnc); %mM
cai_calinit=0.*ones(Msnc,Nsnc); %mM
cal_actinit=0.*ones(Msnc,Nsnc); %mM
casp12init=1.*ones(Msnc,Nsnc); %mM
cal_act_casp12init=0.*ones(Msnc,Nsnc); %mM
casp12_actinit=0.*ones(Msnc,Nsnc); %mM
casp9init=1.*ones(Msnc,Nsnc); %mM
casp12_act_casp9init=0.*ones(Msnc,Nsnc); %mM
casp9_actinit=0.*ones(Msnc,Nsnc); %mM
casp3init=1.*ones(Msnc,Nsnc); %mM
casp9_act_casp3init=0.*ones(Msnc,Nsnc); %mM
casp3_actinit=0.*ones(Msnc,Nsnc); %mM
xmin1=0;xmax1=0.1;
apopinit = xmin1+rand(Msnc,Nsnc)*(xmax1-xmin1);
% apopinit=0.*ones(Msnc,Nsnc); %mM
ROS_mitinit=0.*ones(Msnc,Nsnc);
PTP_mit_actinit=0.*ones(Msnc,Nsnc);
Cytc_mitinit=1.*ones(Msnc,Nsnc);
Cytcinit=0.*ones(Msnc,Nsnc);
Cytc_casp9init=0.*ones(Msnc,Nsnc);
IAPinit=0.*ones(Msnc,Nsnc);
casp9_act_IAPinit=0.*ones(Msnc,Nsnc);
casp3_act_IAPinit=0.*ones(Msnc,Nsnc);
F6Pinit = 0.175883476634895.*ones(Msnc,Nsnc);%0.2
F26Pinit = 0.002191750879602.*ones(Msnc,Nsnc);%0.001
GAPinit = 0.082507126186107.*ones(Msnc,Nsnc);%0.0405
PYRinit = 0.123910489378719.*ones(Msnc,Nsnc);%0.1
LACinit = 0.598605032933119.*ones(Msnc,Nsnc);%0.5
ATPinit = 2.395615876085214.*ones(Msnc,Nsnc);%2.402
PCrinit = 18.044071098085976.*ones(Msnc,Nsnc);%18.14
% apopinit=0.*ones(Msnc,Nsnc); %mM
Ssnc=zeros(Msnc,Nsnc);

ROS_mit=ROS_mitinit;
PTP_mit_act=PTP_mit_actinit;
Cytc_mit=Cytc_mitinit;
Cytc=Cytcinit;
Cytc_casp9=Cytc_casp9init;
IAP=IAPinit;
casp9_act_IAP=casp9_act_IAPinit;
casp3_act_IAP=casp3_act_IAPinit;
F6P=F6Pinit;
F26P=F26Pinit;
GAP=GAPinit;
PYR=PYRinit;
LAC=LACinit;
ATP=ATPinit;
PCr=PCrinit;

V_snc=V_sncinit;m_cal=m_calinit;m_na=m_nainit;
h_na=h_nainit;O_hcn=O_hcninit;Calb=Calbinit;
Cam=Caminit;y_nk=y_nkinit;y_pc=y_pcinit;m_kdr=m_kdrinit;
K_i=K_iinit;Na_i=Na_iinit;Ca_i=Ca_iinit;Ca_er=Ca_erinit;Ca_mt=Ca_mtinit;

cda=cdainit;vda=vdainit;
eda=edainit;Iexts=Iextinit;ATPused=ATPusedinit;
cal=calinit;cai_cal=cai_calinit;cal_act=cal_actinit;casp12=casp12init;
cal_act_casp12=cal_act_casp12init;casp12_act=casp12_actinit;casp9=casp9init;
casp12_act_casp9=casp12_act_casp9init;casp9_act=casp9_actinit;casp3=casp3init;
casp9_act_casp3=casp9_act_casp3init;casp3_act=casp3_actinit;apop=apopinit;

sim_mM=1e-3;
sim_mM_msec=1e-3/3.6e6;
sim_msec=1/3.6e6;
sim_msec_mM=1/((1e-3)*(3.6e6));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
%Soma
R = 8314.472; %mJ/mol. K
T = 310.15; %K
F = 96485.30929; %coulomb/mol.
Ca_o = 1.8; %mM
Na_o = 137; %mM
K_o = 5.4; %mM
vol_pmu = 5; %pl
fr_cyt = 0.5;
C_sp = 0.9e6; %pF/cm2
SVR_pmu = 1.6667e4; %1/cm
% ATP = 0.411; %mM 0.411 - bursting
Calbtot = 0.005; % mM
Camtot = 0.0235; % mM
kcal_1 = 10; %1/mM.ms
kcal_2 = 2e-3; %1/ms
kcam_cd = 0.003; %1/ms
kcam_nd = 3; %1/ms
g_cal = 2101.2; %pA/mM
g_na = 907.68; %pA/mM
A_mna = 1.9651; %1/ms
B_mna = 0.0424; %1/ms
A_hna = 9.566e-5; %1/ms
B_hna = 0.5296; %1/ms
za_mna = 1.7127;
zb_mna = 1.5581;
za_hna = -2.4317;
zb_hna = -1.1868;
g_nalk = 0.0053; %pA/mM
g_nahcn = 51.1; %pA/mM
cAMP = 1e-5; %mM
g_ksk = 2.2515; %pA/mM
g_kdr = 31.237; %nS
g_kir = 13.816; %nS
k_2pc = 0.001; %1/ms
k_3pc = 0.001; %1/ms
k_4pc = 1; %1/ms
K_pco = 2; %mM
k_pmca = 2.233;
dell = 0.35;
k_xm = 0.0166; %pA
k_2nk = 0.04; %1/ms
k_3nk = 0.01; %1/ms
k_4nk = 0.165; %1/ms
K_nknai = 4.05; %mM
K_nknao = 69.8; %mM
K_nkki = 32.88; %mM
K_nkko = 0.258; %mM
k_nk = 1085.7; %pA
V_tau = (R*T)./F;
vol_cyt = fr_cyt*vol_pmu;
P_c = 1.00000./(1.00000+cAMP./0.00116300);
P_o = 1.00000./(1.00000+cAMP./1.45000e-05);
P_E2Spc = 1.00000./(1.00000+K_pco./Ca_o);
A_pmu = (SVR_pmu.*vol_pmu.*0.00100000.*0.00100000.*0.00100000)./1.00000;
P_E2pc = 1.00000-P_E2Spc;
beta_pc = k_2pc.*P_E2Spc+k_4pc.*P_E2pc;

%CICR
% ER
rho_er = 0.01; % rho_er
beta_er = 0.0025; %beta_er
k_pump = 20.0/1000; % 1/ms
k_ch = 3000.0/1000; %1/ms
K1 = 5.0*0.001; %mM
k_leak = 0.05/1000; %1/ms

% Mito
rho_mt = 0.01; % rho_mt
beta_mt = 0.0025; % beta_mt
k_in = 300.0*0.001/1000; % mM/ms
K2 = 0.8*0.001; % mM
k_out = 125.0/1000; % 1/ms
k_m = 0.00625/1000; % 1/ms
K3 = 5.0*0.001; % mM

% DA terminal (Tello-Bravo (2012))
krel = 0.031; %(mM) 0.055
psi = 17.4391793; %(mM/ms)
nRRP = 10; % :RANGE ~ 10-30
Veda_max = 1e-6; %(mM/ms)
Keda = 3e-5; %(mM)
kcomt = 0.0083511; %(1/ms)
%vda = 500; %(mM)
vdao = 500; %(mM)
vdas = 1e-2; %(mM)
dara = 5e-5; %(mM)
dars = 1e-2; %(mM)
Vsynt_max = 250e-5;%(mM/ms)%30.2e-6 %25e-6
Ksynt = 35e-4; %(mM)
Ktyr = 46e-3; %(mM)
TYR = 126e-3; %(mM)
Kicda = 11e-2; %(mM)
Kieda = 46e-3; %(mM)
Vcda_max = 0.2*133.33e-6; %(mM/ms)%133.33e-6*0.03
Kcda = 238e-4; %(mM)
kmao = 0.00016; %(1/ms)

% Apoptosis pathway (Hong et.al., (2012))
k3f=1; % (muM*sec)-1
k3b=1/1e3; % (sec)-1
k4f=1/1e3; % (sec)-1
k5f=1; % (muM*sec)-1
k5b=1/1e3; % (sec)-1
k6f=1/1e3; % (sec)-1
k7f=10; % (muM*sec)-1
k7b=0.5/1e3; % (sec)-1
k8f=1/1e3; % (sec)-1
k9f=10; % (muM*sec)-1
k9b=0.5/1e3; % (sec)-1
k10f=0.1/1e3; % (sec)-1
k11f=1; % (muM*sec)-1

k29f=0.5; % (mM*msec)-1
k30f=0.5; % (mM*msec)-1
k31f=1; % (mM*msec)-1
k27f=1; % (mM*msec)-1
k27b=1/1e3; % (msec)-1
k28f=1/1e3; % (msec)-1
k12f=5; % (mM*msec)-1
k12b=0.0035/1e3; % (msec)-1
k13f=5; % (mM*msec)-1
k13b=0.0035/1e3; % (msec)-1

Mit=1;
Sig_ers=0;%0.0001;
Sig_mts=0;%0.0001;
PTP_mit=1;

% Energy Metabolism
GLCe=1;%mM
Vmax_hk = 2.5/1000;%mM/ms
Km_ATP_hk = 0.5;%mM
KI_F6P = 0.068;%mM
Vmax_pfk = 3.85/1000;%mM/ms
Km_ATP_pfk = 0.05;%mM
Km_F6P_pfk = 0.18;%mM
Km_F26P_pfk = 0.01;%mM
Vmaxf_pfk2 = 2e-04/1000;%mM/ms
Vmaxr_pfk2 = 1.036e-04/1000;%mM/ms
Km_ATP_pfk2 = 0.05;%mM
Km_F6P_pfk2 = 0.01;%mM
Km_F26P_pfk2 = 0.0001;%mM
Vmax_pk = 5.0/1000;%mM/ms
Km_ADP_pk = 0.005;%mM
Km_GAP_pk = 0.4;%mM
Vmax_op = 1.0/1000;%mM/ms
Km_ADP_op = 0.005;%mM
Km_PYR_op = 0.5;%mM
kf_ldh = 12.5/1000;%1/ms
kr_ldh = 2.5355/1000;%1/ms
kf_ck = 3.0/1000;%1/mM.ms
kr_ck = 1.26/1000;%1/mM.ms
PCr_tot = 20.0;%mM
Vmax_ATPase = 0.9355/1000;%mM/ms
Km_ATP = 0.5;%mM
Vlac_0 = 0.355/1000;%mM/ms
K_lac_eff = 0.71/1000;%1/ms
K_lac = 0.641;
ANP = 2.51;%mM
Q_adk = 0.92;
nATP = 0.4;
KI_ATP = 1.0;%mM
nAMP = 0.5;
Ka_AMP = 0.05;%mM
Kamp_pfk2 = 0.005;%mM
nh_amp = 2;
beta_ldh_ros=0.25;
Kldh_ros=10*sim_mM;%muM
eta_op_max=0.995;
beta_op_asyn=0.08;
Kasyn=8.5*sim_mM; %mM
ROS=0.001;%mM
ASYNA=0.001;%mM

%%
snc_firings=[];snc_firings2=[];
DA=0;rss=1.6;nlat=5;

%STN_SNc projections (no. of STN (projXproj) to no. (1) of SNc)
proj=4;
nproj=proj*proj;
idx = randsample(1:Pstn,Pstn);

% per=CL;
% perkill=round(Psnc*per/100);
% kil = randsample(1:Psnc,perkill);

%SNc
h_nmdasnc=zeros(Msnc,Nsnc);
h_ampasnc=zeros(Msnc,Nsnc);

% STN-SNc connections
stsn=[1];
wstnsnc_matrix=wstsn.*normrnd(wstsn,0.5,8,8);%sd=0.1

% Effect of DA on post-synaptic currents
cd2=0.1;CD2=0.1;

wsg1=1;wgs1=20;stt=1000;

da=0.5;
wlatgpe = weightcal_gpe(da);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Stimulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Istim=0; % Phasic bursting Istim=0.00001
del=1000/dt;
dur=200/dt;
sigg1=0;sigg2=0;phier=0;phimt=0;
enfatra=1;
enfamit=1;
eda1=0;

dDA=[];
% dcai=[];datpused=[];dapop=[];
deda=[];
% dV_snc=[];dcaer=[];dcamt=[];dcam=[];dcalb=[];dros_mit=[];
dnc=[];
% dIstnsnc=[];dIsncsnc=[];dIstnstn=[];dIgpestn=[];dIgpegpe=[];dIstngpe=[];dV_snc=[];datp=[];

sncstart=500;count=1;stt=1000;samp_start=stt/dt;
samp_start=stt/dt;%5000
indsapp=[];
% yeda1=zeros(Psnc,samp_start);
% dapop1=zeros(Psnc,samp_start);
% dnc1=zeros(1,samp_start);
% dDA1=zeros(1,samp_start);
% dros_mit1=zeros(Psnc,samp_start);
% dcai1=zeros(Psnc,samp_start);
% datpused1=zeros(Psnc,samp_start);
% datp1=zeros(Psnc,samp_start);
% Istnsnc1=zeros(Psnc,samp_start);
% Isncsnc1=zeros(Psnc,samp_start);
% Istnstn1=zeros(Pstn,samp_start);
% Igpestn1=zeros(Pstn,samp_start);
% Igpegpe1=zeros(Pgpe,samp_start);
% Istngpe1=zeros(Pgpe,samp_start);
% dV_snc1=zeros(Psnc,samp_start);
% dcaer1=zeros(Psnc,samp_start);
% dcamt1=zeros(Psnc,samp_start);
% dcam1=zeros(Psnc,samp_start);
% dcalb1=zeros(Psnc,samp_start);

yeda1=zeros(1,samp_start);
% dapop1=zeros(Psnc,samp_start);
dnc1=zeros(1,samp_start);
dDA1=zeros(1,samp_start);
% dros_mit1=zeros(Psnc,samp_start);
% dcai1=zeros(Psnc,samp_start);
% datpused1=zeros(Psnc,samp_start);
% datp1=zeros(Psnc,samp_start);
% Istnsnc1=zeros(Psnc,samp_start);
% Isncsnc1=zeros(Psnc,samp_start);
% Istnstn1=zeros(Pstn,samp_start);
% Igpestn1=zeros(1,samp_start);
% Igpegpe1=zeros(1,samp_start);
% Istngpe1=zeros(1,samp_start);
% dV_snc1=zeros(Psnc,samp_start);
% dcaer1=zeros(Psnc,samp_start);
% dcamt1=zeros(Psnc,samp_start);
% dcam1=zeros(Psnc,samp_start);
% dcalb1=zeros(Psnc,samp_start);

NEWstn2snc=1;wDA_snc=1;
nolat=1;counttt=1;

lmin=0.00001;lmax=0.00005;
lam=lmin+rand(8,8)*(lmax-lmin);


%%%SELF-KILLING
idxx = randsample(1:Psnc,Psnc);
sst=durr/dt;
ssp=1;
indsapp=[];

%% Energy deficit
enfatra=1.*ones(Msnc,Nsnc);
enfamit=1.*ones(Msnc,Nsnc);
tatp=numel(enfatra);
peratp=round(peren*tatp/100);
idxED = randsample(1:tatp,tatp);
pratp=idxED(1:peratp);
endf=10000/dt;
counttt=1;
indsappcamt=[];
indsappcaer=[];

Sig_mts=0*ones(Msnc,Nsnc);
Sig_ers=0*ones(Msnc,Nsnc);

lmin=0.00001;lmax=0.00005;
lam=lmin+rand(Msnc,Nsnc)*(lmax-lmin);

% Initiating Therapy after certain percentage cell loss
NcellLoss=round((cl*Psnc)/100); % 50% cell loss

%%%%%% Therrapeutics %%%%%%
% Glutamate Inhibitor
gi=1;
Inc=0;
dr=0;
ccb=1;
cbd=0;
asb=1;
%         enfatra=0.2*ones(Msnc,Nsnc);%
%         enfamit
% clit=cell(1,Ttime);
indsappp=[];
lims=0.2;
k2=1;
%%
for k = 1:Ttime
    
    if gion==1
        if k>k2
            if Inc>=NcellLoss
                gi=gi_dose;
            end
        end
    end
    
    if dron==1
        if k>k2
            if Inc>=NcellLoss
                dr=dr_dose;
            end
        end
    end
    
    if ccbon==1
        if k>k2
            if Inc>=NcellLoss
                ccb=ccb_dose;
            end
        end
    end
    
    if cbdon==1
        if k>k2
            if Inc>=NcellLoss
                cbd=cbd_dose;
            end
        end
    end
    
    if asbon==1
        if k>k2
            if Inc>=NcellLoss
                asb=asb_dose;
            end
        end
    end
    
        ATP=2.4*ones(Msnc,Nsnc);
    
        if k==sst
            indsapp=idxx(1:ssp);
            ssp=ssp+1;
            sst=sst+(200/dt);
        end
    
    V_snc(indsapp) = -80.*ones(size(indsapp));
    apop(indsapp)=zeros(size(indsapp));
    Ca_i(indsapp)=zeros(size(indsapp)); %mM
    eda(indsapp) = 26e-6.*zeros(size(indsapp)); %mM
    ATPused(indsapp)=0.*zeros(size(indsapp));
    ATP(indsapp)=0.*zeros(size(indsapp));
    
    if k>endf
        Renfatra=1.*exp(-counttt.*lam);
        Renfamit=1.*exp(-counttt.*lam);
        %         edda=1e-5.*exp(-counttt.*lam);
        counttt=counttt+1;
        
        enfatra(pratp)=Renfatra(pratp);
        enfamit(pratp)=Renfamit(pratp);
        %         eda(pratp)=edda(pratp);
        
        enfatra(enfatra<lims)=lims;
        enfamit(enfamit<lims)=lims;
        %        nATP(nATP<rATP1)=rATP1;
        
    end
    
    %     if(k > 5000/dt)
    %         enfatra=0.2*ones(Msnc,Nsnc);%
    %         enfamit=0.1*ones(Msnc,Nsnc);;
    %     else
    %         enfatra=1*ones(Msnc,Nsnc);;
    %         enfamit=1*ones(Msnc,Nsnc);;
    %     end
    %
    
    %     phier=Ca_i-Ca_er;
    
    %     inds1=Ca_mt(pratp)>camtthr;
    %     inds2=pratp(inds1);
    
    inds2=find(Ca_mt>camtthr);
    indsappcamt=[indsappcamt inds2];
    if isempty(inds2)==0
        indsappcamt=unique(indsappcamt);
        inds2=[];
    end
    Ca_mt(indsappcamt)=zeros(size(indsappcamt)); %mM
    %
    %     indsapcaer=find(phier>0.0);
    %     indsappcaer=[indsappcaer indsapcaer'];
    %     indsappcaer=unique(indsappcaer);
    %
    Sig_mts(indsappcamt)=0.01.*ones(size(indsappcamt));
    %     Sig_ers(indsappcaer)=0.01.*ones(size(indsappcaer));
    
    DA = RescaleRange(eda1,1e-5,2.1e-5,0,1);
    
    if DA<0
        DA=0;
    end
    if DA>1
        DA=1;
    end
    DA=DA+dr;
    
    wda_gpe=1;
    %     ssmax = RescaleRange(DA,0,1,40,0.1);
    wsg=((1-CD2*da))*wsg1;
    wgs=((1-CD2*DA))*wgs1;
    wlatstn = weightcal_stn(DA*pd);
    
    %----------------------------------------SNc-----------------------------------------%
    % Lateral SNc-SNc connections
    wlatsnc = weightcal_snc(DA,rss,nlat);
    Hsnc=1./(1+exp(-(V_snc-20+57)/2));
    Ssncnxt=Ssnc+((2.*Hsnc.*(1-Ssnc))-Ssnc.*0.08).*dt;
    Ssncfin=conv2(Ssncnxt,wlatsnc,'same');
    Igabasnc=1.*0.01.*(V_snc+63.45).*Ssncfin;
    
    % STN-SNc connections
    if (NEWstn2snc==1)
        start=1;stop=nproj;totspk=zeros(1,Psnc);
        for ll=1:Psnc
            %             stn_snccurr1=zeros(proj,proj);
            %             stn_snccurr1=reshape(stncurr_spk(idx(start:stop)),proj,proj);
            totspk(ll)=sum(stncurr_spk(idx(start:stop)));
            start=start+nproj;stop=stop+nproj;
        end
    end
    stn_snccurr=zeros(8,8);
    stn_snccurr=reshape(totspk,8,8);
    
    h_nmdasnc = (1-lam_nmda).* h_nmdasnc + lam_nmda.*stn_snccurr; %psp nmda snc
    h_ampasnc = (1-lam_ampa).* h_ampasnc + lam_ampa.*stn_snccurr;
    
    tmp_nmda_snc = h_nmdasnc.*(Enmda - V_snc);
    tmp_ampa_snc = h_ampasnc.*(Eampa - V_snc);
    
    
    I_nmda_snc = wDA_snc.*gi.*wstnsnc_matrix.*tmp_nmda_snc;
    I_ampa_snc = wDA_snc.*gi.*wstnsnc_matrix.*tmp_ampa_snc;
    
    B = 1./(1 + (mg0/3.57).*exp(-0.062.*V_snc));
    
    Istnsnc=1*(B.*I_nmda_snc + I_ampa_snc);
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Stimulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Ibg=0;
    if(k < del + dur && k > del)
        Iapp = Ibg +Istim;
    else
        Iapp = Ibg;
    end
    Iext=Iapp;
    
    %Soma
    
    %ATP-dependent equations
    k_1pc = 1.00000./(1.00000+0.100000./ATP);
    k_1nk = 0.370000./(1.00000+0.0940000./ATP);
    
    %Soma
    %Membrane potential
    VD = V_snc./V_tau;
    
    % HCN current
    kf_free = 0.00600000./(1.00000+exp((V_snc+87.7000)./6.45000));
    kf_bnd = 0.0268000./(1.00000+exp((V_snc+94.2000)./13.3000));
    kf_hcn = kf_free.*P_c+kf_bnd.*(1.00000-P_c);
    kr_free = 0.0800000./(1.00000+exp(-(V_snc+51.7000)./7.00000));
    kr_bnd = 0.0800000./(1.00000+exp(-(V_snc+35.5000)./7.00000));
    kr_hcn = kr_free.*P_o+kr_bnd.*(1.00000-P_o);
    
    % Calcium binding proteins
    CaCalb = Calbtot-Calb;
    J_calb = kcal_1.*(Calb+cbd).*Ca_i-kcal_2.*CaCalb;
    
    CaCam = Camtot-Cam;
    kcam_cb = 12000.0.*(power(Ca_i, 2.00000));
    kcam_nb = 3.70000e+06.*(power(Ca_i, 2.00000));
    alpha_cam = kcam_cb.*kcam_nb.*(1.00000./(kcam_cb+kcam_nd)+1.00000./(kcam_cd+kcam_nd));
    beta_cam = kcam_cd.*kcam_nd.*(1.00000./(kcam_cb+kcam_nd)+1.00000./(kcam_cd+kcam_nd));
    J_cam = alpha_cam.*(Cam+cbd)-beta_cam.*CaCam;
    
    K_pci = (173.600./(1.00000+CaCam./5.00000e-05)+6.40000).*1.00000e-05;
    P_E1Spc = 1.00000./(1.00000+K_pci./Ca_i);
    P_E1pc = 1.00000-P_E1Spc;
    alpha_pc = k_1pc.*P_E1Spc+k_3pc.*P_E1pc;
    
    V_Ca = 0.500000.*log(Ca_o./Ca_i);
    h_cal = 0.000450000./(0.000450000+Ca_i);
    I_CaL = ccb.*((g_cal.*m_cal.*h_cal.*(power(Ca_i.*Ca_o, 1.0./2)).*sinh(VD-V_Ca))./(sinh(VD)./VD));
    K_pmca = k_pmca.*((10.5600.*CaCam)./(CaCam+5.00000e-05)+1.20000);
    I_pmca = K_pmca.*(k_1pc.*P_E1Spc.*y_pc-k_2pc.*P_E2Spc.*(1.00000-y_pc)).*1.00000;
    Dr = (1.00000+0.00100000.*((power(Na_i, 3.00000)).*Ca_o+(power(Na_o, 3.00000)).*Ca_i)).*(1.00000+Ca_i./0.00690000);
    I_xm = (k_xm.*((power(Na_i, 3.00000)).*Ca_o.*exp(dell.*VD)-(power(Na_o, 3.00000)).*Ca_i.*exp((dell-1.00000).*VD)))./Dr;
    J_ca = (-1.00000./(2.00000.*F.*vol_cyt)).*((I_CaL+2.00000.*I_pmca)-2.00000.*I_xm);
    
    V_Na = log(Na_o./Na_i);
    O_na = (power(m_na, 3.00000)).*h_na;
    I_Na = (g_na.*O_na.*(power(Na_i.*Na_o, 1.0./2)).*sinh(0.500000.*(VD-V_Na)))./(sinh(0.500000.*VD)./(0.500000.*VD));
    I_Nalk = (g_nalk.*(power(Na_i.*Na_o, 1.0./2)).*sinh(0.500000.*(VD-V_Na)))./(sinh(0.500000.*VD)./(0.500000.*VD));
    I_NaHCN = (g_nahcn.*O_hcn.*(power(Na_i.*Na_o, 1.0./2)).*sinh(0.500000.*(VD-V_Na)))./(sinh(0.500000.*VD)./(0.500000.*VD));
    P_E1Snk = 1.00000./(1.00000+(K_nknai./Na_i).*(1.00000+K_i./K_nkki));
    Na_eff = Na_o.*exp(-0.820000.*VD);
    P_E2Snk = 1.00000./(1.00000+(K_nknao./Na_eff).*(1.00000+K_o./K_nkko));
    I_nk = k_nk.*(k_1nk.*P_E1Snk.*y_nk-k_2nk.*P_E2Snk.*(1.00000-y_nk)).*1.00000;
    J_Na = (-1.00000./(F.*vol_cyt)).*(3.00000.*I_nk+3.00000.*I_xm+I_Na+I_Nalk+I_NaHCN);
    
    P_E1Dnk = 1.00000./(1.00000+(K_nkki./K_i).*(1.00000+Na_i./K_nknai));
    alpha_nk = k_1nk.*P_E1Snk+k_3nk.*P_E1Dnk;
    P_E2Dnk = 1.00000./(1.00000+(K_nkko./K_o).*(1.00000+Na_eff./K_nknao));
    beta_nk = k_2nk.*P_E2Snk+k_4nk.*P_E2Dnk;
    
    V_K = log(K_o./K_i);
    O_sk = (power(Ca_i, 4.20000))./(power(0.000350000, 4.20000)+power(Ca_i, 4.20000));
    I_Ksk = (g_ksk.*O_sk.*(power(K_i.*K_o, 1.0./2)).*sinh(0.500000.*(VD-V_K)))./(sinh(0.500000.*VD)./(0.500000.*VD));
    O_kdr = power(m_kdr, 3.00000);
    I_Kdr = g_kdr.*O_kdr.*(V_snc-V_K.*V_tau);
    O_kir = 1.00000./(1.00000+exp((V_snc+85.0000)./12.1000));
    I_Kir = g_kir.*O_kir.*(V_snc-V_K.*V_tau);
    I_K = I_Ksk+I_Kdr+I_Kir;
    J_K = (-1.00000./(F.*vol_cyt)).*(I_K-2.00000.*I_nk);
    
    % ER
    J_pump = k_pump.*Ca_i.*ATP; %J_pump
    J_ch = k_ch.*((power(Ca_i, 2.00000))./(power(K1, 2.00000)+power(Ca_i, 2.00000))).*(Ca_er-Ca_i); %J_ch
    J_leak = k_leak.*(Ca_er-Ca_i); %J_leak
    
    % Mito
    J_out = (k_out.*((power(Ca_i, 2.00000))./(power(K3, 2.00000)+power(Ca_i, 2.00000)))+k_m).*Ca_mt; % J_out
    J_in = k_in.*((power(Ca_i, 8.00000))./(power(K2, 8.00000)+power(Ca_i, 8.00000)));%.*ATP; % J_in
    
    % Calcium dynamics
    J_Ca = J_ca-(J_calb+4.00000.*J_cam)-J_pump+J_ch+J_leak-J_in+J_out;
    
    adca=0;
    %Terminal
    
    % ATP-dependent DA packing
    %     ada=RescaleRange(ATP,0.2,2.3,0.001,1);
    ada=0.001.*(exp(3.*ATP));
    
    % ATP-dependent vescile recycling
    %     nRRP=5;
    nRRP=1.*(exp(1.*ATP));
    
    Vsynt = Vsynt_max./(((Ksynt./(adca+Ca_i))^4)+1);
    jsynt = (Vsynt./(1+((Ktyr./TYR).*(1+(cda./Kicda)+(eda./Kieda)))));
    
    jvmat = ada.*Vcda_max .* MM_kin(cda,Kcda,1);%.*ATP;
    
    jida = kmao .* cda;
    
    %     nRRP = (40./((1+exp(-(vda-vdao)./vdas)).*(1+exp((eda-dara)./dars))));
    prob = 0.14 .* MM_kin((adca+Ca_i),krel,4);
    jrel = psi .* nRRP .* prob;
    
    jdat = Veda_max .* MM_kin(eda,Keda,1);
    
    jeda = kcomt .* eda;
    
    % Energy metabolism
    % Energy consumed by active pumps
    % V_pumps=0;
    V_pumps1 = 1.*(1.00000./(F.*vol_cyt)).*(I_nk+I_pmca);
    V_pumps2 = 1.*(jvmat);
    V_pumps3 = 100.*jrel;
    
    v_stim=0;
    % v_stim1=(1.00000/(F*vol_cyt))*(I_nk+I_pmca);
    v_stim1=0.0.*(I_nk+I_pmca);
    J_er=(beta_er./rho_er).*(J_pump);
    
    %     V_id(k)=V_pumps1;
    %     V_dp(k)=V_pumps2;
    %     V_er(k)=J_er;
    
    V_pumps=V_pumps1+V_pumps2+J_er+V_pumps3;
    
    uADP = power(Q_adk, 2.00000)+4.00000.*Q_adk.*(ANP./ATP-1.00000);
    ADP = (ATP./2.00000).*(-Q_adk+power(uADP, 1.0./2));
    Cr = PCr_tot-PCr;
    V_ck = 0.*(kf_ck.*PCr.*ADP-kr_ck.*Cr.*ATP);
    
    ATP_inh = power((1.00000+nATP.*(ATP./KI_ATP))./(1.00000+ATP./KI_ATP), 4.00000);
    V_pk = enfatra.*Vmax_pk.*(GAP./(GAP+Km_GAP_pk)).*(ADP./(ADP+Km_ADP_pk)).*ATP_inh;
    pa=1;
    V_op = enfamit.*Vmax_op.*((pa.*PYR)./((pa.*PYR)+Km_PYR_op)).*(ADP./(ADP+Km_ADP_op)).*(1.00000./(1.00000+0.100000.*(ATP./ADP)));
    
    AMP = ANP-(ATP+ADP);
    AMP_act = power((1.00000+AMP./Ka_AMP)./(1.00000+nAMP.*(AMP./Ka_AMP)), 4.00000);
    V_pfk = Vmax_pfk.*(F6P./(F6P+Km_F6P_pfk)).*(ATP./(ATP+Km_ATP_pfk)).*(F26P./(F26P+Km_F26P_pfk)).*ATP_inh.*AMP_act;
    
    eta_ldh=1-beta_ldh_ros.*((ROS^4)./((ROS^4)+(Kldh_ros^4)));
    V_ldh = 1.*eta_ldh.*(kf_ldh.*PYR-kr_ldh.*LAC);
    V_lac = Vlac_0.*(1.00000+v_stim1.*K_lac)-K_lac_eff.*LAC;
    
    V_hk = Vmax_hk.*(ATP./(ATP+Km_ATP_hk)).*(power(1.00000+power(F6P./KI_F6P, 4.00000), -1.00000)).*GLCe;
    AMP_pfk2 = (power(AMP./Kamp_pfk2, nh_amp))./(1.00000+power(AMP./Kamp_pfk2, nh_amp));
    V_pfk2 = Vmaxf_pfk2.*(ATP./(ATP+Km_ATP_pfk2)).*(F6P./(F6P+Km_F6P_pfk2)).*AMP_pfk2-Vmaxr_pfk2.*(F26P./(F26P+Km_F26P_pfk2));
    
    V_ATPase = Vmax_ATPase.*(ATP./(ATP+Km_ATP)).*(1.00000+v_stim);
    dAMP_dATP = -1.00000+Q_adk./2.00000+-(0.500000.*(power(uADP, 1.0./2)))+Q_adk.*(ANP./(ATP.*(power(uADP, 1.0./2))));
    
    eta_op=eta_op_max-beta_op_asyn.*(((ASYNA^4)./((ASYNA^4)+(Kasyn^4))));
    %%
    %%%%%%%%%%%%%%%%%%%%%% Differential equations %%%%%%%%%%%%%%%%%%%%%%%%%
    
    V_sncnxt = V_snc + (((F.*vol_cyt)./(C_sp.*A_pmu)).*(J_Na+J_K+2.00000.*J_Ca+Iext-Igabasnc+(scfa*Istnsnc))).*dt;
    Ca_inxt = Ca_i + (J_Ca).*dt;
    Na_inxt = Na_i + (J_Na).*dt;
    K_inxt = K_i + (J_K).*dt;
    Calbnxt = Calb + (-J_calb).*dt;
    Camnxt = Cam + (-J_cam).*dt;
    m_calnxt = m_cal + ((1.00000./(1.00000+exp(-(V_snc+15.0000)./7.00000))-m_cal)./(7.68000.*exp(-(power((V_snc+65.0000)./17.3300, 2.00000)))+0.723100)).*dt;
    m_nanxt = m_na + (A_mna.*exp(za_mna.*VD).*(1.00000-m_na)-B_mna.*exp(-zb_mna.*VD).*m_na).*dt;
    h_nanxt = h_na + (A_hna.*exp(za_hna.*VD).*(1.00000-h_na)-B_hna.*exp(-zb_hna.*VD).*h_na).*dt;
    O_hcnnxt = O_hcn + (kf_hcn.*(1.00000-O_hcn)-kr_hcn.*O_hcn).*dt;
    m_kdrnxt = m_kdr + ((1.00000./(1.00000+exp(-(V_snc+25.0000)./12.0000))-m_kdr)./(18.0000./(1.00000+exp((V_snc+39.0000)./8.00000))+1.00000)).*dt;
    y_pcnxt = y_pc + (beta_pc.*(1.00000-y_pc)-alpha_pc.*y_pc).*dt;
    y_nknxt = y_nk + (beta_nk.*(1.00000-y_nk)-alpha_nk.*y_nk).*dt;
    ATPusednxt = ATPused+(-ATPused+(1.00000./(F.*vol_cyt)).*(I_nk+I_pmca)).*dt;%ATPused
    Ca_ernxt = Ca_er + ((beta_er./rho_er).*(J_pump-(J_ch+J_leak))).*dt;
    Ca_mtnxt = Ca_mt + ((beta_mt/rho_mt).*(J_in-J_out)).*dt;
    mCatrajectorySNC(k, 1) = Ca_mtnxt(randi([1 8],1),randi([1 8],1));
    mCatrajectorySNC(k, 2) = Ca_mtnxt(randi([1 8],1),randi([1 8],1));
    mCatrajectorySNC(k, 3) = Ca_mtnxt(randi([1 8],1),randi([1 8],1));
    mCatrajectorySNC(k, 4) = Ca_mtnxt(randi([1 8],1),randi([1 8],1));
    mCatrajectorySNC(k, 5) = Ca_mtnxt(randi([1 8],1),randi([1 8],1));
    cdanxt = cda+(jsynt + jdat - jvmat - jida).*dt;%cda
    vdanxt = vda+(jvmat - jrel).*dt;%vda
    edanxt = eda+(jrel - jdat - jeda).*dt;%eda
    calnxt = cal+(-k3f.*(Sig_ers.*cal)+k3b.*(cai_cal)).*dt;%cal
    cai_calnxt = cai_cal+(k3f.*(Sig_ers.*cal)-k3b.*(cai_cal)-k4f.*(cai_cal)).*dt;%cai_cal
    cal_actnxt = cal_act+(k4f.*(cai_cal)-k5f.*(cal_act.*casp12)+k5b.*(cal_act_casp12)).*dt;%cal_act
    casp12nxt = casp12+(-k5f.*(cal_act.*casp12)+k5b.*(cal_act_casp12)).*dt;%casp12
    cal_act_casp12nxt = cal_act_casp12+(k5f.*(cal_act.*casp12)-k5b.*(cal_act_casp12)-k6f.*(cal_act_casp12)).*dt;%cal_act_casp12
    casp12_actnxt = casp12_act+(k6f.*(cal_act_casp12)-k7f.*(casp12_act.*casp9)+k7b.*(casp12_act_casp9)).*dt;%casp12_act
    casp9nxt = casp9+(-k7f.*(casp12_act.*casp9)+k7b.*(casp12_act_casp9)).*dt;%casp9
    casp12_act_casp9nxt = casp12_act_casp9+(k7f.*(casp12_act.*casp9)-k7b.*(casp12_act_casp9)-k8f.*(casp12_act_casp9)).*dt;%casp12_act_casp9
    casp9_actnxt = casp9_act+(k8f.*(1.*casp12_act_casp9)+k9b.*(casp9_act_casp3)-k9f.*(casp9_act.*casp3)+1.*k28f.*Cytc_casp9-k12f.*casp9_act.*IAP+k12b.*casp9_act_IAP).*dt;%casp9_act
    casp3nxt = casp3+(-k9f.*(casp9_act.*casp3)+k9b.*(casp9_act_casp3)).*dt;%casp3
    casp9_act_casp3nxt = casp9_act_casp3+(-k10f.*(casp9_act_casp3)-k9b.*(casp9_act_casp3)+k9f.*(casp9_act.*casp3)).*dt;%casp9_act_casp3
    casp3_actnxt = casp3_act+(k10f.*(casp9_act_casp3)-k11f.*(casp9_act.*casp3_act)-k13f.*casp3_act.*IAP+k13b.*casp3_act_IAP).*dt;%casp3_act
    apopnxt =  apop+(asb.*k11f.*(casp9_act.*casp3_act)).*dt;%apop
    
    F6Pnxt = F6P+(V_hk-(V_pfk-V_pfk2)).*dt;%F6P
    F26Pnxt = F26P+(V_pfk2).*dt;%F26P
    GAPnxt = GAP+(2.00000.*V_pfk-V_pk).*dt;%GAP
    PYRnxt = PYR+(V_pk-(V_op+V_ldh)).*dt;%PYR
    LACnxt = LAC+(2.25000.*V_ldh+V_lac).*dt;%LAC
    ATPnxt = ATP+(((1.*(1.*2.00000.*V_pk+15.0000.*eta_op.*V_op+V_ck))-(V_hk+V_pfk+V_pfk2+V_ATPase+V_pumps)).*(power(1.00000-dAMP_dATP, -1.00000))).*dt;%ATP
    PCrnxt = PCr+(-V_ck).*dt;%PCr
    
    ROS_mitnxt = ROS_mit+(k29f.*Sig_mts.*Mit).*dt;%ROS_mit
    PTP_mit_actnxt = PTP_mit_act+(k30f.*ROS_mit.*PTP_mit).*dt;%PTP_mit_act
    Cytc_mitnxt = Cytc_mit+(-k31f.*PTP_mit_act.*Cytc_mit).*dt;%Cytc_mit
    Cytcnxt = Cytc+(-k27f.*Cytc.*casp9+k27b.*Cytc_casp9+k31f.*PTP_mit_act.*Cytc_mit).*dt;%Cytc
    Cytc_casp9nxt = Cytc_casp9+(k27f.*Cytc.*casp9-k27b.*Cytc_casp9-k28f.*Cytc_casp9).*dt;%Cytc_casp9
    IAPnxt = IAP+(-k12f.*casp9_act.*IAP+k12b.*casp9_act_IAP-k13f.*casp3_act.*IAP+k13b.*casp3_act_IAP).*dt;%IAP
    casp9_act_IAPnxt = casp9_act_IAP+(k12f.*casp9_act.*IAP-k12b.*casp9_act_IAP).*dt;%casp9_act_IAP
    casp3_act_IAPnxt = casp3_act_IAP+(k13f.*casp3_act.*IAP-k13b.*casp3_act_IAP).*dt;%casp3_act_IAP
    
    
    %     inds=find(V <=80 & V >=0);
    %     snc_firings=[snc_firings; k+0*inds,inds+0*inds];
    
    %     if k>sncstart/dt
    %         [snc_firings]=ConvertAPtoST(snc_firings,Psnc);
    %         sncstart=sncstart+500;
    %     end
    
    V_snc=V_sncnxt;m_cal=m_calnxt;m_kdr=m_kdrnxt;m_na=m_nanxt;
    VtrajectorySNC(k,:) = V_snc(1,:);
    h_na=h_nanxt;O_hcn=O_hcnnxt;Calb=Calbnxt;
    Cam=Camnxt;y_nk=y_nknxt;y_pc=y_pcnxt;
    K_i=K_inxt;Na_i=Na_inxt;Ca_i=Ca_inxt;
    Ca_trajectorySNC(k,:) = Ca_i(1,:); % keep track of the internal calcium concentration

    Ca_er=Ca_ernxt;Ca_mt=Ca_mtnxt;
    er_catrajectorySNC(k,:) = Ca_er(1,:);
    mt_catrajectorySNC(k,:) = Ca_mt(1,:);
%     erCatrajectorySNC(k, 1) = Ca_ernxt(randi([1 8],1),randi([1 8],1));
%     erCatrajectorySNC(k, 2) = Ca_ernxt(randi([1 8],1),randi([1 8],1));
%     erCatrajectorySNC(k, 3) = Ca_ernxt(randi([1 8],1),randi([1 8],1));
%     erCatrajectorySNC(k, 4) = Ca_ernxt(randi([1 8],1),randi([1 8],1));
%     erCatrajectorySNC(k, 5) = Ca_ernxt(randi([1 8],1),randi([1 8],1));
    cda=cdanxt;vda=vdanxt;eda=edanxt;ATPused=ATPusednxt;
    cal=calnxt;cai_cal=cai_calnxt;cal_act=cal_actnxt;
    casp12=casp12nxt;cal_act_casp12=cal_act_casp12nxt;casp12_act=casp12_actnxt;
    casp9=casp9nxt;casp12_act_casp9=casp12_act_casp9nxt;casp9_act=casp9_actnxt;
    casp3=casp3nxt;casp9_act_casp3=casp9_act_casp3nxt;casp3_act=casp3_actnxt;
    apop=apopnxt;
    ROS_mit=ROS_mitnxt;
    PTP_mit_act=PTP_mit_actnxt;
    Cytc_mit=Cytc_mitnxt;
    Cytc=Cytcnxt;
    Cytc_casp9=Cytc_casp9nxt;
    IAP=IAPnxt;
    casp9_act_IAP=casp9_act_IAPnxt;
    casp3_act_IAP=casp3_act_IAPnxt;
    F6P=F6Pnxt;
    F26P=F26Pnxt;
    GAP=GAPnxt;
    PYR=PYRnxt;
    LAC=LACnxt;
    ATP=ATPnxt;
    PCr=PCrnxt;
    Ssnc=Ssncnxt;
    
    %     inds=find(V_snc_array(k) >= -20 && V_snc_array(k) > V_snc_array(k-1) && V_snc_array(k) > V_snc_array(k+1));
    %     inds=find(V_snc <=80 & V_snc >=-20);
    %     snc_firings=[snc_firings; k+0*inds,inds+0*inds];
    %%
    %----------------------------------------STN-----------------------------------------%
    %---------------------------Input from stn to stn(laterals)--------------------------%
    
    % psp variable
    h_nmdastn = (1-lam_nmda).* h_nmdastn + lam_nmda.*stncurr_spk; %psp nmda stn lat
    h_ampastn = (1-lam_ampa).* h_ampastn + lam_ampa.*stncurr_spk; %psp ampa stn lat
    
    tmplat_nmda_stn = h_nmdastn.*(Enmda - Vstn);
    tmplat_ampa_stn = h_ampastn.*(Eampa - Vstn);
    
    Ilat_nmda_stn = conv2(tmplat_nmda_stn, wlatstn, 'same'); % lat nmda stn
    Ilat_ampa_stn = conv2(tmplat_ampa_stn, wlatstn, 'same');% lat ampa stn
    
    B = 1./(1 + (mg0/3.57).*exp(-0.062.*Vstn));
    
    %----------------------------------------Input from gpe to stn---------------------%
    % psp variable
    h_gs = (1-lam_gaba).* h_gs + lam_gaba.*gpecurr_spk;% input from gpe to nmda stn
    % gaba current
    I_gs = wgs.*h_gs.*(Egaba - Vstn);
    
    Istnstn=B.*Ilat_nmda_stn + Ilat_ampa_stn;
    % total current stn recieves
    Itmpstn =  Istnstn + I_gs; % total currents
    
    % V,U updated
    dvstn = taustn.*(((0.04.*Vstn.*Vstn)+5.*Vstn - Ustn +140 + Istn + Itmpstn)./Cstn);%+wgn(Mstn,Nstn,lnoise_st,imp_st);
    dustn = taustn.*astn.*(bstn.*(Vstn) - Ustn);
    Vstn_nxt = Vstn + dvstn;
    Ustn_nxt = Ustn + dustn;
    
    inds = find(Vstn_nxt > vpeak_stn);
    %     stn_firings=[stn_firings; k+0*inds,inds];
    
    Vstn_nxt(inds) = cstn.*ones(size(inds));
    Ustn_nxt(inds) = Ustn(inds) + dstn.*ones(size(inds));
    stncurr_spk = stn_zeros;
    stncurr_spk(inds) = ones(size(inds));
    
    Vstn = Vstn_nxt;
    Ustn = Ustn_nxt;
    %for adding stress, it might need to be cumulative.
    
    STNcells2track = [2, 4; 1, 8; 4, 3;];
    for ncell = 1:size(STNcells2track,1)      
        VtrajectorySTN(k,ncell,1) = Vstn(STNcells2track(ncell,1),STNcells2track(ncell,2));
    end


    
    %%
    %----------------------------------------GPe-----------------------------------------%
    %------------- -------- Input from stn to gpe----------------------------------%
    % psp variable
    h_nmdagpe = (1-lam_nmda).* h_nmdagpe + lam_nmda.*stncurr_spk;
    h_ampagpe = (1-lam_ampa).* h_ampagpe + lam_ampa.*stncurr_spk;
    
    % nmda and ampa currents
    Inmdagpe = wsg.*h_nmdagpe.*(Enmda - Vgpe);
    Iampagpe = wsg.*h_ampagpe.*(Eampa - Vgpe);
    
    %-------------------------- Input from gpe to gpe(laterals)-------------------%
    % psp variable
    hlat_gaba_gpe = (1-lam_gaba).* hlat_gaba_gpe + lam_gaba.*gpecurr_spk;
    tmplat_gaba_gpe = wda_gpe.*hlat_gaba_gpe.*(Egaba - Vgpe);
    
    % lateral currents
    Ilatgpe = conv2(tmplat_gaba_gpe, wlatgpe, 'same');
    
    B = 1./(1 + (mg0/3.57).*exp(-0.062.*Vgpe));
    I_sg=B.*Inmdagpe+Iampagpe;
    Itmpgpe = Ilatgpe + I_sg;
    
    % LFP
    % LFP_GPe_exc(k)=sum(sum(I_sg));
    % LFP_GPe_inh(k)=sum(sum(Ilatgpe));
    % LFP_GPe_tot(k)=LFP_GPe_exc(k)+LFP_GPe_inh(k);
    
    % V ,U updated
    dvgpe = taugpe.*((0.04.*Vgpe.*Vgpe)+5.*Vgpe+140 - Ugpe + Itmpgpe+Igpe)./Cgpe;
    dugpe = taugpe.*agpe.*(bgpe.*(Vgpe) - Ugpe);
    Vgpe_nxt = Vgpe + dvgpe;
    Ugpe_nxt = Ugpe + dugpe;
    
    inds = find(Vgpe_nxt > vpeak_gpe);
    %     gpe_firings=[gpe_firings; k+0*inds,inds];
    
    Vgpe_nxt(inds) = cgpe.*ones(size(inds));
    Ugpe_nxt(inds) = Ugpe(inds) + dgpe.*ones(size(inds));
    gpecurr_spk = gpe_zeros;
    gpecurr_spk(inds) = ones(size(inds));
    
    Vgpe = Vgpe_nxt;
    Ugpe = Ugpe_nxt;
    VtrajectoryGPE(k,:) = Vgpe(1,:);

    % Killing based on apoptosis threshold
    indsap=find(apop>apopthr);
    indsapp=[indsapp indsap'];
    if isempty(indsap)==0
        indsapp=unique(indsapp);
        indsap=[];
    end
    V_snc(indsapp) = -80.*ones(size(indsapp));
    apop(indsapp)=zeros(size(indsapp));
    Ca_i(indsapp)=zeros(size(indsapp));
    cda(indsapp) = 1e-4.*zeros(size(indsapp)); %mM
    vda(indsapp) = 500.*zeros(size(indsapp)); %mM
    eda(indsapp) = 26e-6.*zeros(size(indsapp)); %mM
    ATPused(indsapp)=0.*zeros(size(indsapp));
    ATP(indsapp)=0.*zeros(size(indsapp));
    
    % Storing pattern of SNc cell loss
    %     indsnum = find(V_snc < -70);
    %     if isempty(indsnum)==0
    %         indsnums=setdiff(indsnum,indsappp);
    %         if isempty(indsnums)==0
    %             indsappp=[indsappp; k+0*indsnums,indsnums];
    %         end
    %     end
    %
    
    %     clit{k}=setdiff(indsapp,indsappp);
    
    %     indsappp=indsapp;
    
    eda1=sum(sum(eda))/(Psnc);
    
    yeda1(count)=sum(sum(eda))/(Psnc);
    dDA1(count)=DA;
    dnc1(count)=numel(indsapp);
    
    Inc=(numel(indsapp));
    
    count=count+1;
    % Sub-sampling iterative
    if (k==samp_start)
        
        if gpuon==1
            yeda0=sub_sampling_GPU(yeda1,dt);
            dnc0=sub_sampling_GPU(dnc1,dt);
            dDA0=sub_sampling_GPU(dDA1,dt);
        elseif gpuon==0
            yeda0=sub_sampling(yeda1,dt);
            dnc0=sub_sampling(dnc1,dt);
            dDA0=sub_sampling(dDA1,dt);
        end
        
        deda=[deda yeda0];
        dnc=[dnc dnc0];
        dDA=[dDA dDA0];
        
        samp_start=samp_start+(stt/dt);
        count=1;
        disp(['Current time is ',num2str(k*dt),' ms'])
    end
    
    
    
end

kid=((Psnc)-dnc);
% out_cell={};
% out_cell{1,1}=clit;
% out_cell{2,1}=kid;
% out_cell{3,1}=deda;
% out_cell{4,1}=dDA;

stt=toc;
if stt>60 && stt<=3600
    stt1=stt/60;
    f1=['Simulation time is ',num2str(stt1),' minutes_',num2str(asb)];
    %     disp(f1);
elseif stt>3600
    stt2=stt/(60*60);
    f1=['Simulation time is ',num2str(stt2),' hours_',num2str(asb)];
    %     disp(f1);
else
    f1=['Simulation time is ',num2str(stt),' seconds_',num2str(asb)];
    %     disp(f1);
end
simtime=f1;
