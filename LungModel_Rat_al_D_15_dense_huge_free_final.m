% Adapted from: Yu JY, Rosania GR. Cell-based multiscale computational modeling of small molecule absorption and retention in the lungs. Pharm Res. 2010;27(3):457–467.
% This implementation is an independent adaptation and extension of the published framework, rather than the original code released by the authors.
% Major methodological extensions and improvements:
% (1) Stabilized Nernst–Planck factor
% Original form: phi(N) = N / (exp(N) - 1)
% Modified with a numerically stable implementation:
% phi(N) ≈ 1 - N/2 + N^2/12 for |N| < 1e-6
% phi(N) = N / (exp(N) - 1) otherwise
% This avoids numerical instability when N is close to zero.
% (2) Explicit handling of neutral species (z = 0)
% Neutral molecules bypass ionization partitioning in the implementation.
% The neutral-only branch sets ionized fractions to zero and avoids non-physical ionization-related transport artifacts.
% (3) Extension to a high-throughput BB-PTK platform
% This implementation extends the original transport framework into a burden-based pulmonary tissue kinetics (BB-PTK) platform by:
% converting concentration-based outputs to compartmental amounts using compartment volumes
% enabling automated batch processing across multiple compounds
% supporting automated extraction of dosimetry-oriented outputs such as lung burden (mol), plasma burden (mol), Tmax, and AUP
% Lung burden definition: Z_lung = Z(1)+Z(2)+Z(3)+Z(4)+Z(5)+Z(6)+Z(7), where Z denotes compartmental amount after converting concentration-based state variables using compartment volumes.
% Author: Wanjun Zhang
% Year: 2026
% -------------------------------------------------------------------------
% 8 compartment aEp, imEp, cEp, int, imInt, sm, cEd, p
% Alveolar region
% Given dose, logP and pKa, charge, calculate the absorption profile
% all geometric parameters units are based on meter m m2 m3
% membrane potential: V
% time : seconds
% amount : mol concentration: mol/m3 or mM

function [t_half_a,t_half_d,t,Y,M_v,M_v_free_n,M_v_free_d,M_v_free_t, ...
          PLung_h,PSL_h,PMlung_max,TLung_max,PMimEp_max,PMcEp_max, ...
          PMint_max,PMsm_max,PMimInt_max,PMcEd_max] = ...
          LungModel_Rat_al_D_15_dense_huge_free_final(logPN,pKa,z,dose)
 
% ===================== NUMERICAL STABILIZATION PATCHES =====================
% 1) Avoid ode15s failure at t=0 caused by extreme tolerances + ultra-dense output grid:
%    - Use tspan = [0 604800] and let ode15s choose its internal steps.
%    - Use realistic strict tolerances: RelTol=1e-8, AbsTol=1e-10.
% 2) Use expm1(N) for N/(exp(N)-1) to prevent catastrophic cancellation at small N.
% 3) For z==0, bypass HH ionization split entirely (neutral-only) to avoid unstable intermediates.
% 4) Sanity check M/G for NaN/Inf before ODE solve.
% 5) Whole-lung burden used for model¨Cdata comparison now includes aEp.
% ===========================================================================
 
%% ------------------ Constants and parameters ------------------
T = 273.15+37; 
R = 8.314; 
F = 96485.34;
 
LaEp = 0.95; LimEp = 0; LcEp = 0.05; Lint = 0.05; LimInt = 0; Lsm = 0; LcEd = 0.05; Lp = 0;
WaEp = 1 - LaEp; WimEp = 1 - LimEp; WcEp = 1 - LcEp; Wint = 1 - Lint; WimInt = 1 - LimInt;
Wsm  = 1 - Lsm;  WcEd = 1 - LcEd;   Wp   = 1 - Lp;
 
GaEpN = 1; GaEpD = 1;
GimEpN = 1.23; GimEpD = 0.74;
GcEpN = 1.23; GcEpD = 0.74;
GintN = 1; GintD = 1;
GimIntN = 1.23; GimIntD = 0.74;
GsmN = 1.23; GsmD = 0.74;
GcEdN = 1.23; GcEdD = 0.74;
GpN = 1; GpD = 1;
 
%% ------------------ Areas (m2) ------------------
AaEp = 0.387; 
AbEp = AaEp;
AimEp = 3.14e-10 * 0.89e9 * 3/100/2;
AimInt = AimEp/10; 
Asm = 0; 
AbEd = 0.452; 
AaEd = 0.452;
 
%% ------------------ Volumes (m3) ------------------
ASL = 5; 
VaEp = AaEp * ASL * 1e-6;
VcEp = AaEp * 0.384e-6;
VimEp = 0.89e9 * 3/100 * 1058e-18;
Vint = AaEp * 0.693e-6;
VimInt = VimEp / 20;
Vsm = VcEp * 1e-12;
VcEd = AbEd * 0.358e-6;
Vp = 5;
 
M_v = diag([VaEp,VimEp,VcEp,Vint,Vsm,VimInt,VcEd,Vp]);
 
%% ------------------ Membrane potentials (V) ------------------
EaEp = -0.0093; 
EbEp = 0.0119; 
EimEp = -0.06; 
EimInt = -0.06;
Esm  = -0.06;   
EbEd = -0.06;  
EaEd  = -0.03;
 
%% ------------------ pH values ------------------
pHaEp = 7.4; 
pHimEp = 7.0; 
pHcEp = 7.0; 
pHint = 7.0; 
pHimInt = 7.0;
pHsm  = 7.0; 
pHcEd  = 7.0; 
pHp   = 7.4;
 
%% ------------------ Charge handling & stable Nernst¨CPlanck helper ------------------
isNeutralMol = isfinite(z) && abs(z) < 1e-12;
 
% Stable phi(N) = N/(exp(N)-1), with correct limit phi(0)=1
phiNP = @(N) phiNP_safe(N);
 
%% ------------------ Physicochemical properties ------------------
logPD = logPN - 3.7;
 
if abs(z - 1) <= 1e-6
    logP_nlipT = 0.33*logPN + 2.2;
    logP_dlipT = 0.37*logPD + 2;
elseif abs(z + 1) <= 1e-6
    logP_nlipT = 0.37*logPN + 2.2;
    logP_dlipT = 0.33*logPD + 2.6;
else
    logP_nlipT = 0.33*logPN + 2.2;
    logP_dlipT = 0.33*logPD + 2.2;
end
 
logP_n = round(logP_nlipT*100)/100;
logP_d = round(logP_dlipT*100)/100;
 
Pn = 10^(logP_n - 6.7);
Pd = 10^(logP_d - 6.7);
 
N = 1.22 * 10^(logP_n);
D = 1.22 * 10^(logP_d);
 
%% ------------------ Partition coefficients ------------------
KaEpN = N * LaEp; KaEpD = D * LaEp;
KimEpN = N * LimEp; KimEpD = D * LimEp;
KcEpN = N * LcEp; KcEpD = D * LcEp;
KintN = N * Lint; KintD = D * Lint;
KimIntN = N * LimInt; KimIntD = D * LimInt;
KsmN = N * Lsm; KsmD = D * Lsm;
KcEdN = N * LcEd; KcEdD = D * LcEd;
KpN = N * Lp; KpD = D * Lp;
 
%% ------------------ Flux potentials ------------------
NaEp   = z * EaEp * F / (R*T);
NbEp   = z * -EbEp * F / (R*T);
NimEp  = z * EimEp * F / (R*T);
NimInt = z * EimInt * F / (R*T);
Nsm    = z * Esm * F / (R*T);
NbEd   = z * EbEd * F / (R*T);
NaEd   = z * -EaEd * F / (R*T);
 
%% ------------------ Distribution fractions (fn/fd for 8 compartments) ------------------
if isNeutralMol
    faEpN   = 1/(WaEp/GaEpN + KaEpN/GaEpN);                 faEpD   = 0;
    fimEpN  = 1/(WimEp/GimEpN + KimEpN/GimEpN);             fimEpD  = 0;
    fcEpN   = 1/(WcEp/GcEpN + KcEpN/GcEpN);                 fcEpD   = 0;
    fintN   = 1/(Wint/GintN + KintN/GintN);                 fintD   = 0;
    fimIntN = 1/(WimInt/GimIntN + KimIntN/GimIntN);         fimIntD = 0;
    fsmN    = 1/(Wsm/GsmN + KsmN/GsmN);                     fsmD    = 0;
    fcEdN   = 1/(WcEd/GcEdN + KcEdN/GcEdN);                 fcEdD   = 0;
    fpN     = 1/(Wp/GpN + KpN/GpN);                         fpD     = 0;
else
    ichg = -sign(z);
 
    faEpN = 1/(WaEp/GaEpN + KaEpN/GaEpN + WaEp*10^(ichg*(pHaEp-pKa))/GaEpD + KaEpD*10^(ichg*(pHaEp-pKa))/GaEpD);
    faEpD = faEpN * 10^(ichg*(pHaEp-pKa));
 
    fimEpN = 1/(WimEp/GimEpN + KimEpN/GimEpN + WimEp*10^(ichg*(pHimEp-pKa))/GimEpD + KimEpD*10^(ichg*(pHimEp-pKa))/GimEpD);
    fimEpD = fimEpN * 10^(ichg*(pHimEp-pKa));
 
    fcEpN = 1/(WcEp/GcEpN + KcEpN/GcEpN + WcEp*10^(ichg*(pHcEp-pKa))/GcEpD + KcEpD*10^(ichg*(pHcEp-pKa))/GcEpD);
    fcEpD = fcEpN * 10^(ichg*(pHcEp-pKa));
 
    fintN = 1/(Wint/GintN + KintN/GintN + Wint*10^(ichg*(pHint-pKa))/GintD + KintD*10^(ichg*(pHint-pKa))/GintD);
    fintD = fintN * 10^(ichg*(pHint-pKa));
 
    fimIntN = 1/(WimInt/GimIntN + KimIntN/GimIntN + WimInt*10^(ichg*(pHimInt-pKa))/GimIntD + KimIntD*10^(ichg*(pHimInt-pKa))/GimIntD);
    fimIntD = fimIntN * 10^(ichg*(pHimInt-pKa));
 
    fsmN = 1/(Wsm/GsmN + KsmN/GsmN + Wsm*10^(ichg*(pHsm-pKa))/GsmD + KsmD*10^(ichg*(pHsm-pKa))/GsmD);
    fsmD = fsmN * 10^(ichg*(pHsm-pKa));
 
    fcEdN = 1/(WcEd/GcEdN + KcEdN/GcEdN + WcEd*10^(ichg*(pHcEd-pKa))/GcEdD + KcEdD*10^(ichg*(pHcEd-pKa))/GcEdD);
    fcEdD = fcEdN * 10^(ichg*(pHcEd-pKa));
 
    fpN = 1/(Wp/GpN + KpN/GpN + Wp*10^(ichg*(pHp-pKa))/GpD + KpD*10^(ichg*(pHp-pKa))/GpD);
    fpD = fpN * 10^(ichg*(pHp-pKa));
end
 
%% ------------------ Free-volume matrices ------------------
M_v_free_n = diag([VaEp*faEpN, VimEp*fimEpN, VcEp*fcEpN, Vint*fintN, Vsm*fsmN, VimInt*fimIntN, VcEd*fcEdN, Vp*fpN]);
M_v_free_d = diag([VaEp*faEpD, VimEp*fimEpD, VcEp*fcEpD, Vint*fintD, Vsm*fsmD, VimInt*fimIntD, VcEd*fcEdD, Vp*fpD]);
M_v_free_t = M_v_free_n + M_v_free_d;
 
%% ------------------ Coefficient matrix for ODEs ------------------
KaEp_aEp = -AaEp/VaEp*(Pn*faEpN + Pd*phiNP(NaEp)*faEpD) ...
           -AimEp/VaEp*(Pn*faEpN + Pd*phiNP(NimEp)*faEpD);
KaEp_imEp = -AimEp/VaEp*(Pn*(-fimEpN) + Pd*phiNP(NimEp)*(-fimEpD)*exp(NimEp));
KaEp_cEp  = -AaEp/VaEp*(Pn*(-fcEpN)  + Pd*phiNP(NaEp)*(-fcEpD)*exp(NaEp));
KaEp_int  = 0; KaEp_sm = 0; KaEp_imInt = 0; KaEp_cEd = 0; KaEp_p = 0;
SaEp = 0;
 
KimEp_aEp  = AimEp/VimEp*(Pn*faEpN + Pd*phiNP(NimEp)*faEpD);
KimEp_imEp = AimEp/VimEp*(Pn*(-fimEpN) + Pd*phiNP(NimEp)*(-fimEpD)*exp(NimEp));
KimEp_cEp  = 0; KimEp_int = 0; KimEp_sm = 0; KimEp_imInt = 0; KimEp_cEd = 0; KimEp_p = 0;
SimEp = 0;
 
KcEp_aEp = AaEp/VcEp*(Pn*faEpN + Pd*phiNP(NaEp)*faEpD);
KcEp_imEp = 0;
KcEp_cEp = AaEp/VcEp*(Pn*(-fcEpN) + Pd*phiNP(NaEp)*(-fcEpD)*exp(NaEp)) ...
          -AbEp/VcEp*(Pn*fcEpN + Pd*phiNP(NbEp)*fcEpD);
KcEp_int = -AbEp/VcEp*(Pn*(-fintN) + Pd*phiNP(NbEp)*(-fintD)*exp(NbEp));
KcEp_sm = 0; KcEp_imInt = 0; KcEp_cEd = 0; KcEp_p = 0;
ScEp = 0;
 
Kint_aEp = 0; Kint_imEp = 0;
Kint_cEp = AbEp/Vint*(Pn*fcEpN + Pd*phiNP(NbEp)*fcEpD);
Kint_int = AbEp/Vint*(Pn*(-fintN) + Pd*phiNP(NbEp)*(-fintD)*exp(NbEp)) ...
         -Asm/Vint*(Pn*fintN + Pd*phiNP(Nsm)*fintD) ...
         -AimInt/Vint*(Pn*fintN + Pd*phiNP(NimInt)*fintD) ...
         -AbEd/Vint*(Pn*fintN + Pd*phiNP(NbEd)*fintD);
Kint_sm = -Asm/Vint*(Pn*(-fsmN) + Pd*phiNP(Nsm)*(-fsmD)*exp(Nsm));
Kint_imInt = -AimInt/Vint*(Pn*(-fimIntN) + Pd*phiNP(NimInt)*(-fimIntD)*exp(NimInt));
Kint_cEd  = -AbEd/Vint*(Pn*(-fcEdN) + Pd*phiNP(NbEd)*(-fcEdD)*exp(NbEd));
Kint_p = 0;
Sint = 0;
 
Ksm_aEp = 0; Ksm_imEp = 0; Ksm_cEp = 0;
Ksm_int = Asm/Vsm*(Pn*fintN + Pd*phiNP(Nsm)*fintD);
Ksm_sm  = Asm/Vsm*(Pn*(-fsmN) + Pd*phiNP(Nsm)*(-fsmD)*exp(Nsm));
Ksm_imInt = 0; Ksm_cEd = 0; Ksm_p = 0;
Ssm = 0;
 
KimInt_aEp = 0; KimInt_imEp = 0; KimInt_cEp = 0;
KimInt_int = AimInt/VimInt*(Pn*fintN + Pd*phiNP(NimInt)*fintD);
KimInt_sm  = 0;
KimInt_imInt = AimInt/VimInt*(Pn*(-fimIntN) + Pd*phiNP(NimInt)*(-fimIntD)*exp(NimInt));
KimInt_cEd = 0; KimInt_p = 0;
SimInt = 0;
 
KcEd_aEp = 0; KcEd_imEp = 0; KcEd_cEp = 0;
KcEd_int = AbEd/VcEd*(Pn*fintN + Pd*phiNP(NbEd)*fintD);
KcEd_sm = 0; KcEd_imInt = 0;
KcEd_cEd = AbEd/VcEd*(Pn*(-fcEdN) + Pd*phiNP(NbEd)*(-fcEdD)*exp(NbEd)) ...
         -AaEd/VcEd*(Pn*fcEdN + Pd*phiNP(NaEd)*fcEdD);
KcEd_p = -AaEd/VcEd*(Pn*(-fpN) + Pd*phiNP(NaEd)*(-fpD)*exp(NaEd));
ScEd = 0;
 
Kp_aEp = 0; Kp_imEp = 0; Kp_cEp = 0; Kp_int = 0; Kp_sm = 0; Kp_imInt = 0;
Kp_cEd = AaEd/Vp*(Pn*fcEdN + Pd*phiNP(NaEd)*fcEdD);
Kp_p   = AaEd/Vp*(Pn*(-fpN) + Pd*phiNP(NaEd)*(-fpD)*exp(NaEd));
Sp = 0;
 
M = [KaEp_aEp,KaEp_imEp,KaEp_cEp,KaEp_int,KaEp_sm,KaEp_imInt,KaEp_cEd,KaEp_p;...
     KimEp_aEp,KimEp_imEp,KimEp_cEp,KimEp_int,KimEp_sm,KimEp_imInt,KimEp_cEd,KimEp_p;...
     KcEp_aEp,KcEp_imEp,KcEp_cEp,KcEp_int,KcEp_sm,KcEp_imInt,KcEp_cEd,KcEp_p;...
     Kint_aEp,Kint_imEp,Kint_cEp,Kint_int,Kint_sm,Kint_imInt,Kint_cEd,Kint_p;...
     Ksm_aEp,Ksm_imEp,Ksm_cEp,Ksm_int,Ksm_sm,Ksm_imInt,Ksm_cEd,Ksm_p;...
     KimInt_aEp,KimInt_imEp,KimInt_cEp,KimInt_int,KimInt_sm,KimInt_imInt,KimInt_cEd,KimInt_p;...
     KcEd_aEp,KcEd_imEp,KcEd_cEp,KcEd_int,KcEd_sm,KcEd_imInt,KcEd_cEd,KcEd_p;...
     Kp_aEp,Kp_imEp,Kp_cEp,Kp_int,Kp_sm,Kp_imInt,Kp_cEd,Kp_p];
 
G = [SaEp,SimEp,ScEp,Sint,Ssm,SimInt,ScEd,Sp]';
 
% ---------- Sanity check ----------
if any(~isfinite(M(:))) || any(~isfinite(G(:)))
    [r,c] = find(~isfinite(M), 1);
    if ~isempty(r)
        fprintf('First non-finite M at (%d,%d): %g\n', r, c, M(r,c));
    end
    error('Non-finite entries detected in M or G (NaN/Inf). Check parameters/logP/pKa/z.');
end
 
%% ------------------ Solve ODE system ------------------
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
 
CaEp0 = dose/VaEp; % initial concentration in apical epithelium
X0 = [CaEp0;0;0;0;0;0;0;0];
 
tspan = [0 604800];   % 7 days in seconds
[t,Y] = ode15s(@LungModelOde, tspan, X0, options, M, G);
 
%% ------------------ PK metrics ------------------
j = find(Y(:,8) >= (0.5*dose/Vp), 1, 'first');
if isempty(j), j = numel(t); end
t_half_a = t(j);
 
idx_half_d = find(Y(:,1) <= (0.5*CaEp0), 1, 'first');
if isempty(idx_half_d), idx_half_d = numel(t); end
t_half_d = t(idx_half_d);
 
% ------------------ Lung burden definitions ------------------
% Tissue-associated burden excluding aEp
MLung_tissue = Y(:,2)*VimEp + Y(:,3)*VcEp + Y(:,4)*Vint + ...
               Y(:,5)*Vsm + Y(:,6)*VimInt + Y(:,7)*VcEd;
 
% Total lung burden including aEp
MLung_total = Y(:,1)*VaEp + Y(:,2)*VimEp + Y(:,3)*VcEp + Y(:,4)*Vint + ...
              Y(:,5)*Vsm + Y(:,6)*VimInt + Y(:,7)*VcEd;
 
% Peak burden metrics now based on total lung burden including aEp
[MLung_max,index_max] = max(MLung_total);
 
PMimEp_max  = max(Y(:,2)*VimEp)/dose*100;
PMcEp_max   = max(Y(:,3)*VcEp)/dose*100;
PMint_max   = max(Y(:,4)*Vint)/dose*100;
PMsm_max    = max(Y(:,5)*Vsm)/dose*100;
PMimInt_max = max(Y(:,6)*VimInt)/dose*100;
PMcEd_max   = max(Y(:,7)*VcEd)/dose*100;
 
% Lung burden at plasma half-appearance time, including aEp
MLung_h = Y(j,1)*VaEp + Y(j,2)*VimEp + Y(j,3)*VcEp + Y(j,4)*Vint + ...
          Y(j,5)*Vsm + Y(j,6)*VimInt + Y(j,7)*VcEd;
 
% Surface-layer fraction (aEp only)
PSL_h = Y(j,1)*VaEp/dose*100;
 
% Whole-lung burden including aEp
PLung_h = MLung_h/dose*100;
PMlung_max = MLung_max/dose*100;
TLung_max = t(index_max);
 
end
 
function y = phiNP_safe(N)
    % Stable phi(N) = N/(exp(N)-1) with correct limit phi(0)=1
    y = N ./ expm1(N);
    small = abs(N) < 1e-6;
    y(small) = 1 - N(small)/2 + (N(small).^2)/12;
    y(N==0) = 1;
end
