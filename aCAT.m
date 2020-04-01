%-------------------------------------------------------------------------
%% Advanced Compartmental Absorption and Transit Model (aCAT)
%-------------------------------------------------------------------------
% Function to calculate derivatives of mass transfer equations

% Lawrence X. Yu and Gordon L. Amidon. A compartmental absorption and
% transit model for estimating oral drug absorption. International
% Journal of Pharmaceutics 186 (1999) 119-125.

function deriv = aCAT(t, x, param)
%--------------------------------------------------------------------------
%% Definitions
%--------------------------------------------------------------------------
Mlt     = param.Mlt;        ksi   = param.ksi;      comptVmax = param.Vmax;        
age     = param.age;        kge = 1 / param.Tge;    comptKm   = param.Km;
enz     = param.enz;        kabs  = param.kabs;     
run    = param.run;         kexcr = param.kexcr;    
                            kcat = param.kcat;
phe = 165.189;
tspan = param.interval;

Vi   = param.Vint;
Vb   = param.Vblo;
Vlum = param.Vlum;  

Sto = x(1);
SIcompt = zeros(7,1);
    for i = 1:7
        SIcompt(i) = x(i+1);
    end
Col = x(9); tCA = x(10); Excr = x(11); Blo = x(12);

%% Stomach
if age > 5                              % Gastric emptying stochastic model
    sGEprime = param.sGEprime;
    dSto_dt = ppval(sGEprime, t);
else
    dSto_dt = -kge * Sto;               % Linear gastric emptying
end
%% Flux and conversion factors
JuAPmax = 16.1;     JuBLmax = 1.07;     Cimax = 13.26; % Jflux (nmol/mgprotein.min)
KmAP    = 2.7;      KmBL    = 0.18;     KmCi  = 4.203; % C and Km (mM)
fuAP    = 0.2;      fuBL    = 0.2;      Kmt = 0.00072; % f constant (nmol/mgprotein.min/mM)
feAP    = 0.289;    feBL    = 0.476;
 
Clum = zeros(7,1);
    for i = 1:7
        Clum(i) = SIcompt(i) / Vlum / phe;
    end
Cblo = Blo / Vb / phe;

Jap = zeros(7,1);     Jtmax  = zeros(7,1);	Jconv = zeros(7,1); 
Jbl = zeros(7,1);     Jtrans = zeros(7,1);	conv  = zeros(7,1); % Jtrans (mmol/30s.mgwetcell)
Cinter  = zeros(7,1); alpha  = zeros(7,1);

for i = 1:7
    Cinter(i) = ...
        Cimax * Clum(i) / (KmCi + Clum(i));
    Jap(i) = ...
        (JuAPmax * Clum(i) / (KmAP + Clum(i)) + fuAP * Clum(i)) - feAP * Cinter(i);
    Jbl(i) = ...
        feBL * Cinter(i) - (JuBLmax * Cblo / (KmBL + Cblo) + fuBL * Cblo);
    alpha(i) = ...
        Jbl(i) / Jap(i);    
    Jtmax(i) = ...                      % Convert to mg/min.compt
        0.75 * 2 * 10^-15 * Mlt(i) * phe; 
    Jtrans(i) = ...
        Jtmax(i) * Clum(i) / (Kmt + Clum(i));
    if enz == 0
        Jconv(i) = comptVmax(i);
    else
        Jconv(i) = ...
        kcat * Clum(i) * Vlum * phe;
    end
    if Jtrans(i) > Jconv(i)
        conv(i) = Jtrans(i);
    elseif enz == 0
        conv(i) = comptVmax(i);
    else
        comptKm = comptKm / phe / Vlum; % Convert from mg/compt to mM
        conv(i) = comptVmax(i) * Clum(i) / (comptKm + Clum(i));
    end
end

%% Small Intestine
dSBL_dt = zeros(7,1);
dSBL_dt(1) = ...
    -dSto_dt - ksi * SIcompt(1) - alpha(1) * kabs*SIcompt(1) - conv(1);
for i = 2:7
    dSBL_dt(i) = ...
        ksi * (SIcompt(i-1) - SIcompt(i)) - alpha(i) * kabs * SIcompt(i) - conv(i);
end

%% Colon
dCol_dt = ksi * SIcompt(7) - kabs * Col;

%% Absorption into plasma
blophe = 0;
    for i = 1:7
        blophe = alpha(i) * SIcompt(i) + blophe;
    end
dBlo_dt = kabs * blophe + kabs * Col;

%% Conversion by therapy
dtCA_dt = sum(conv) - kexcr * tCA;
dExcr_dt = kexcr * tCA;

%% Output
deriv = [dSto_dt; dSBL_dt; dCol_dt; dtCA_dt; dExcr_dt; dBlo_dt];

%% Status
if t == 0
    fprintf('Beginning meal %d of run %d...',param.P,run);
end
end