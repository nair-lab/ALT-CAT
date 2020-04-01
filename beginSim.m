%-------------------------------------------------------------------------
%% Altered for Living Therapeutics Compartmental Absorption and Transit (ALT-CAT) Model
%-------------------------------------------------------------------------
% An adapted compartmental absorption and transit model to estimate the
% rate of phenylalanine absorption with the presence of a enzyme-substitution
% therapy delivered with a probiotic. The model considers simultaneous small
% intestinal transit flow, phenylalanine conversion by a microbiome-bot
% theraputic, and phenylalanine absorption into the plasma. The stomach is
% modeled as spondaneous bolus injections into the duodennum.

function [final, intbuild, param, blovol] = beginSim
%-------------------------------------------------------------------------
%% User inputs
%-------------------------------------------------------------------------
runs = 1;         param.runs = runs; % number of simulations (remove Plot section if > 1)
age  = 9;         param.age  = age;  % 1 through 10 for diff age groups (see popGut)
enz  = 0;         param.enz  = enz;  % 0 through 3 to change PAL enzyme (see popGut)
therapy0 = 0;                        % change therapeutic coverage, alpha_c

runtime = 24 * 60;
dt   = 0.1;
%-------------------------------------------------------------------------
%% Definition of initial conditions
%-------------------------------------------------------------------------
runtspan = 0:dt:runtime;        param.runtime = runtime;
druntime = length(runtspan);    param.dt      = dt;

phe = 165.189;

plateR = 24.5 / 2;                      % 6-well plate diameter (mm); Holds 10^6 epithelial cells
plateA = pi * plateR^2 / 100;           % well area (cm2)
plateV = 14.7 / plateA / 1000000;       % Cell internal volume conversion (L/cm2)

mu = therapy0; sigma = 0.25;
therbuild = (exp(sigma .* randn(1,7) + mu))./100;
fprintf('Therapy set\n');
param.therapy = therbuild;              % Percent of cell population that are therapeutic

[intbuild, param] = popGut(param);      % Populating the intestines

phestart = param.mealPKU;               % Phe at start of meal (mg)
mealnum = param.mealnum;                % Number of meals over day
blovol = param.Vblo;                    % Blood volume (L)

sleep = param.sleep;                    % Hours of sleep over day
sleeptspan = 0:dt:sleep*60;
interval = (24 - sleep)/mealnum * 60;   % Digestion period for each meal (min)
param.interval = interval;

MLT = param.Mlt;                        % Cell concentration of therapy (cell/compt)
Vmax = param.Vmax;                      % Conversion of phe (mg/min.compt)
Km = param.Km;                          % Half max conc (mg/compt)

final = {};
for n = 1:runs
    param.run = n;
    day = []; 
    mealtspan = 0:dt:interval;
    param.mealtspan = mealtspan;
    param.Mlt = MLT(n,:);
    param.mealPKU = phestart(n,:);
    param.Vmax = Vmax(n,:);
    param.Km = (n);
    param.Vblo = blovol(n);
    
    fSA = (120 - 60) * rand + 60;       % Surface area scale factor in SI    
    param.Vint = plateV * intbuild.dimSum.sa(n) / 7 * fSA;
    param.Vlum = intbuild.SI.ComptVol(n) / 1000; % Compartment volume in L

    Cb0 = 0.5 * param.Vblo * phe;       % Phe in blood at start of day (mg)
    phelum0 = ones(1,7) .* 0.5 .* param.Vlum .* phe; % Phe in compartment at start of meal (mg)
    phelum0 = [phestart(n,1), phelum0, 0, 0, 0, Cb0];    
    for P = 1:mealnum + 1
        param.P = P;
        if (P <= mealnum && P > 1)
            initial.Sto = phestart(n,P) + phelum0(1);
        else
            initial.Sto = phelum0(1);
        end
        initial.SI1 = phelum0(2);       initial.SI2 = phelum0(3);
        initial.SI3 = phelum0(4);       initial.SI4 = phelum0(5);
        initial.SI5 = phelum0(6);       initial.SI6 = phelum0(7);
        initial.SI7 = phelum0(8);       initial.Col = phelum0(9);
        initial.tCA = phelum0(10);      initial.Excr = phelum0(11);
        initial.Blood = phelum0(12);
        if age > 5                      % Sleep time is significant
              if P > mealnum            % Sleep conversion iteration
                  mealtspan = sleeptspan;
                  param.mealtspan = sleeptspan;
                  phestart(n,P) = phelum0(1);
                  param.mealPKU(P) = phelum0(1);
              end
%-------------------------------------------------------------------------
%% Modeling the stomach
%-------------------------------------------------------------------------
            flag1 = 1;
            flag2 = 0;
            while (flag1 == 1 || flag2 == 0)
                [stobuild, param] = sGEM(param);
                sGE = spline(mealtspan, stobuild(1:length(mealtspan)));
                param.sGEprime = fnder(sGE);
%               [t,y] = ode23(@(t,x) stocheck(t,x,param), mealtspan, param.mealPKU(P));
                [t,y] = ode45(@(t,x) stocheck(t,x,param), mealtspan, param.mealPKU(P));
                flag1 = any(y > 1.05*param.mealPKU(P));
                flag2 = any(y(end) < 50);
            end
        end
        fprintf('Done modeling stomach\n');
%-------------------------------------------------------------------------
%% Solve the ODE system
%------------------------------------------------------------------------
        var_names = fieldnames(initial);
        y0 = [];
            for i = 1:length(var_names)
                y0 = [y0; initial.(var_names{i})];
            end
        opts = odeset('NonNegative',12);
%       [t,y] = ode23(@(t,x) aCAT(t,x,param), mealtspan, y0,opts);
        [t,y] = ode45(@(t,x) aCAT(t,x,param), mealtspan, y0,opts);
        fprintf('DONE \n');
        phelum0 = y(end,:);             % Next meal start concentration
        day = [day; y];
    end
    final{n} = day(1:druntime,:);
%-------------------------------------------------------------------------
%% Plot the results
%-------------------------------------------------------------------------
legend_texts = cell(12, 1);%cell(10,1); %
    for i = 1:12
        text = [var_names{i} '(t)'];
        legend_texts{i} = text;
    end
model_title = 'aCAT model solution space';
figure(n)
    p = plot(runtspan, final{n}(:,1:12));
    set(p, 'LineWidth', 3);
    set(gca, 'FontSize',18);
    xlabel('Time (min)');
    ylabel('Phe (mg)');
    title(model_title);
   
    legend(legend_texts);
end
end