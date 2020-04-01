function [intbuild, param] = popGut(param)
%-------------------------------------------------------------------------
phe    = 165.189;                   % mol wt phe (mg/mmol)
cellwt = 280 * 1E-12;               % dry cell wt E coli (mg)
Peff = 5E-4 * 60;                   % Observed effective permeability of phe (cm/min)

runs = param.runs;
age  = param.age;
enz  = param.enz;
therapy = param.therapy;

%-------------------------------------------------------------------------
%% Intestine data
%-------------------------------------------------------------------------
S = table;
S.section = ...
    {'Duodenum'; 'Jejenum'; 'Ilium'};
S.percent = ...
    [0.05; 0.4; 0.55];
S.bounds = ...                      % Cell pop boundaries in 0.1mL segments
    [1E4 1E5; 1E5 1E6; 1E6 1E9];

%-------------------------------------------------------------------------
%% Enzyme data
%-------------------------------------------------------------------------
E = table;
E.enzyme = {                        % Enzyme name
    'PaPAM'; 'AvPAL'; 'PcPAL'
    };
E.kcat = [                          % Turnover (1/s)
    0.08736; 4.4; 22
    ];
E.size = [                          % Molecular size (kDa) ****CHECK PcPAL****
    59; 62; 38.9
    ];
E.Km = [                            % Half max substrate concentration (mM)
    0.032; 0.06; 0.015
    ];
E.wtr = ...                         % Assumed enzyme expression (mg pro/mg cell)
    zeros(3,1);
    for i = 1:3
        E.wtr(i) = cellwt * 0.55 * 0.3;
    end

if enz == 0
    Jenz = 3.5 * phe / 1000 / 60 / 4 / 10^9; % Reported activity by Synlogic (mg/min.cell)
    Km = 1; param.kcat = 1;         % Spacer value not used in calculations
else
    Km = E.Km(enz) / 1000 * phe;    % Converted to mg/mL
    param.kcat = E.kcat(enz) * 60;  % Converted to 1/min
    Jenz = E.wtr(enz) / E.size(enz) * (phe / 1000) * param.kcat / 4;
end

%-------------------------------------------------------------------------    
%% Age data
%-------------------------------------------------------------------------
A = table;
A.group = {                         % Age range for stats
    'newborn'; '3 day'; '1 week'; '1 month'; '1-6 month'
    '7-12 month'; '1-6 year'; '8-13 year'; 'adult'; 'max adult'
    };
A.length = [                        % Average and stdev, small intesine length (cm)
    246.5 56.2; 246.5 56.2; 246.5 56.2; 246.5 56.2; 351.4 97.8;
    359.5 75.7; 423.0 99.0; 464.2 103.1; 595.5 101.8; 1000 200
    ];

D = table;
mu = A.length(age,1);
sigma = A.length(age,2);
D.length = sigma .* randn(runs,1) + mu;

A.id = [                            % Average, small intestine internal diameter (cm)
    1.5; 1.5; 1.5; 1.5; 1.5; 1.5; 2.0; 2.5; 2.5; 3
    ];

r = A.id(age) / 2;

A.weight = [                        % Average and stdev, person weight (kg)
    3.4 0.2; 3.4 0.2; 3.6 0.2; 4.4 0.3; 6.2 0.3;
    7.9 0.4; 17.5 1.7; 35 5.0; 71.0 5.1; 96.0 12.2
    ];

mu = A.weight(age,1);
sigma = A.weight(age,2);
weight = sigma .* randn(runs,1) + mu;
param.Vblo = 0.075 .* weight;

A.mealnum = [                       % Average, number of meals per day
    12; 12; 10; 8; 8; 6; 5; 4; 4; 3
    ];

A.pheReg = [                        % Average, phe in regular diet (mg/kg.day)
    83.3; 83.3; 83.3; 58.0; 40.0; 31.0; 46.0; 14.7; 9.3; 11.0
    ];

mealReg = A.pheReg(age) * A.weight(age,1) / A.mealnum(age);

A.phePKU = [                        % Average and stdev, phe in PKU diet (mg/kg.day)
    58 18; 58 18; 58 18; 58 18; 40 10; 31 8.5; 20 1; 14 1; 14 1; 9.3 1   
    ];

mu = A.phePKU(age,1);
sigma = A.phePKU(age,2);
mealPKU = zeros(runs,A.mealnum(age));
for i = 1:A.mealnum(age)
    mealPKU(:,i) = sigma .* randn(runs,1) + mu;
end

A.Tge = [                           % Average, gastric half emptying time (min)
    35; 35; 35; 48; 60; 90; 90; 90; 120; 120
    ];

Tge = A.Tge(age) * 2;

A.sleep = [                         % Average, time slept at night (h)
    0; 0; 0; 0; 0; 10; 11; 11; 9; 8
    ];  

m = runs;
n = height(S);

D.vol = zeros(m,n);                 % Total volume of small intestine (mL)
    for j = 1:n
        for i = 1:m
            D.vol(i,j) = S.percent(j) * D.length(i);
            D.vol(i,j) = round(pi * r^2 * D.vol(i,j), 1);
        end
    end
    
D.sa = zeros(m,1);                  % Cylindrical surface area of small intestine (cm2)
    for i = 1:m
        D.sa(i) = ...
        (2 * pi * r/2 * D.length(i)) + (2 * pi * r^2);
    end

%-------------------------------------------------------------------------
%% Sectioning the small intestines
%-------------------------------------------------------------------------
% Populating the small intestines with APT
SI = {};
for i = m:-1:1
    for j = 1:n
        k = D.vol(i,j) * 10;        % Number of 0.1mL segments in each section
        lb = S.bounds(j,1);
        ub = S.bounds(j,2);
        SI{i,j} = (ub - lb) .* rand(k,1) + lb; % Random population of bacteria
    end
    SI{i,4} = [
        SI{i,1}; SI{i,2}; SI{i,3}
        ];

% Compartmentalizing the small intestines
    compt = SI{i,4}; SI{i,5} = compt;
    p = floor(length(compt) / 7);   % Evenly dividing the SI into 7 compartments
    SI{i,6} = p / 10;
    for q = 6:-1:1
        cellnum(i,1) = sum(compt(1:p*q)) * therapy(q);
        Vmax(i,1) = cellnum(i,1) .* Jenz;
        cellnum(i,q+1) = sum(compt((p*q)+1:p*(q+1))) * therapy(q+1);
        Vmax(i,q+1) = cellnum(i,q+1) .* Jenz;
    end
    SI{i,5} = cellnum(i,:);         % Cells in each compt
    SI{i,7} = Vmax(i,:);            % Converted Vmax (mg phe/min.compt)
    SI{i,8} = Km * SI{i,6};         % Converted Km (mg phe/compt)
end
SI = cell2table(SI,'VariableNames',{
    'Duodenum','Jejunum','Ileum','FullSI','Cells','ComptVol','Vmax','Km'
    });

%-------------------------------------------------------------------------
%% Definition of model parameters
%-------------------------------------------------------------------------
if age > 5                          % Average, small intestinal transit time (min)
    Tsi = 3.5 * 60;
else
    Tsi = 1.5 * 60;
end
param.Tge = Tge;

param.ksi = 7 / Tsi;                % Average, SI transit rate (1/min)
param.kabs = 2 * Peff / r;          % Observed absorbance rate of phe (1/min)
param.kexcr = 0.0452;               % Reported excretion rate of cinnamic acid (1/min)

param.Vmax = SI.Vmax;
param.Km = SI.Km;

param.mealPKU = mealPKU .* weight / A.mealnum(age); % Converted PKU diet (mg/meal)
param.mealReg = mealReg;
param.mealnum = A.mealnum(age);

param.sleep = A.sleep(age);
param.Mlt = SI.Cells;

intbuild = struct;
intbuild.sectSum = S;
intbuild.agesSum = A;
intbuild.dimSum = D;
intbuild.SI = SI;

fprintf('Done populating gut\n');
end