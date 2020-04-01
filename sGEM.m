%-------------------------------------------------------------------------
%% Stochastic Gastric Emptying Model (sGEM)
%-------------------------------------------------------------------------
% Fuction to mimic the stomach emptying into the intestines

% Yokrattanasak J, De Gaetano A, Panunzi S, Satiracoo P, Lawton WM,
% Lenbury Y. A Simple, Realistic Stochastic Model of Gastric Emptying.
% Thomas DM, ed. PLoS ONE. 2016;11(4):e0153297.
% doi:10.1371/journal.pone.0153297.

function [stomach, param] = sGEM(param)
Tge = param.Tge;                        % Average, gastric emptying time (min)
dt = param.dt;
tspan = 0:dt:Tge;
n = length(tspan);

P = param.P;                            % Meal number of the day
phe0 = param.mealPKU(P);                % Phe at start of meal (mg)

mealtspan = param.mealtspan;            % Time between meals spaced by dt

stomach = ones(1,length(mealtspan));
a = 1E-6; b = 6;
dW = zeros(1, n); S = ones(1, n);
while S(n) > 0.001
    for j = 1:n
        dW(j) = sqrt(dt) * randn(1);
        W = cumsum(dW);
        PEW = exp(abs(W) .^ b .* -a);
        S(j) = min(PEW);
    end
end
stomach(1:n) = S .* phe0;               % Phe in stomach across gastric emptying (mg)
end