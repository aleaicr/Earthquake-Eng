function [mu_ln_SAH,sigma_ln_SAH] = campbell_1997_SAH(T,mu_ln_PGA,sigma_ln_PGA,Mw,r,S_SR,S_HR,D,f_SA_D)
% Calcular mu y sigma con GMM de campbell 1997 considerando los parámetros
% de Villaverde 2009

% Coeficientes de Villaverde 2009 tabla 7.3
T_villaverde = [0.05; 0.075; 0.1; 0.15; 0.2; 0.3; 0.5; 0.75; 1;1.5;2;3;4];
c1 = [0.05;0.27;0.48;0.72;0.79;0.77;-0.28;-0.108;-1.79;-2.65;-3.28;-4.07;-4.26];
c2 = [0;0;0;0;0;0;0.74;1.23;1.59;1.98;2.23;2.39;2.03];
c3 = [0;0;0;0;0;0;0.66;0.66;0.66;0.66;0.66;0.66;0.66];
c4 = [-0.0011;-0.0024;-0.0024;-0.0010;0.0011;0.0035;0.0068;0.0077;0.0085;0.0094;0.0100;0.0108;0.0112];
c5 = [0.000055;0.000095;0.000007;-0.00027;-0.00053;-0.00072;-0.00100;-0.00100;-0.00100;-0.00100;-0.00100;-0.00100;-0.00100];
c6 = [0.20;0.22;0.14;-0.02;-0.18;-0.40;-0.42;-0.44;-0.38;-0.32;-0.36;-0.22;-0.30];
c7 = [0;0;0;0;0;0;0.25;0.37;0.57;0.72;0.83;0.86;1.05];
c8 = [0;0;0;0;0;0;0.62;0.62;0.62;0.62;0.62;0.62;0.62];

% Interpolar si es necesario (por si piden periodos que no hay valores de ci
if ~isequal(T,T_villaverde)
    c1 = interp1(T_villaverde,c1,T,"linear","extrap");
    c2 = interp1(T_villaverde,c2,T,"linear","extrap");
    c3 = interp1(T_villaverde,c3,T,"linear","extrap");
    c4 = interp1(T_villaverde,c4,T,"linear","extrap");
    c5 = interp1(T_villaverde,c5,T,"linear","extrap");
end

% GMPE Campbell 1997 para la media y desviación de la distribución
% lognormal de la IM

% Convertir señales vectoriales a vectores
nT1 = ones(length(T),1);
mu_ln_PGA = mu_ln_PGA*nT1;
Mw = Mw*nT1;
r = r*nT1;
S_SR = S_SR*nT1;
S_HR = S_HR*nT1;
D = D*nT1;
f_SA_D = f_SA_D*nT1; %%% Revisar este valor (asumo que D > 1km)

% Media
mu_ln_SAH = mu_ln_PGA + c1 + c2.*tanh(c3.*(Mw-4.7*nT1)) + (c4 + c5.*Mw).*r + 0.5*c6.*S_SR + c6.*S_HR + c7.*tanh(c8.*D).*(nT1 - S_HR) + f_SA_D;

% Desviación estándar
sigma_ln_SAH = sqrt(sigma_ln_PGA.^2 + nT1*0.27^2);

end