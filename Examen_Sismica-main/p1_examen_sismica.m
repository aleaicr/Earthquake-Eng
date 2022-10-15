%% Examen Final Ingeniería Sísmica
% Cristobal Adasme - Alexis Contreras R.

%% Inicializar
clear variables
close all
clc

%% Pregunta 1
g = 1; % cm/s2  %si g = 1 se obtienen resultados con g?
R = 10;    % km                                                             % Distancia del sitio a fuente sismogénica (modelado puntual)
lambda_5 = 0.002;                                                           % Frecuencia de ocurrencia por año

% Parámetro Gutenberg-Richter
b = 1.0;
beta = 2.303;
Mmax = 8.0;
m0 = 5.0;

% Modelo de movimiento s´simico de Cornell et al (1979)
% mu_lnPGA = -0.152 + 0.859*M - 1.803*ln(R+25)
% sigma_lnPGA = 0.57;

%% Estimar CURVA DE AMENAZA SÍSMICA del PGA en el sitio

% Calcular la frecuencia anual de excedencia lambda_PGA para PGA = 0.1g a
% 2.0g

Delta_PGA = 0.1*g;
PGA_init = 0.10*g;
PGA_final = 2.0*g;
PGA_vect = (PGA_init:Delta_PGA:PGA_final)';


% Para integración numérica
Delta_m = 0.1;
% nM = (mmax-m0)/Delta_m;
M_vect = (m0:Delta_m:Mmax)'; % length(M_vect = nM)                             % Utilizar momentos 
nM = length(M_vect);
lambdas = zeros(length(PGA_vect),1);
R_vect = 10;                                                                   % Fijo, para todas las integraciones
nR = length(R_vect);
sigma_lnPGA_cornell = 0.57;

for i = 1:length(PGA_vect)
    % N° de fuentes = 1 -> se elimina la prima sumatoria
    % for i = 1:nFuentes
    sum_nM = 0;
    PM = zeros(nM,1);
    PIMxmr = zeros(nM,1);
    for j = 1:nM
        % mj = M_vect(j);
        % Solo hay una distancia -> nR = 1 y P(Ri = 10km) = 1              
        % for k = 1:nR
        %   rk = R_vect(k);
        mu_lnPGa_cornell = -0.152 + 0.859*M_vect(j,1) - 1.803*log(R+25);        % Cornell et al (1979), todos tienen el mismo r = 10km
%         sigma_lnPGA_cornell = 0.57;                                           % Cornell et al (1979)
        [CDF,PDF] = gutenberg_richter_acotada(M_vect(j,1),beta,m0,Mmax);
        PM(j,1) = CDF;                                                          % Ley recurrencia Gutenberg-Richter acotada (arriba y abajo)
        PIMxmr(j,1) = 1 - normcdf((log(PGA_vect(i,1))-mu_lnPGa_cornell)/sigma_lnPGA_cornell); % IP(IM>x|mj,rk)
        sum_nM = sum_nM + PIMxmr(j,1)*PM(j,1);
    end
    lambdas(i,1) = lambda_5*sum_nM;                                             % lambda(IM>x) = lambda(Mi>m0)*Sum_nM(P(IM>x|mj,rk)_normcdf*P(Mi=mj)_gut_richt*P(Ri = rk)_=1)
end

% Lambdas
% = lambda(IM>x) = lambda(PGA> x*g), ordenados

figure
loglog(PGA_vect,lambdas)
xlabel('PGA [g]')
ylabel('\lambda(PGA > x[g])')
title('Curva de Amenaza Sísmica')
grid on

tabla = table();
tabla.PGA_g = PGA_vect;
tabla.lambda = lambdas;
disp(tabla)
