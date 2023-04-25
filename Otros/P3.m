%% Pregunta 3
% Considerando dos fuentes puntuales que solo generan una Magnitud de
% Momento Mw cada una.
%% Inicializar
clear variables
close all
clc

%% a)
% Grafique la curva de amenaza sísmica para el PGA. La curva debe graficarse
% en escala logarítmica y debe cubrir un rango de PGA entre 0.1g y 2.0g, 
% calculados cada 0.1g.

% Inputs fuentes
lam_m_A = 0.01;
lam_m_B = 0.002;
Mw_A = 6.5;
Mw_B = 7.5;
r_A = 10;
r_B = 20;

% Inputs GMPES   (revisar)
F_A = 0;  % Strike-slip
S_SR_A = 1; % ? ni idea
S_HR_A = 1; % hard-rock sites

F_B = 0;
S_SR_B = 1;
S_HR_B = 1;

% Rango de PGA
PGA_range  = (0.001:0.001:4).';  % en [g]   rango de PGA para la curva de amenaza

% GMPE utilizando Campbell 1997
[mu_ln_PGA_A,sigma_ln_PGA_A] = campbell_1997_AH(Mw_A,r_A,F_A,S_SR_A,S_HR_A);
[mu_ln_PGA_B,sigma_ln_PGA_B] = campbell_1997_AH(Mw_B,r_B,F_B,S_SR_B,S_HR_B);

% Calculo de probabilidades
probabilidad_A = zeros(length(PGA_range),1);
probabilidad_B = zeros(length(PGA_range),1);

for PGA_index = 1:length(PGA_range)
    PGA_value = PGA_range(PGA_index,1);
    probabilidad_A(PGA_index,1) = probs(PGA_value,mu_ln_PGA_A,sigma_ln_PGA_A);
    probabilidad_B(PGA_index,1) = probs(PGA_value,mu_ln_PGA_B,sigma_ln_PGA_B);
end

% Calculo de la amenaza
lambda_A = lam_m_A*probabilidad_A;
lambda_B = lam_m_B*probabilidad_B;
lambda_IM = lambda_A + lambda_B;

figure
loglog(PGA_range,lambda_IM)
hold on
loglog(PGA_range,lambda_A)
loglog(PGA_range,lambda_B)
hold off
grid on
xlabel('PGA [g]')
ylabel('\lambda_{IM}(IM>im) [1/yr]')
title('Curva de amenaza')
legend('A+B','A','B')

%% b Generar espectro de amenaza uniforme (UHS)
% Extra Inputs
S_SR_A = 1;
S_SR_B = 1;

S_HR_A = 1;
S_HR_B = 1; 

D_A = 2;
D_B = 2;

f_SA_D_A = 0;
f_SA_D_B = 0;

SAH_range = (0.001:0.001:4).';                % SAH, es la aceleración espectral Sa(T), esta es la nomenclatura que usa Campbell
T_range = [0.05; 0.075; 0.1; 0.15; 0.2; 0.3; 0.5; 0.75; 1;1.5;2;3;4];       % Rango de periodos a observar el espectro de amenaza uniforme

% mu y sigma son vectores donde cada fila corresponde a un periodo.
[mu_ln_SAH_A,sigma_ln_SAH_A] = campbell_1997_SAH(T_range,mu_ln_PGA_A,sigma_ln_PGA_A,Mw_A,r_A,S_SR_A,S_HR_A,D_A,f_SA_D_A);
[mu_ln_SAH_B,sigma_ln_SAH_B] = campbell_1997_SAH(T_range,mu_ln_PGA_B,sigma_ln_PGA_B,Mw_B,r_B,S_SR_B,S_HR_B,D_B,f_SA_D_B);

% Calcular probabilidades
probabilidad_A = zeros(length(SAH_range),1);
probabilidad_B = zeros(length(SAH_range),1);

for T_index = 1:length(T_range)
    for SAH_index = 1:length(SAH_range)
        SAH_value = SAH_range(SAH_index);
        probabilidad_A(SAH_index,T_index) = probs(SAH_value,mu_ln_SAH_A(T_index,1),sigma_ln_SAH_A(T_index,1));
        probabilidad_B(SAH_index,T_index) = probs(SAH_value,mu_ln_SAH_B(T_index,1),sigma_ln_SAH_B(T_index,1));
    end
end

lambda_SAH_A = lam_m_A*probabilidad_A;
lambda_SAH_B = lam_m_B*probabilidad_B;
lambda_SAH = lambda_SAH_A + lambda_SAH_B;

% Figuras
figure
for T_index = 1:length(T_range)
    loglog(SAH_range,lambda_SAH(:,T_index))
    hold on
end
loglog(SAH_range,1/475*ones(length(SAH_range),1),'linewidth',2,'color','k')
hold off
grid on
xlabel('Sa(T) [g]')
ylabel('\lambda_{IM}(IM>im) [1/yr]')
title('Curva de amenaza')
legend([convertStringsToChars("T = " + string(T_range) + " [s]");'\lambda_{IM} = 1/475'])

% % Generar UHS  ************************************** HE AQUÍ EL PROBLEMA, NO FUNCIONA INTERP1
% IM_UHS = zeros(length(T_range),1);
% for T_index = 1:length(T_range)
%     IM_UHS(T_index,1) = interp1(lambda_SAH(:,T_index),SAH_range,1/475,'linear','extrap');
% end
% % *************************************************** HE AQUÍ EL PROBLEMA

% Generar UHS  (alternativa a INTERP1, buscar desnivel menor)
IM_UHS = zeros(length(T_range),1);
desnivel_menor = 1000*ones(length(T_range),1);
index = 0;
for T_index = 1:length(T_range)
    for Sa_index = 1:length(SAH_range)
        desnivel = abs(1/475 - lambda_SAH(Sa_index,T_index)); 
        if desnivel < desnivel_menor(T_index)
            desnivel_menor(T_index) = desnivel;
            index = Sa_index;
        end
    end
    IM_UHS(T_index,1) = SAH_range(index);
end

% Figura espectro de amenaza unforme
figure
plot(T_range,IM_UHS)
xlabel('T [s]')
ylabel('Sa(T) [g]')
title('Espectro de Amenaza uniforme')
grid on
xlim([0 T_range(end,1)])

tabla = table();
tabla.T = T_range;
tabla.Sa = IM_UHS;
tabla.desnivel = desnivel_menor;
disp(tabla)


