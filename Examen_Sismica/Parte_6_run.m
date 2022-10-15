%% Examen Final Ingeniería Sísmica
% Cristóbal Adasme - Alexis Contreras

%% Inicializar
clear variables
close all
clc

%% PREGUNTA 6
% Variamos Periodo del primer modo a buscar
% Condiciones: (0 si cumplen, si 1 no cumple -> se retorna periodo
% anterior)
% Condición 1: Vb_SRSSmax > Qmin
% Condición 2: Vb_SRSSmax < Qmax
% Condición 3: rdp_SRSSmax < delta

T1_init = 0.01;
T1_step = 0.01;
T1_final = 3;
T1_rango = T1_init:T1_step:T1_final;

contador = 0;

for i = 1:length(T1_rango)
    T1 = T1_rango(i);
    fprintf('T1 = %2.2f \n',T1)
    run p2_examen_sismica_parte6.m
    fprintf('Qmin = %2.3f \n',Qmin)
    fprintf('Qmax = %2.3f \n',Qmax)
    fprintf('Vb_SRSSmax = %2.3f \n',Vb_SRSSmax)
    fprintf('rdp counter = %2.f \n \n',counter)
    if cond3 == 1                                                           % cond1 == 1 || cond2 == 1 || cond3 == 1
        break
    end
    % Mientras no se rompa la condición 3, esta será la solución óptima
    EI_opt = EI;
    T1_opt = T1;
    R_ast_opt = R_ast;
    Vb_opt = Vb_SRSSmax;
    Vmax = Qmax;
    Vmin = Qmin;
    FE_opt = FE;
    rdp_SRSSmax_anterior = rdp_SRSSmax;
    u_SRSSmax_oficial_anterior = u_SRSSmax_oficial;
end

tabla = table();
tabla.EI_kNcm2 = EI_opt*100^2;
tabla.T1_sec = T1_opt;
tabla.R_ast = R_ast_opt;
tabla.Vb_kN = Vb_opt;
tabla.Vmax_kN = Vmax;
tabla.Vmin_kN = Vmin;
tabla.FE = FE;
disp(tabla)
clear tabla

figure
plot(u_SRSSmax_oficial_anterior*100,Pisos,'-o')
xlabel('Desplazamiento (U_{SRSS,max}) [cm]')
ylabel('Pisos')
grid on

delta = 0.002*ones(cant_pisos,1);
figure
plot(rdp_SRSSmax_anterior,Pisos,'-o')
hold on
plot(delta,Pisos,'--')
hold off
xlim([0.8*10^-3 2.1*10^-3])
grid on
xlabel('Desplazamiento')
ylabel('Pisos')
legend('rdp_{SRSS} \delta','\delta_{max} = 0.002')





