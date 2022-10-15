%% Examen Final Ingeniería Sísmica
% Cristóbal Adasme - Alexis Contreras

%% Mala
% Ya que controlaba Qmin, Qmax, deltamax; solo había que controlar deltamax
% = 0.002

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
T1_step = -0.01;
T1_final = 0.35;
T1_rango = T1_final:T1_step:T1_init;

contador = 0;
EI_vect = zeros(length(T1_rango),1);

for i = 1:length(T1_rango)
    T1 = T1_rango(i);
    fprintf('T1 = %2.2f \n',T1)
    run p2_examen_sismica_parte6.m
    fprintf('Qmin = %2.3f \n',Qmin)
    fprintf('Qmax = %2.3f \n',Qmax)
    fprintf('Vb_SRSSmax = %2.3f \n',Vb_SRSSmax)
    fprintf('rdp counter = %2.f \n \n',counter)
    contador = contador + 1;
    
    if cond1 == 1 || cond2  == 1 || cond3 == 1
        EI_opt = EI;
        break
    end
    EI_opt = EI;
    T1_opt = T1;
    R_ast_opt = R_ast;
    Vb_opt = Vb_SRSSmax;
    Vmax = Qmax;
    Vmin = Qmin;
    FE_opt = FE;
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




