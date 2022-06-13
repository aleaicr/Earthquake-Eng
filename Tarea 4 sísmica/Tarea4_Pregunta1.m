%% Tarea 4 Pregunta 1
% Ingeniería Sísmica
% Alexis Contreras R. - Cristobal Adasme 


%% Inicializar
clear variables
close all
clc

%% Datos
load('el_centro.mat');                                                      % dt = 0.0200; uppG = 1560x1 double
load('resultados_P1.mat');                                                  % fs_normalizado; t; u

% Renombramos
u_resultados = u;

% Parámetros
alfa = -0.15;
xi = 0.05;
Fy = 0.08;
Tn = 2;
R = 10^(-10);
T = transpose(0:dt:(size(uppG,1)-1)*dt);
beta = 1/6;

% Condiciones iniciales
u0 = 0;
u0p = 0;

%% PARTE (i)
[u,ud,udd,fs,fluencia,ductilidad] = Newmark_No_Lineal_Bilineal(alfa,beta,Tn,xi,dt,u0,u0p,uppG,Fy,R,T);
figure
hold on
plot(T,u*980)
plot(t,u_resultados)
% A = max(abs(u))*980;
% B = max(abs(u_resultados));
% (A-B)/B;
xlabel('Tiempo [s]')
ylabel('Desplazamiento [cm]')
legend("Aproximación Numérica","Real")
hold off

%% PARTE (ii)
figure
hold on
plot(u*980,fs*980)                                                          
plot(u_resultados,fs_normalizada)
xlabel('Desplazamiento [cm]')
ylabel('Fuerza Restitutiva')
legend("Aproximación Numérica","Real")
hold off
