%% Tarea 3 - Pregunta 2 - Ingeniería Sísmica
% Cristóbal Adasme |  Alexis Contreras

%% Inicializar

clear variables
close all
clc

%% Importar Registro
% load() guarda el registro en variables uppG (Vector de aceleraciones) y
% dt (paso temporal)

load('el_centro.mat');
n = length(uppG);                                                           % 15601

%% Datos

% Parámetros
Tn = 0.5;                                                                   % Periodo 
beta = 1/4;                                                                 % Factor de método newmark-beta
xi = 0;                                                                     % Razón de amortiguamiento
wn = 2*pi/Tn;                                                               % Frecuencia natural
kt = wn^2;                                                                  % kt

% Condiciones iniciales
ui = 0;
udi = 0;

Sd = Newmark_Lineal(beta,xi,dt,ui,udi,uppG);

for i = 1:size(Sd,1)                                                        
    if Sd(i,1) == Tn                                                        % Si estoy en el periodo
        u_max = Sd(i,2);                                                    % Desplazamiento máximo del oscilador de periodo Tn
    end
end

t = transpose(0:dt:(length(uppG)-1)*dt);                                    % Vector de tiempos
f0 = kt*u_max;                                                              % Fuerza resistente máxima
Fy = 0.125*f0;                                                              % Resistencia lateral
R = 10^(-5);                                                                % Error admisible

%% Parte /I) y (II)
[u,ud,udd,fs,fluencia,ductilidad] = Newmark_No_Lineal(beta,Tn,xi,dt,ui,udi,uppG,Fy,R,t);

g = 386.06;  %in/s2
u = u*g;       

% Graficamos (con subplot para que se paresca a lo del libro de chopra)
figure
hold on
subplot(3,1,1)
plot(t(:,1),u(:,1))
title('Desplazamientos')
ylabel('Desplazamientos [in]')
xlim([0 10])

subplot(3,1,2)
plot(t(:,1),fs(:,1))
title('Fuerza Restitutiva')
ylabel('f_s/W')
xlim([0 10])

subplot(3,1,3)
plot(t(:,1),fluencia(:,1))
title('Fluencia') 
xlabel('Tiempo [s]')
ylabel('Fluencia')
xlim([0 10])
hold off

% figure
% plot(t(:,1),u(:,1))
% title('Desplazamientos')
% xlabel('Tiempo [s]')
% ylabel('Desplazamientos [in]')
% 
% figure
% plot(t(:,1),fs(:,1))
% title('Fuerza Restitutiva')
% xlabel('Tiempo [s]')
% ylabel('f_s/W')
% 
% figure
% plot(t(:,1),fluencia(:,1))
% title('Fluencia') 
% xlabel('Tiempo [s]')
% ylabel('Fluencia')


%% Parte (III)
j = 0;
rango = [1.368; 2.05];
for i = 1:length(u)
    if t(i,1) <= rango(2,1) && t(i,1) >= rango(1,1)
        j = j + 1;
        fs_new(j,1) = fs(i,1);
        u_new(j,1) = u(i,1);
    end
%     if t(i,1) > rango(2,1)
%         break
%     end
end

figure
plot(u_new(:,1),fs_new(:,1))
title('Fuerza Restitutiva vs Desplazamientos') 
xlabel('Desplazamientos [in]')
ylabel('f_s/W')
