%% Tarea 3 - Pregunta 1 - Ingeniería Sísmica
% Cristóbal Adasme |  Alexis Contreras

%% Inicializar

clc
clear varibales
close all

%% Pregunta 1

%% Importar datos
Registro = importdata('Northridge-Tarzana-EW.txt', ' ',1);                  % Registro de Northridge 1994
dt = 0.02;                                                                  % dt sampling: Pasto Temporal

%% Datos / Variables
xi = [0; 0.01; 0.05: 0.10; 0.20];                                           % Razones de amortiguamiento
n = length(xi);                                                             % Cantidad de razones de amortiguamiento para analizar
beta = 1/4;                                                                 % Coeficiente del método de Newmark beta

eje_y = ["Espectro de Desplazamientos [m]","Espectro de Velocidades [m/s]","Espectro de Aceleraciones [m/s^2]","Pseudo-Espectro de Velocidades [m/s]","Pseudo-Espectro de Aceleraciones [m/s^2]"];
titulazos = ["Espectro de Desplazamientos","Espectro de Velocidades","Espectro de Aceleraciones","Pseudo-Espectro de Velocidades","Pseudo-Espectro de Aceleraciones"];
%%  Condiciones Iniciales
ui = 0;
udi = 0;

%% Aplicación del método de Newmark y graficar
% Newmark_Lineal es una función en Newmark_Lineal.m 
for j = 1:5
    figure (j)
    hold on
    for i = 1:n
        [Sd,Sv,Sa,Psv,Psa] = Newmark_Lineal(beta,xi(i,1),dt,ui,udi,Registro.data); % Ineficiente, calcula cuando no necesita
        Spectr = [Sd Sv Sa Psv Psa];
        plot(Spectr(:,2*(j-1)+1),Spectr(:,2*(j-1)+2))
    end
    hold off
    xlabel('Tiempo (t) [sec]')
    ylabel(eje_y(j))
    title(titulazos(j))
end

%% Transformada de Fourier

[UPPG,T] = T_Fourier(0.02,Registro.data);

UG = -UPPG./(2*pi./T).^2;
UPG = UPPG.*2.*pi./T.*sqrt(-1);

% Acotamos los datos a los que se quieren analizar
count = 1;
rango = [0; 5];
for i = 1:length(T)
    if (T(i,1))>rango(1,1) && (T(i,1))<=rango(2,1)
        UG_aco(count,1) = UG(i,1);
        UPG_aco(count,1) = UPG(i,1);
        UPPG_aco(count,1) = UPPG(i,1);
        T_aco(count,1) = T(i,1);
        count = count+1;
    end
%     if T(i,1)>5
%         break
%     end
end

%% Ordenar Fourier -> Para graficar
UG_new = fftshift(UG_aco);                                                  % Para dar vuelta la transformada
UPG_new = fftshift(UPG_aco);                                                 % y que coincida con el rango de periodo
UPPG_new = fftshift(UPPG_aco);                                          

%% Obtenición de espectro por Newmark Lineal
xi_0 = xi(1,1);                     
[Sd,Sv,Sa,Psv,Psa] = Newmark_Lineal(beta,xi_0,dt,ui,udi,Registro.data);     

%% Graficamos
% Para desplazamiento
max_UG = max(abs(UG_new));
Ug_norm = UG_new/(max_UG/max(abs(Sd(:,2))));
figure
hold on
plot(abs(T_aco),abs(Ug_norm))
plot(Sd(:,1),Sd(:,2))
xlabel('Periodo T [s]')
ylabel('Desplazamiento [cm]')
legend('F|(\omega)|','Sd')
hold off

% Para velocidades
max_UPG = max(abs(UG_new));
UPG_norm = UG_new/(max_UPG/max(abs(Sv(:,2))));
figure                           
hold on
plot(abs(T_aco),abs(UPG_norm))
plot(Sv(:,1),Sv(:,2))
xlabel('Periodo T [s]')
ylabel('Valocidad [cm/s]')
legend('F|(\omega)|','Sv')
hold off

% Para aceleraciones
max_UPPG = max(abs(UPPG_new));
UPPG_norm = UPPG_new/(max_UPPG/max(abs(Sa(:,2))));
figure
hold on
plot(abs(T_aco),abs(UPPG_norm))
plot(Sa(:,1),Sa(:,2))
xlabel('Periodo T [s]')
ylabel('Aceleración [cm/s^2]')
legend('F|(\omega)|','Sa')
hold off