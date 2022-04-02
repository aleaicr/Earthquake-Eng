%% Inicializar
clear variable
close all
clc

%% Comentarios
% CCC = Christmas Canyon China Lake
% CICCC_raw.txt : Registro Sísmico en crudo
% CICCC_acc_procesado, contiene el registro procesado, con correción de linea base

%% Importar datos

Regi_raw = readmatrix('CICCC_raw'); % g, dt = 0.01
Regi_acc_procesado = readmatrix('CICCC_acc_procesado.txt'); % cm/s2, dt = 0.02
Regi_disp_procesado = readmatrix('CICCC_disp_procesado'); % cm, dt = 0.02

%% Convertir Regi_raw a cm/s2
Regi_raw2 = Regi_raw;
Regi_raw  = Regi_raw2*981; % g = 9.81m/s2 = 981cm/s2
clear Regi_raw2

%% Procesamiento
mr = length(Regi_raw);                      % Tamaño del registro
mapr = length(Regi_acc_procesado);
mdpr = length(Regi_disp_procesado);

dtr = 0.01;                                 % Paso temporal del registro
dtapr = 0.02;
dtdpr = 0.02;

tr = (0:dtr:(mr-1)*dtr).';                  % Rango de Frecuencias
tapr = (0:dtapr:(mapr-1)*dtapr).';
tdpr = (0:dtdpr:(mdpr-1)*dtdpr).';

Fsr = 1/dtr;                                % Frecuencia de Muestreo
Fsapr = 1/dtapr;
Fsdpr = 1/dtdpr;

sr = -Fsr/2:Fsr/mr:(Fsr/2-Fsr/mr);          % Rango de frecuencias
sapr = -Fsapr/2:Fsapr/mapr:(Fsapr/2-Fsapr/mapr);    
sdpr = -Fsdpr/2:Fsdpr/mdpr:(Fsdpr/2-Fsdpr/mdpr);  

wr= 2*pi*sr;                                % Rango de frecuencias angulares
wapr = 2*pi*sapr; 
wdpr = 2*pi*sdpr;

Tr = 2*pi./wr;                              % Rango de periodos
Tapr = 2*pi./wapr;  
Tdpr = 2*pi./wdpr;  

%%4.1
% Gráfico Ac. Vs Tiempo
figure
plot(tr,Regi_raw)
xlabel('Acc [cm/s2]')
ylabel('t [s]')
title('Raw vs t (Dom.Tiempo)')

figure
plot(tapr,Regi_acc_procesado)
xlabel('t [s]')
ylabel('Acc [cm/s2]')
title('Acc_Procesado vs t (Dom.Tiempo)')

figure
hold on
plot(tr,Regi_raw)
plot(tapr,Regi_acc_procesado)
xlabel('t [s]')
ylabel('Acc [cm/s2]')
title('Raw vs Acc_Procesado (Dom.Tiempo)')
hold off

%% 4.2
% Gráfico Ac. vs Frecuencia
figure
plot(sr,Regi_raw)
ylabel('Acc [cm/s2]')
xlabel('Frec [hz]')
title('Raw vs s (Dom.Frecuencia)')

figure
plot(sapr,Regi_acc_procesado)
xlabel('Frec [s]')
ylabel('Acc [cm/s2]')
title('Acc_Procesado vs s (Dom.Frecuencia)')

figure
hold on
plot(sr,Regi_raw)
plot(sapr,Regi_acc_procesado)
xlabel('Frec [hz]')
ylabel('Acc [cm/s2]')
title('Raw vs Acc Procesado (Dom.Frecuencia)')
hold off

%% 4.3
% Integrar dos veces el registro crudo de acleraciones para obtener el registro de desplazamientos

%% 4.3.1 Aproximación desplazamiento medieante trapecios
Regi_raw_vel = cumtrapz(Regi_raw);
Regi_raw_disp = cumtrapz(Regi_raw_vel);

%% 4.3.2 u(T) = i.Fourier(-F(upp)/w^2)
disp_fourier = ifft(-fft(Regi_raw)./wr.^2);
UPPw = fft(Regi_raw);
UPPw2 = UPPw./wr.^2
furier_asd = ifft(-Uppw2);


figure
hold on
plot(sr,disp_fourier)
plot(sr,Regi_raw_disp)
plot(sdpr,Regi_disp_procesado)
hold off

%% NOTE TENGO 18GB de RAM para poder 