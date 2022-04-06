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

Regi_raw  = Regi_raw*981; % g = 9.81m/s2 = 981cm/s2

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
xlabel('Aceleración [cm/s2]')
ylabel('t [s]')
title('Raw vs t')

figure
plot(tapr,Regi_acc_procesado)
xlabel('t [s]')
ylabel('Aceleración [cm/s2]')
title('Procesado vs t')

figure
hold on
plot(tr,Regi_raw)
plot(tapr,Regi_acc_procesado)
xlabel('t [s]')
ylabel('Aceleración [cm/s2]')
title('Raw vs Procesado (Dom.Tiempo)')
hold off

%% 4.2
% Gráfico Ac. vs Frecuencia

fourier_raw = abs(fft(Regi_raw)).^2;
fourier_procesado = abs(fft(Regi_acc_procesado)).^2;

figure
plot(sr,fftshift(fourier_raw))
ylabel('Aceleración [cm/s2]')
xlabel('Frecuencia [hz]')
title('Raw vs s')

figure
plot(sapr,fftshift(fourier_procesado))
xlabel('Frecuencia [s]')
ylabel('|H_{up}(\omega)|^2')
title('Procesado vs s')

figure
hold on
plot(sr,fftshift(fourier_raw))
plot(sapr,fourier_procesado)
xlabel('Frecuencia [hz]')
ylabel('|H_{up}(\omega)|^2')
title('Raw vs Procesado (Dom.Frecuencia)')
hold off

%% 4.3
% Integrar dos veces el registro crudo de acleraciones para obtener el registro de desplazamientos

%% 4.3.1 Aproximación desplazamiento medieante trapecios

Regi_raw_vel = cumtrapz(dtr,Regi_raw);
Regi_raw_disp = cumtrapz(dtr,Regi_raw_vel);

Regi_proc_vel = cumtrapz(dtapr,Regi_acc_procesado);
Regi_proc_disp = cumtrapz(dtapr,Regi_proc_vel);


figure
plot(tr,Regi_raw_vel)
xlabel('Tiempo [s]')
ylabel('Velocidad [cm/s]')
title('Velocidad integrando Aceleración Raw')

figure
plot(tr,Regi_raw_disp)
xlabel('Tiempo [s]')
ylabel('Desplazamiento [cm]')
title('Desplazamiento integrando Velocidad desde Raw')

figure
plot(tapr,Regi_proc_vel)
xlabel('Tiempo [s]')
ylabel('Velocidad [cm/s]')
title('Velocidad integrando Aceleración Procesada')

figure
plot(tapr,Regi_proc_disp)
xlabel('Tiempo [s]')
ylabel('Desplazamiento [cm]')
title('Desplazamiento integrando Velocidad desde Procesada')

%% 4.3.2 Transformada Inversa de Fourier
% u(t) = F^-1(-F(upp)/w^2)

% Hay que arreglar lo de abajo, no me funciona

disp_fourier = ifft(-fftshift(fft(Regi_raw))./(wr.^2).');
% disp_fourier_proc = ifft(-fft(Regi_acc_procesado)./wr.^2);

figure
plot(sr,real(disp_fourier))
xlabel('t [s]')
ylabel('u(t)')
title('Desplazamiento obtenido desde transformada de fourier del Registo')

% figure
% plot(tapr,disp_fourier_proc)
% xlabel('t [s]')
% ylabel('u(t)')
% title('Desplazamiento obtenido desde transformada de fourier del Registo')