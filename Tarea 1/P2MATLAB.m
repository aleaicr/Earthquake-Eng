%% Inicializar
clear variables
close all
clc

%% Importar datos
% Todos los datos tienen paso temporal de 0.005
% Todos tienen unidad cm/s2 para la aceleración del registro excepto pica2005_ew.txt

Concepcion = readmatrix("Concepcion2010-L.txt");
mc = length(Concepcion);

Santiago = readmatrix("stgocentro2010-L.txt");
ms = length(Santiago);

Talca = readmatrix("Talca2010-L.txt");
mt = length(Talca);

Valparaiso = readmatrix("valparaisoUTFSM2010-L.txt");
mv = length(Valparaiso);

Pica = readmatrix("pica2005_ew.txt");
mp = length(Pica);

Northridge = readmatrix("Northridge-Tarzana-EW.txt"); 
mn = length(Northridge);

%% Vector de tiempo;

dt = 0.005; % general para todos
tc = 0:dt:(mc-1)*dt; %Vector de tiempos para Concepción
ts = 0:dt:(ms-1)*dt; % Vector de tiempos para Santiago
tt = 0:dt:(mt-1)*dt; % Vector de tiempos para TAlca
tv = 0:dt:(mv-1)*dt; % Vector de tiempos para valparaíso
tp = 0:dt:(mp-1)*dt; % Vector de tiempos para Pica
tn = 0:dt:(mn-1)*dt; % Vector de tiempos apra Northridge - Tarzana

%% Conversión a cm/s2 para Pica

Pica_cms2 = Pica*981;  % g = 9.81 m/s2 = 981 cm/s2

%% Generación de Gráficos

figure
plot(tc,Concepcion);
figure
plot(ts,Santiago);
figure
plot(tt, Talca);
figure
plot(tv, Valparaiso);
figure
plot(tp, Pica_cms2);
figure
plot(tn, Northridge);

%% Calcular PGA
PGAc = max(Concepcion);
PGAs = max(Santiago);
PGAt = max(Talca);
PGAv = max(Valparaiso);
PGAp = max(Pica_cms2);
PGAn = max(northridge);

%% Maximos Sostenidos
% Para concepción

maximos_sostenidos = []; % Se genera una lista vacía
max = Concepcion(1,1); % Valor inicial de aceleración en el registro
j = 1;   %inicio en el tramo 1

for i = 2:mc
    if sign(Concepcion(i,1))*Concepcion(i-1,1)) == -1   %Si es negativo entonces hay un nuevo tramo
        % maximos_sostenidos = [maximos_sostenidos; max];
        maximos_sostenidos(j,1) = max;
        % if maximos_sostenidos(j,1) < max || isempty(maximos_sostenidos)
        %     maximos_sostenidos = [lista_maximos sostenidos; last_max];
        % end
        max = 0;
        p = j;
        j = p + 1;    % Siguiente tramo (no me funciona j = j + 1;)
        clear p % No es de interés guardar la variable p
    end
    if abs(Concepcion(i,1)) > abs(max)
        max = Concepcion(i,1);
    end
end

% Ordenar las aceleraciones sostenidas en orden descendente (la primera es la mayor aceleración)
sort(maximos_sostenidos,'descend');

%Mostrar las aceleraciones sostenidas
disp(maximos_sostenidos(3,1))
disp(maximos_sostenidos(5,1))
disp(maximos_sostenidos(7,1))
disp(maximos_sostenidos(9,1))

%% Duración del movimiento fuerte
g = 9.81; %m/s^2
umbral = 0.05*g;

ti = 0; % sec, tiempo inicial
tf = 0; % sec, tiempo final

for i = 2:mc
    if   
end

dur = tf-ti;