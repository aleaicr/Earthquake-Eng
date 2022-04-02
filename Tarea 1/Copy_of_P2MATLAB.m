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
tc = (0:dt:(mc-1)*dt).'; %Vector de tiempos para Concepción
ts = (0:dt:(ms-1)*dt).'; % Vector de tiempos para Santiago
tt = (0:dt:(mt-1)*dt).'; % Vector de tiempos para TAlca
tv = (0:dt:(mv-1)*dt).'; % Vector de tiempos para valparaíso
tp = (0:dt:(mp-1)*dt).'; % Vector de tiempos para Pica
tn = (0:dt:(mn-1)*dt).'; % Vector de tiempos apra Northridge - Tarzana

%% Ordenar Pica (registro viene con muchos NaN)
Pica2 = [];
for i = 1:mp
    if Pica(i,3) > 0
        Pica2 = Pica(i,1);
    end
    if Pica(i,2) < 0
        Pica2 = Pica(i,1);
    end
end
Pica = Pica2*981;  % g = 9.81 m/s2 = 981 cm/s2

%% Generación de Gráficos

% figure
% plot(tc,Concepcion)
% xlabel("Tiempo [s]")
% ylabel("Aceleración [cm/s^2]")
% figure
% plot(ts,Santiago);
% xlabel("Tiempo [s]")
% ylabel("Aceleración [cm/s^2]")
% figure
% plot(tt, Talca);
% xlabel("Tiempo [s]")
% ylabel("Aceleración [cm/s^2]")
% figure
% plot(tv, Valparaiso);
% xlabel("Tiempo [s]")
% ylabel("Aceleración [cm/s^2]")
% figure
% plot(tp, Pica);
% xlabel("Tiempo [s]")
% ylabel("Aceleración [cm/s^2]")
% figure
% plot(tn, Northridge);
% xlabel("Tiempo [s]")
% ylabel("Aceleración [cm/s^2]")

%% Calcular PGA
PGAc = max(Concepcion);
PGAs = max(Santiago);
PGAt = max(Talca);
PGAv = max(Valparaiso);
PGAp = max(Pica);
PGAn = max(Northridge);

%% Maximos Sostenidos
% Para concepción

maximos_sostenidos = []; % Se genera una lista vacía
max_array = [Concepcion(1,1)]; % Valor inicial de aceleración en el registro
j = 1;   %inicio en el tramo 1
signo = [];

for i = 2:mc
    if sign(Concepcion(i,1)*Concepcion(i-1,1)) == -1   %Si es negativo entonces hay un nuevo tramo
        % maximos_sostenidos = [maximos_sostenidos; max];
        maximos_sostenidos(j,1) = abs(max_array(j,1));
        maximos_sostenidos(j,2) = signo(j,1);
        % if maximos_sostenidos(j,1) < max || isempty(maximos_sostenidos)
        %     maximos_sostenidos = [lista_maximos sostenidos; last_max];
        % end
        p = j;
        j = p + 1;    % Siguiente tramo (no me funciona j = j + 1;)
        max_array(j,1) = 0;
        clear p % No es de interés guardar la variable p
    end
    if abs(Concepcion(i,1)) > abs(max_array(j,1))
        max_array(j,1) = Concepcion(i,1);
        signo(j,1) = sign(Concepcion(i,1));
    end
end

% Ordenar las aceleraciones sostenidas en orden descendente (la primera es la mayor aceleración)
sortrows(maximos_sostenidos,'descend')
for i = 1:j-1
    var = maximos_sostenidos(i,1);
    maximos_sostenidos(i,1) = var*maximos_sostenidos(i,2);
end
%Mostrar las aceleraciones sostenidas
disp(maximos_sostenidos(3,1))
disp(maximos_sostenidos(5,1))
disp(maximos_sostenidos(7,1))
disp(maximos_sostenidos(9,1))

%% Duración del movimiento fuerte
g = 9.81; %m/s^2
u = 0.05; % porcentaje del máximo admitido
umbral = u*g;

for i = 1:mc
    if abs(Concepcion(i,1)) > umbral
        ti = tc(i,1);      % guardamos la primera vez que supera ese umbral
        break %paramos para que no 
    end
end
for i = 1:mc
    if abs(Concepcion(i,1)) > umbral
        tf = tc(i,1);  % guardamos la última vez que se supera ese umbral
    end
end
dur = tf-ti;

%% Duración con Método 5% a 95% de Intensidad de Arias

% Generamos la "gráfica" de la Integral de Arias I_A

Ia = [];
Ia_new = 0; 
for i = 1:mc
    new_add = pi/(2*g)*abs(Concepcion(i,1))^2*dt;
    p = Ia_new;
    Ia_new = p + new_add;
    clear p
    Ia = [Ia; Ia_new];
end
clear Ia_new
% Máxima Ia (Para saber el porcentaje del total)

Iamax = max(Ia);

% Determinamos cuando se genera el 5 y 95 por ciento

for i = 1:mc
    if Ia(i,1) > 0.05*Iamax
        ti = tc(i,1);
    end
    if Ia(i,1) > 0.95*Iamax
        tf = tc(i,1);
    end
end

dur_mov_fuerte = tf - ti;