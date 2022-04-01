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
tc = 0:dt:(mc-1)*dt;
ts = 0:dt:(ms-1)*dt;
tt = 0:dt:(mt-1)*dt;
tv = 0:dt:(mv-1)*dt;
tp = 0:dt:(mp-1)*dt;
tn = 0:dt:(mn-1)*dt;

%% Conversión a cm/s2 para Pica

Pica_cms2 = Pica*981;  % g = 9.81 m/s2 = 981 cm/s2

%% Generación de Gráficos

figure
plot(tc,Concepcion);
figure
plot(Santiago,ts);
figure
plot(Talca,ts);
figure
plot(Valparaiso,tv);
figure
plot(Pica_cms2,tp);
figure
plot(Northridge,tn);

%% Calcular PGA
PGAc = max(Concepcion);
PGAs = max(Santiago);
PGAt = max(Talca);
PGAv = max(Valparaiso);
PGAp = max(Pica_cms2);
PGAn = max(northridge);

%% Maximos Sostenidos
elem_inicial = 0;
% Anotar los máximos sostenidos

lista_maximos_sostenidos = [];
max = Concepcion(1,1); %Valor inicial
j = 0;   %inicio en el tramo 0

for i = 2:mc
    if sign(Concepcion(i,1))*Concepcion(i-1,1)) != 1   %Si es negativo entonces hay un nuevo tramo
        signo_nuevo = sign(Concepcion(i,1)*Concepcion(i-1,1));
        max = new_max; 
        j = j + 1;    % Siguiente tramo
    end
    if Concepcion(i,1) > max
        new_max = Concepcion(i,1);
        if lista_maximos_sostenidos(j,1) < new_max
            lista_maximos_sostenidos = [lista_maximos sostenidos; new_max];
        end
end
