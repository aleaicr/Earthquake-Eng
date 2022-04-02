%% Inicializar
clear variables
close all
clc
disp('<a href="https://github.com/aleaicr/Tareas_Sismica">Github - Tareas_Sismica - Alexis Contreas R.</a>')
disp('Si no corre el registro de Pica es por que hay que quitar los comentarios')
disp('LINEA 24 Y 25')
%% Importar datos
% Todos los datos tienen paso temporal de 0.005
% Todos tienen unidad cm/s2 para la aceleración del registro excepto pica2005_ew.txt

Regi = readmatrix("pica2005_ew_arreglado.txt");
Nombre = "Pica";
mc = length(Regi);

%% Vector de tiempo;

dt = 0.005; % Paso temporal de la toma de datos
t_reg = (0:dt:(mc-1)*dt).'; % Vector que contiene los tiempos


%% Ordenar Pica (registro viene con muchos NaN)
% QUITAR COMENTARIO SI SE ESTÁ VIENDO EL REGISTRO DE PICA
% Regi2 = Regi;
% Regi = Regi2*981;  % g = 9.81 m/s2 = 981 cm/s2

%% Generación de Gráficos

figure
plot(t_reg,Regi)
xlabel("Tiempo [s]")
ylabel("Aceleración [cm/s^2]")
title(Nombre)
grid on


%% Calcular PGA
PGA_reg = max(Regi);

%% Maximos Sostenidos

maximos_sostenidos = []; % Se genera una lista vacía
max_array = [Regi(1,1)]; % Valor inicial de aceleración en el registro
j = 1;   % Inicio en el tramo 1
signo = [];

for i = 2:mc
    if sign(Regi(i,1)*Regi(i-1,1)) == -1   %Si es negativo entonces hay un nuevo tramo
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
    if abs(Regi(i,1)) > abs(max_array(j,1))
        max_array(j,1) = Regi(i,1);
        signo(j,1) = sign(Regi(i,1));
    end
end

% Ordenar las aceleraciones sostenidas en orden descendente (la primera es la mayor aceleración)

max_ordenados = sortrows(maximos_sostenidos,'descend');

for i = 1:j-1
    maximos_sostenidos_ordenados(i,1) = max_ordenados(i,1)*maximos_sostenidos(i,2);
end

%Mostrar las aceleraciones sostenidas
sprintf('La 3ra aceleracion máxima sostenida es %d cm/s2',maximos_sostenidos_ordenados(3,1))
sprintf('La 5ta aceleracion máxima sostenida es %d cm/s2',maximos_sostenidos_ordenados(5,1))
sprintf('La 7ma aceleracion máxima sostenida es %d cm/s2',maximos_sostenidos_ordenados(7,1))
sprintf('La 9na aceleracion máxima sostenida es %d cm/s2',maximos_sostenidos_ordenados(9,1))

%% Duración del movimiento fuerte
g = 9.81; %m/s2
u = 0.05; % porcentaje del máximo admitido
umbral = u*g*100;

for i = 1:mc
    if abs(Regi(i,1)) > umbral
        ti = t_reg(i,1);      % guardamos la primera vez que supera ese umbral
        break % paramos para que no siga buscando 
    end
end
for i = 1:mc
    if abs(Regi(i,1)) > umbral
        tf = t_reg(i,1);  % Guardamos la última vez que se supera ese umbral
    end
end
dur = tf - ti; % Duración del movimiento fuerte
sprintf('La duración del movimiento fuerte es de %d [sec]',dur)
%% Duración con Método 5% a 95% de Intensidad de Arias

% Generamos la "gráfica" de la Integral de Arias I_A

Ia(1,1) = 0;

for i = 2:mc+1
    new_add = pi/(2*g*100)*Regi(i-1,1)^2*dt; % Integrando
    Ia(i,1) = Ia(i-1,1) + new_add; %Integral de arias (se puede graficar)
end

clear new_add p

% Máximo valor de la Integral de Arias

Iamax = max(Ia); % Para saber que valor es 5 y 95 % del total

% Determinamos cuando se genera el 5 y 95 por ciento

for i = 1:mc
    if Ia(i,1) > 0.05*Iamax
        ti = t_reg(i,1); %Encontramos el primer
        break % Paramos de buscar
    end
end
for i = 1:mc
    if Ia(i,1) > 0.95*Iamax
        tf = t_reg(i,1); %Encontramos hasta el último, hasta finalizar el registro
    end
end

dur_arias = tf - ti; % Duración con Método de la Intensidad de Arias

sprintf('La duración con el método de la Intensidad de Arias es %d [sec]',dur_arias)
clear ans