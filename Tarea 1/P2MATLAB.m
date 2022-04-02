%% Inicializar
clear variables
close all
clc

%% Importar datos
% Todos los datos tienen paso temporal de 0.005
% Todos tienen unidad cm/s2 para la aceleración del registro excepto pica2005_ew.txt

Reg = readmatrix("Concepcion2010-L.txt");
Nombre = "Concepción"
mc = length(Reg);

%% Vector de tiempo;

dt = 0.005; % Paso temporal de la toma de datos
t_reg = (0:dt:(mc-1)*dt).'; % Vector que contiene los tiempos


%% Ordenar Pica (registro viene con muchos NaN)
% QUITAR COMENTARIO SI SE ESTÁ VIENDO EL REGISTRO DE PICA
% Pica2 = [];
% for i = 1:mp
%     if Pica(i,3) > 0
%         Pica2 = Pica(i,1);
%     end
%     if Pica(i,2) < 0
%         Pica2 = Pica(i,1);
%     end
% end
% Pica = Pica2*981;  % g = 9.81 m/s2 = 981 cm/s2

%% Generación de Gráficos

figure
plot(t_reg,reg)
xlabel("Tiempo [s]")
ylabel("Aceleración [cm/s^2]")


%% Calcular PGA
PGA_reg = max(Reg);

%% Maximos Sostenidos

maximos_sostenidos = []; % Se genera una lista vacía
max_array = [Reg(1,1)]; % Valor inicial de aceleración en el registro
j = 1;   % Inicio en el tramo 1
signo = [];

for i = 2:mc
    if sign(Reg(i,1)*Reg(i-1,1)) == -1   %Si es negativo entonces hay un nuevo tramo
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
    if abs(Reg(i,1)) > abs(max_array(j,1))
        max_array(j,1) = Reg(i,1);
        signo(j,1) = sign(Reg(i,1));
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
    if abs(Reg(i,1)) > umbral
        ti = t_reg(i,1);      % guardamos la primera vez que supera ese umbral
        break %paramos para que no 
    end
end
for i = 1:mc
    if abs(Reg(i,1)) > umbral
        tf = t_reg(i,1);  % guardamos la última vez que se supera ese umbral
    end
end
dur = tf - ti; % Duración del movimiento fuerte

%% Duración con Método 5% a 95% de Intensidad de Arias

% Generamos la "gráfica" de la Integral de Arias I_A

Ia = [];
Ia_new = 0; 

for i = 1:mc
    new_add = pi/(2*g)*abs(Reg(i,1))^2*dt;
    p = Ia_new;
    Ia_new = p + new_add;
    clear p
    Ia = [Ia; Ia_new];
end
clear Ia_new

% Máximo valor de la Integral de Arias

Iamax = max(Ia); % Para saber que valor es 5 y 95 % del total

% Determinamos cuando se genera el 5 y 95 por ciento

for i = 1:mc
    if Ia(i,1) > 0.05*Iamax
        ti = t_reg(i,1);
    end
    if Ia(i,1) > 0.95*Iamax
        tf = t_reg(i,1);
    end
end

dur_arias = tf - ti; % Duración con Método de La Integral de Arias