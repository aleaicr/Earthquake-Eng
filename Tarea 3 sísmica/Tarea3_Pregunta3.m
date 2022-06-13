%% Tarea 3 - Pregunta 3 - Ingeniería Sísmica
% Cristóbal Adasme |  Alexis Contreras


%% Inicializar
clear variables
close all
clc

%% Datos y parámetros
% Se asume comportamiento elastoplástico perfecto
Registro = importdata('stgocentro2010-L.txt', ' ', 1);
uppg = Registro.data/100;
dt = 0.005;                                                                 % Paso temporal del sampling del registro uppg
xi = 0.05;                                                                  % Razón de amortiguamiento de 5%
T = [0.15; 0.75; 2.5];                                                      % Periodos para las estructuras de 3,10 y 30 pisos
N = length(T);
W = 9.81;                                                                   % W = m*g, m = 1 => W = g
beta = 1/4;
ui = 0;                                                                     % Condición inicial de desplazamiento   % Se parte del reposo
udi = 0;                                                                    % Condición inicial de velocidad
Error = 10^(-5);                                                            % Error admisible
t = transpose(0:dt:(length(uppg)-1)*dt);

%% Parte (I)
[Sd,Sv,Sa,PSv,PSa] = Newmark_Lineal(beta,xi,dt,ui,udi,uppg);

Sa(:,2) = Sa(:,2)/(9.81);
Ce = zeros(N,1);
for i = 1:N
    for j = 1:size(Sa,1)
        if T(i,1) == Sa(j,1)
            Ce(i,1) = Sa(j,2);
        end
    end
end
disp('Los resultados de Ce son los siguientes para los period T')
table(T,Ce)

%% Parte (II)
R = 6;                                                                      % Factor de Reducción de la Respuesta R = 6 para las tres estructuras
Cy = Ce/R;
Fy = Cy*W;
k = zeros(N,1);
uy = zeros(N,1);
mu = zeros(N,1);
for i = 1:N
    [u,ud,udd,fs,fluencia,ductilidad] = Newmark_No_Lineal(beta,T(i,1),xi,dt,ui,udi,uppg,Fy(i,1),Error,t);
    % Graficamos la historia de desplazamientos relativos
    figure
    plot(t,u)
    xlabel('Tiempo [s]')
    ylabel('Desplazamiento [m]')
    % Graficamos el comportamiento histerético fs vs u
    figure
    plot(u,fs)
    xlabel('Desplazamiento [m]')
    ylabel('Fuerza Restitutiva')
    
    % Calcular la demanda de desplazamiento mu
    k(i,1) = 4*pi^2/T(i,1)^2;
    uy(i,1) = Fy(i,1)/k(i,1);
    mu(i,1) = ductilidad;
    fprintf('Demanda de ductilidad para la estructura %d \n', i)
end

table(T,mu)
disp('Se generaron los gráficos de historia de desplazamiento y del comportamiento histerético')

disp('Tabla T, Cy_II')
fprintf('R = %d\n',6)
table(T,Cy)
Cy_ii = Cy;
%% Parte (III)
mu_max = 4;                                                                 % Demanda de ductilidad debe ser menor al mu_max = 4
Cy_min = zeros(N,1);
figure
hold on
for i = 1:N
    Cy_min(i,1) = Cy_vs_mu(beta,T(i,1),xi,dt,ui,udi,uppg,Error,t,mu_max);   % Grafica, pero retorna los Cy_minimos para un mu_max = 4
%     k = 4*pi^2/T(i,1)^2;
%     Cy = transpose(0.005:0.001:0.5);
%     
%     for j=1:size(Cy)
%         Fy = Cy(j,1)*9.81;
%         uy = Fy/k;
%         [u,ud,udd,fs,fluencia] = Newmark_No_Lineal(beta,T(i,1),xi,dt,ui,udi,uppg,Fy,Error,t);
%         mu(j,1) = max(abs(u))/uy;
%     end
%     
%     plot(mu,Cy)
%     xlabel('Demanda de Ductilidad')
%     ylabel('Coeficiente Sismico Inelastico')
%     mu = flip(mu,1);
%     Cy = flip(Cy,1);
%     c = 0;
%     % Determinar Cy_min
%     for j = 1:size(mu)
%         if c == 1
%             c = 1;
%         elseif mu(j,1) > mu_max
%             Cy_min = Cy(j-1,1);
%             c = 1;
%         elseif mu(i,1) == mu_max
%             Cy_min = Cy(j,1);
%             c = 1;
%         end
%     end
end
xlabel('Demanda de Ductilidad (\mu)')
ylabel('Coeficiente Sismico Inelastico (C_y)')
legend('T=0.15[s]','T=0.75 [s]','T=2.5[s]')
hold off

disp('Se grafica el coeficiente sísmico inelástico vs demanda de ductilidad ')
disp('Se aprecia que para las estructuras de periodo menor, el coeficiente sísmico disminuye a una menor tasa mientras aumenta la demanda de dictuilidad')
disp('Cy aumenta mucho para bajas demandas de ductilidad')

%% Parte (IV)
% Calcular la demanda sísmica inelástica para que la demanda de ductilidad
% de desplazamientos mu sea mayor que 4
table(T,Cy_min)
Cy_iv = Cy_min;

% Para la última pregunta
disp_iv = zeros(N,1);
for i = 1:N
    [u,ud,udd,fs,fluencia] = Newmark_No_Lineal(beta,T(i,1),xi,dt,ui,udi,uppg,Cy_min(i,1)*W,Error,t);
    disp_iv(i,1) = max(abs(u));
end


%% Parte (V)
Rmu = zeros(N,1);
Phi = zeros(N,1);
Cy = zeros(N,1);

for i = 1:N
    [Rmu(i,1),Phi(i,1)] = R_Miranda_Aluvial(T(i,1),mu_max);                % Factor de Reducción de la Respuesta propuestos por Miranda (1993)
    % Utilizando la Demanda sísmica elástica calculada Ce en Parte (I) 
    Cy(i,1) = Ce(i,1)/Rmu(i,1);                                            % Coeficiente de Demanda sísmica inelástica con R de Miranda (1993)
end

Fy = Cy*W;                                                                 % Fy con Cy con R de Miranda (1993) para cada estructura
disp('A continuación se muestra lo solicitado')
table(T,Phi,Rmu,Cy,Fy)
Cy_v = Cy;

%% Parte (VI)
I = 1.2;                                                                    % Categoria III: Universidad
A0 = 0.3*9.81;                                                              % Zona Sísmica 2: Santiago
S = 1.05;                                                                   % Suelo Tipo C
T0 = 0.4;
T_p = 0.45;
n = 1.4;
p = 1.6;
Ro = 11;                                                                    % Hormigón Armado
R = 7;                                                                      % Hormigón Armado
Sa = zeros(N,1);
R_ast = zeros(N,1);
alfa = zeros(N,1);
Ve = zeros(N,1);
Vd = zeros(N,1);
disp_vi = zeros(N,1);
Cy_vi = zeros(N,1);

for i = 1:N
    R_ast(i,1) = 1 + T(i,1)/(0.1*T0+T(i,1)/Ro);
    alfa(i,1) = (1+4.5*(T(i,1)/T0)^p)/(1+(T(i,1)/T0)^3);
    Sa(i,1) = S*A0*alfa(i,1)/(R_ast(i,1)/I);
    Ve(i,1) = Sa(i,1)/9.81*W;
    Vd(i,1) = I/R_ast(i,1)*Sa(i,1)/9.81*W;
    [u,ud,udd,fs,fluencia] = Newmark_No_Lineal(beta,T(i,1),xi,dt,ui,udi,uppg,Vd(i,1),Error,t);
    disp_vi(i,1) = max(abs(u));
    Cy_vi(i,1) = 2.75*S*A0/(9.81*R)*(T_p/T(i,1))^n;
end

disp('Factor de amplificación alpha')
table(T,alfa,R_ast)
disp("Demanda Sísmica Elástica, Inelástica, desplazamientos asociados, coeficiente sísmico")
table(T,Ve,Vd,disp_vi,Cy_vi)


%% P7
% Realizar un resumen de los métodos
table(T,Cy_ii,Cy_iv,Cy_v,Cy_vi)
table(T,disp_iv,disp_vi)

close all


