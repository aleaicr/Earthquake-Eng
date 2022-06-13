%% Tarea 4 Pregunta 2
% Ingeniería Sísmica
% Alexis Contreras R. - Cristobal Adasme 

% Para los 4 registros del terremoto de Northridge 1994 (Burbank, Los
% Ángeles, North Hollywood y Sylmar) calcular, para cada registro (4):

% I) El mínimo coeficiente sísmico Cy que garantiza que un sistema bilineal
% de 1GDL con T = 0.5s y xi = 5% no colapse.

% II) El máximo factor de reducción de resistencia que evita la
% inestabilidad dinámica Rc.

% III) Comparar resultados para Rc con los que se obtienen utilizando la
% fórmula empírica desarrollada por Akkar y Miranda (2003), vista en clases

% En todos los puntos compara y comentar resultados


%% Inicializar
clear variables
close all
clc
%% Parámetros
% Parámetros
alpha = -0.15;                                                               % Pendiente de degradación modelo bilineal
xi = 0.05;                                                                  % Razón de amortiguamiento de la estructura
Tn = 0.5;                                                                     % Periodo de la estructura
R = 10^(-10);                                                               % Criterio de parada
beta = 1/6;                                                                 % Beta (Método de Newmark)

% Condiciones iniciales
u_i = 0;    
ud_i = 0;

% Condición de colapso
u_colapso = 100;                                                            

%% Importando Registros
% BURBANK
Registro_Burbank = importdata('Northridge_BURBANK_EW.txt');
dt_B = 0.020;
T_B = transpose(0:dt_B:(size(Registro_Burbank.data(:,2),1)-1)*dt_B);                     % 

% LOSANGELESCA
Registro_LosAngelesCA = importdata('Northridge_LOSANGELESCA_EW.txt');
dt_LCA = 0.010;
T_LCA = transpose(0:dt_LCA:(size(Registro_LosAngelesCA.data(:,2),1)-1)*dt_LCA);

% NORTHHOLLYWOOD
Registro_NorthHollywood = importdata('Northridge_NORTHHOLLYWOOD_EW.txt');
dt_NH = 0.020;
T_NH = transpose(0:dt_NH:(size(Registro_NorthHollywood.data(:,2),1)-1)*dt_NH);

% SYLMAR
Registro_Sylmar = importdata('Northridge_SYLMAR_EW.txt');
dt_S = 0.020;
T_S = transpose(0:dt_S:(size(Registro_Sylmar.data(:,2),1)-1)*dt_S);

figure
subplot(4,1,1)
plot(Registro_Burbank.data(:,1),Registro_Burbank.data(:,2))
title('Registros de aceleración')
ylabel('Burbank [cm/s^2]')
subplot(4,1,2)
plot(Registro_LosAngelesCA.data(:,1),Registro_LosAngelesCA.data(:,2))
ylabel('Los Ángeles California [cm/s^2]')
subplot(4,1,3)
plot(Registro_NorthHollywood.data(:,1),Registro_NorthHollywood.data(:,2))
ylabel('North Hollywood [cm/s^2]')
subplot(4,1,4)
plot(Registro_Sylmar.data(:,1),Registro_Sylmar.data(:,2))
ylabel('Sylmar [cm/s^2]')
xlabel('Tiempo (t) [s]')

%% Parte (I)
Cy_barrido = 0:0.0001:1;                                                            % Barrido de Cy
n = length(Cy_barrido);
% Determinando Cy_minimo para que el sistema con comportamiento biliean no
% colapse

% Burbank
for i = 1:n
    [u_B,ud_B,udd_B,fs,fluencia,ductilidad] = Newmark_No_Lineal_Bilineal(alpha,beta,Tn,xi,dt_B,u_i,ud_i,Registro_Burbank.data(:,2),Cy_barrido(1,n-i+1)*980,R,T_B);
    if max(abs(u_B)) > u_colapso
        Cy_B = Cy_barrido(1,n-i+2);
        break
    end
end

% Los Ángeles CA
for i = 1:n
    [u_LCA,ud_LCA,udd_LCA,fs,fluencia,ductilidad] = Newmark_No_Lineal_Bilineal(alpha,beta,Tn,xi,dt_LCA,u_i,ud_i,Registro_LosAngelesCA.data(:,2),Cy_barrido(1,n-i+1)*980,R,T_LCA);
    if max(abs(u_LCA)) > u_colapso
        Cy_LCA = Cy_barrido(1,n-i+2);
        break
    end
end
% North Hollywood
for i = 1:n
    [u_NH,ud_NH,udd_NH,fs,fluencia,ductilidad] = Newmark_No_Lineal_Bilineal(alpha,beta,Tn,xi,dt_NH,u_i,ud_i,Registro_NorthHollywood.data(:,2),Cy_barrido(1,n-i+1)*980,R,T_NH);
    if max(abs(u_NH)) > u_colapso
        Cy_NH = Cy_barrido(1,n-i+2);
        break
    end
end
% Sylmar
for i = 1:n
    [u_S,ud_S,udd_S,fs,fluencia,ductilidad] = Newmark_No_Lineal_Bilineal(alpha,beta,Tn,xi,dt_S,u_i,ud_i,Registro_Sylmar.data(:,2),Cy_barrido(1,n-i+1)*980,R,T_S);
    if max(abs(u_S)) > u_colapso
        Cy_S = Cy_barrido(1,n-i+2);
        break
    end
end

fprintf('Los mínimos coeficientes sísmicos que garantiza que el sistema no lineal\n')
fprintf('con comportamiento bilineal no colapse son: \n')
tabla = table();
tabla.Registros = ["Burbank"; "Los Angeles CA"; "North Hollywood"; "Sylmar"];
tabla.Cy = [Cy_B; Cy_LCA; Cy_NH; Cy_S];
disp(tabla)
clear tabla
%% Parte (II)
% Determinación del máximo  factor de reducción de resistencia que
% evita la inestabilidad dinámica

% Rc = Ce/Cy_min , Cy_min es el Cy menor que no produce colapso
% Determinando Ce -> Newmark_lineak -> Ce = Sa/g
[Sd,Sv,Sa_B,PSv,PSa] = Newmark_Lineal(beta,xi,dt_B,u_i,ud_i,Registro_Burbank.data(:,2));
[Sd,Sv,Sa_LCA,PSv,PSa] = Newmark_Lineal(beta,xi,dt_LCA,u_i,ud_i,Registro_LosAngelesCA.data(:,2));
[Sd,Sv,Sa_NH,PSv,PSa] = Newmark_Lineal(beta,xi,dt_NH,u_i,ud_i,Registro_NorthHollywood.data(:,2));
[Sd,Sv,Sa_S,PSv,PSa] = Newmark_Lineal(beta,xi,dt_S,u_i,ud_i,Registro_Sylmar.data(:,2));


% Buscamos el Sa(Tn), recordando que nuestro Tn = 1 sec
for i = 1:size(Sa_B,1)
    if Sa_B(i,1) == Tn
        Sa_Burbank = Sa_B(i,2);
    end
end
for i = 1:size(Sa_LCA,1)
    if Sa_LCA(i,1) == Tn
        Sa_LosAngelesCA = Sa_LCA(i,2);
    end
end
for i = 1:size(Sa_NH,1)
    if Sa_NH(i,1) == Tn
        Sa_NorthHollywood = Sa_NH(i,2);
    end
end
for i = 1:size(Sa_S,1)
    if Sa_S(i,1) == Tn
        Sa_Sylmar = Sa_S(i,2);
    end
end
disp('Calculamos las respuestas según el espectro de respuesta para aceleraciones Tn [s] Sa [cm/s2]')
tabla = table();
tabla.Registros = ["Burbank"; "Los Angeles CA"; "North Hollywood"; "Sylmar"];
tabla.Tn = [Tn; Tn; Tn; Tn];
tabla.Sa = [Sa_Burbank; Sa_LosAngelesCA; Sa_NorthHollywood; Sa_Sylmar];
disp(tabla)
clear tabla

% Luego Rc = Ce/Cy_min
Ce_B = Sa_Burbank/980;
Rc_B = Ce_B/Cy_B;
Ce_LCA = Sa_LosAngelesCA/980;
Rc_LCA = Ce_LCA/Cy_LCA;
Ce_NH = Sa_NorthHollywood/980;
Rc_NH = Ce_NH/Cy_NH;
Ce_S = Sa_Sylmar/980;
Rc_S = Ce_S/Cy_S;

disp('Calculamos Ce y Rc para cada registro')
tabla = table();
tabla.Registros = ["Burbank"; "Los Angeles CA"; "North Hollywood"; "Sylmar"];
tabla.Ce = [Ce_B;Ce_LCA;Ce_NH;Ce_S];
tabla.Rc = [Rc_B;Rc_LCA;Rc_NH;Rc_S];
disp(tabla)
clear tabla


%% Parte (III)
% Akkar & Miranda (2003)
% Rc = 1 + a*(-alpha)^b
% a y b dependen del periodo de la estructura

a = 0.26*(1-exp(-7.5*Tn));
b = 0.89 + 0.04*Tn + 0.15*log(Tn);
Rc_AM = 1+a*(-alpha)^(-b);

fprintf('Aplicando Akker & Miranda (2003) para determinar R de colapso \n')
fprintf('a = %f; b = %f\n',a,b)
fprintf('Rc = %f\n',Rc_AM)
fprintf('La diferencia porcentual entre métodos\n')
tabla = table();
tabla.Rc = ["Ce/Cy"; "Akkar-Miranda (2003)"; "Ce/Cy / AM(2003)"];
tabla.Burbank = [Rc_B;Rc_AM;Rc_B/Rc_AM];
tabla.LosAngelesCA = [Rc_LCA;Rc_AM;Rc_LCA/Rc_AM];
tabla.NorthHollywood = [Rc_NH;Rc_AM;Rc_NH/Rc_AM];
tabla.Sylmar = [Rc_S;Rc_AM;Rc_S/Rc_AM];
disp(tabla)
clear tabla
fitlm([Rc_B;Rc_LCA;Rc_NH;Rc_S],[Rc_AM;Rc_AM;Rc_AM;Rc_AM])


alpha_plot = (-0.01:-0.01:-1).';
alpha_length = length(alpha_plot);
Rc_am_plot = zeros(alpha_length,1);
period = (0.01:0.01:1).';
period_length = length(period);

for j  = 1:period_length
    a = 0.26*(1-exp(-7.5*period(j,1)));
    b = 0.89 + 0.04*period(j,1) + 0.15*log(period(j,1));
    for i = 1:alpha_length
        Rc_am_plot(j,i) = 1+a*(-alpha_plot(i,1))^b;
    end
end

figure
surf(alpha_plot,period,Rc_am_plot)
xlabel('\alpha')
ylabel('Tn [sec]')
zlabel('Rc Akkar & Miranda (2003)')
title('a = 0.26*(1-exp(-7.5*Tn)); b = 0.89 + 0.04*Tn + 0.15*log(Tn); Rc_{AM} = 1+a*(-alpha)^(-b);')