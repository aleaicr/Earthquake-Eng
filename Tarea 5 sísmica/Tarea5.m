%% Tarea 5 Ingeniería Sísmica 2022-1
% Cristóbal Adasme - Alexis Contreras

%% Inicializar
clear variables
close all
clc

%% Pregunta 1

% Parámetros
g = 980; % cm/s2
htotal = 40; % m
Wtotal = 9000; %tonf
cant_pisos = 10;                                                            % Cantidad de pisos
h = htotal/cant_pisos;                                                      % Altura de cada piso
W = Wtotal/cant_pisos; % tonf                                               % Peso W
k = 500; %tonf/cm
Pisos = (cant_pisos:-1:1)';
Modos = flip(Pisos);

k1 = 0.2*k;
k2 = 0.2*k;
k3 = 0.4*k;
k4 = 0.4*k;
k5 = 0.6*k;
k6 = 0.6*k;
k7 = 0.8*k;
k8 = 0.8*k;
k9 = k;
k10 = k;

m1 = 0.5*W/g;
m2 = 0.8*W/g;
m3 = W/g;
m4 = W/g;
m5 = W/g;
m6 = W/g;
m7 = W/g;
m8 = W/g;
m9 = W/g;
m10 = 1.4*W/g;

M = diag([m1 m2 m3 m4 m5 m6 m7 m8 m9 m10]);
K = [k1 -k1 0 0 0 0 0 0 0 0;
    -k1 k2+k1 -k2 0 0 0 0 0 0 0;
    0 -k2 k3+k2 -k3 0 0 0 0 0 0;
    0 0 -k3 k4+k3 -k4 0 0 0 0 0;
    0 0 0 -k4 k5+k4 -k5 0 0 0 0;
    0 0 0 0 -k5 k6+k5 -k6 0 0 0;
    0 0 0 0 0 -k6 k7+k6 -k7 0 0;
    0 0 0 0 0 0 -k7 k8+k7 -k8 0;
    0 0 0 0 0 0 0 -k8 k9+k8 -k9;
    0 0 0 0 0 0 0 0 -k9 k10+k9];
G = [-m1; -m2; -m3; -m4; -m5; -m6; -m7; -m8; -m9; -m10];
r = diag(ones(cant_pisos));
disp('Matriz de masa [M]')
disp(M)
disp('Matriz de rigidez [K]')
disp(K)

%% P1 a)
[Phi, lambda] = eig(K,M);                                                   % Problema de valores y vectores propios
wn = diag(lambda.^0.5);                                                             % Frecuencia de cada modo
Tn = 2*pi./wn;                                                              % Periodo de cada modo

tabla = table();
tabla.Modo_n = Modos;
tabla.Periodos_Tn = Tn;
tabla.FreqAngular_wn = wn;
disp(tabla)
clear tabla

%% P1 b)
% Caluclar factores de participación modal Gamma_n y las formas modales {Phi_n}
% Graficar Gamma_n{Phi_n} para los primeros 5 primeros modos
Mn = diag(Phi.'*M*Phi);
Kn = diag(Phi.'*K*Phi);

% Gamma_n
Gamma = Phi.'*M*r./Mn;

% Comentarios
tabla = table(Modos,Gamma);
disp(tabla)
clear tabla

% Graficamos Gamma_n*phi_n para 5 modos
cant_modos_dibujar = 5;                                                     % Cantidad de modos que queremos dibujar (cambiar si se quieren dibujar más)
figure
for i = 1:cant_modos_dibujar
    subplot(1,cant_modos_dibujar,i)
    plot([Gamma(i)*Phi(:,i); 0],[Pisos; 0], '-o')
    xlabel(['\Gamma_', num2str(Modos(i)), '\{\phi_{', num2str(Modos(i)), '}\}'])
    xlim([-1.5 1.5])
    if i == 1
        ylabel('Piso')
        hold on
    end
end
hold off
sgtitle('\Gamma_n \{\phi_n\}')

% Todos juntos
figure
hold on
for i = 1:cant_modos_dibujar
    plot([Gamma(i)*Phi(:,i); 0],[Pisos; 0], '-o')
end
hold off
xlabel('\Gamma_n \{\Phi_n\}')
ylabel('Piso')
legend(strcat('Modo ',string(1:cant_modos_dibujar)))
% legend('Modo 1', 'Modo 2', 'Modo 3', 'Modo 4', 'Modo 5')
grid on

%% P1 c)
% Calcular Masa modal Equivalente
% y determine el número de modos necesario para que la suma de la 
% masa modal equivalente sea al menos un 90% de la masa total del edificio

Mn_ast = zeros(cant_pisos,1);
for n = 1:cant_pisos
    for j = 1:cant_pisos
        Mn_ast(n,1) = Mn_ast(n,1) + M(j,j)*Gamma(n)*Phi(j,n);
    end
end

PMass = Mn_ast/(Wtotal/g);                                                  % Porcentaje de Masa modal Equivalente de cada modo
PMass_acum = zeros(cant_pisos,1);                                           % Porcentaje de Masa modal equivalente acumulado
PMass_acum(1,1) = PMass(1,1);                                               
ifprinted = 0;                                                                % Contador para ver si ya encontró la solución
for i = 2:cant_pisos
    PMass_acum(i,1) = PMass_acum(i-1,1) + PMass(i,1);
    if PMass_acum(i,1) > 0.9 && ifprinted == 0
        fprintf('Solo se requieren los primeros %i modos \n \n', Modos(i,1))
        ifprinted = 1;
    end
end

% Comentarios consola
tabla = table();
tabla.Modos = Modos;
tabla.MasaModalEquivalente = Mn_ast;
tabla.PorcntajeDeMasa = PMass;
tabla.PorcntajeDeMasaAcum = PMass_acum;
disp(tabla)
clear tabla
% Para verificar
% Sum(Mn_Ast) = Sum(diag(M))

%% P1 d)
% Amortiguamiento del primer modo utilizando fórmula propuesta por Cruz & Miranda (2021)
xi = zeros(cant_pisos,1);
% Para el primer modo fórmula Cruz & Miranda (2021)
xi(1,1) = 0.28*htotal^(-0.52);

%% P1 e)
% Amortiguamiento para los modos superiores utilizando fórmula propuesta por Cruz & Miranda (2017)
for i = 2:cant_pisos
    xi(i,1) = xi(1,1)*(1+0.12*(wn(i,1)/wn(1,1)-1));
end

%% P1 f)
% Graficar xi_n vs wn
figure
plot(wn,xi,'-o')
xlabel('\omega_n')
ylabel('\xi_n')
title('\xi_n vs \omega_n')
grid on

% Da una recta como es de esperar viendo la ecuación de xi para modos superiores

% Figura para ver como varía Xi para cada modo (no es una recta)
% figure
% plot(Modos,xi)
% xlabel('Modos')
% ylabel('\xi_n')

%% Pregunta 2
% Para el registro de Santiago Centro del terremoto Maule 2010
Registro = importdata('stgocentro2010-L.txt');

%% P2 a)
% Calcular pseudo-espectro de aceleraciones para el primer modo
% Parámetros
modoi = 1;
Tn_ini = 0.01;
Tn_step = 0.01;
Tn_fin = 6;
udi = 0;
ui = 0;
beta = 1/6;
dt = 0.02;
% Buscar PSa
T_xlabel = Tn_ini:Tn_step:Tn_fin;
[Sd,Sv,Sa,PSv,PSa] = Newmark_Lineal(beta,xi(modoi),0.005,ui,udi,Registro.data,T_xlabel,Mn(1)); % cm y segundos
% PSa [cm/s2] ... Psa [g] => PSa *980[cm/s2] = PSa [g]
% Figura PSa para cada Tn
figure
plot(T_xlabel,PSa/g)
xlabel('Periodo T (sec)')
ylabel('PSa_1 [g]')
title('Pseudo-Espectro de Aceleraciones Modo 1')
grid on

%% P2 b)
% Calcular coeficiente sísmico asociado a demanda elástica del primer modo

% Encontrando periodo
n = round(Tn(1)/Tn_step);

% Determinación de Ce1
Ce1_Sa = Sa(n)/g;
Ce1_PSa = PSa(n)/g;

%Comentario consola
fprintf('El valor del coeficiente sísmico elástico Ce_1 = %f \n \n', Ce1_PSa)

%% P2 c)
% Suponiendo coeficiente de Importancia I = 1, calcule el corte basal
% asociado al primer modo en tonf
% Qo = C*I*P = Ce1*1*Wtotal 
Qo1 = Ce1_PSa*Wtotal;    % tonf

fprintf('Corte Basal del primer modo Qo_1 = C*I*P = Ce_1*Wtotal = %.4f [tonf] \n \n',Qo1)

%% P2 d)
% Utilizando la distribución en altura de la norma NCh433, calcule las
% cargas laterales equivalentes en cada piso, en tonf

% Ak    
Ak = zeros(cant_pisos,1);                                                   % Forma de distribución de cargas
Zk = (h:h:htotal).';                                                        % Altura de cada piso
Ak(1) = sqrt(1)-sqrt(1-Zk(1)/htotal);                                       % Forma de distribución de cargas del primer piso

for i = 2:cant_pisos
    Ak(i,1) = sqrt(1-Zk(i-1)/htotal)-sqrt(1-Zk(i)/htotal);
end
Ak = flip(Ak);                                                              % Damos vuelta para que quede con pisos de arriba a abajo

% Fk
Pk = diag(M)*g;                                                             % Peso de cada piso (de arriba a abajo)
sumAjPj = sum(Ak.*Pk);                                                      % Suma_(j=1)^N (Aj*Pj)
Fk1 = zeros(cant_pisos,1);
for i = 1:cant_pisos
    Fk1(i,1) = Ak(i,1).*Pk(i,1)/sumAjPj*Qo1;                                % Fuerza de cada piso (de arriba a abajo)
end

figure
plot(Fk1,Pisos,'-o')
xlabel('Fk_1 [tonf]')
ylabel('Piso')
grid on

%% P2 e)
% Para cada piso calcular:
% i) Esfuerzos de corte (en tonf)
% ii) Desplazamientos laterales (en cm)
% iii) Razón de derivas de piso (Expresadas como porcentaje)

% P2 e) i) Esfuerzos de corte (en tonf)

corte_pisos1 = zeros(cant_pisos,1);
for i = 1:cant_pisos
    corte_pisos1(i,1) = sum(Fk1(1:i));    %tonf
end

%P2 e) ii) 
% Desplazamientos laterales (en cm)
% {dxe} = {K}^-1*{Fk}
dxe1 = K^-1*Fk1;  %cm

% P2 e) iii)
% Razón de derivadas de piso (Expresadas como porcentaje)

drift1 = dxe1 - [dxe1(2:cant_pisos); 0];   %cm
razon_drift1 = drift1/h;  %cm/cm

% comentarios y figuras
figure
plot(corte_pisos1,Pisos,'-o')
xlabel('V_{Pisos,1} [tonf]')
ylabel('Pisos')
grid on
disp('Corte de cada piso (modo 1)')
disp(corte_pisos1)

figure
plot(dxe1,Pisos,'-o')
xlabel('Desplazamaientos laterales de cada piso [cm] (modo 1)')
ylabel('Pisos')
grid on
disp('Desplazamientos laterales de cada piso [cm] (modo 1)')
disp(dxe1)

figure
plot(drift1,Pisos,'-o')
xlabel('Derivas de piso (drift) [cm] (modo1)')
ylabel('Pisos')
grid on
disp('Derivas de piso [cm] (modo 1)')
disp(drift1)

figure
plot(razon_drift1,Pisos,'-o')
xlabel('Razón de derivas de piso (modo 1)')
ylabel('Pisos')
grid on
disp('Razón de derivas de piso [%] (modo 1)')
disp(razon_drift1*100)
%% P3
% Análisis modal espectral, registro SantiagoCentro_Maule2010
% Utilizar 4 modos
n_modos = 4;                                                                      % Cantidad de modos a utilizar
fprintf('Cantidad de modos a considerar  n = %0.f \n',n_modos)

figure
hold on
for i = 1:n_modos
    [Sd,Sv,Sa,PSv,PSa] = Newmark_Lineal(beta,xi(i),dt,ui,udi,Registro.data,T_xlabel,Mn(i)); % cm y segundos
    modo.PSa(:,i) = PSa;    % cm/s2
    modo.PSv(:,i) = PSv;    % cm/s
    modo.Sa(:,i) = Sa;      % cm/s2
    modo.Sv(:,i) = Sv;      % cm/s
    modo.Sd(:,i) = Sd;      % cm
    plot(T_xlabel,PSa/g)
end
hold off
xlabel('Tn')
ylabel('PSa_n')
legend(strcat('Modo ',string(1:n_modos)))
grid on

%% P3 a)
% Corte Basal Máximo 
% Lo calculamos  en P3 c)
%% P3 b)
% Fuerzas laterales de piso máximas

% Primero determinamos las fuerzas laterales máximas de cada piso
% La fuerza lateral máxima para cada modo para cada piso ocurre cuando
% Dn(t) = Sd(Tn)
% fs_max(j,n) = Gamma_n(n)*phi(j,n)*mj*PSa(T(n))
fs_max = zeros(cant_pisos,n_modos);
figure
hold on
for n = 1:n_modos
    for j = 1:cant_pisos
        fs_max(j,n) = Gamma(n)*Phi(j,n)*M(j,j)*PSa(round(Tn(1)/Tn_step));
    end
    plot(fs_max(:,n),Pisos,'-o')
end
hold off
xlabel('fs_max')
ylabel('Pisos')
legend(strcat('Modo ',string(1:n_modos)))
grid on

% Para encontrar la fuerza máxima de cada piso (combinando todos los modos)
% Combinación modal -> SRSS, CQC o ABSSUM

% SRSS
fs_SRSSmax = zeros(cant_pisos,1);
for j = 1:cant_pisos
    fs_SRSSmax(j,1) = sqrt(sumsqr(fs_max(j,:)));
end
figure
plot(fs_SRSSmax,Pisos,'-o')
xlabel('fs_{SRSS}^{max}')
ylabel('Pisos')
grid on

%% P3 c)
% Máximo esfuerzo de corte de piso
% Primero obtenemos el corte de piso para cada modo

V = zeros(cant_pisos,n_modos);
figure 
hold on
for n = 1:n_modos
    for j = 1:cant_pisos
        V(j,n) = sum(fs_max(1:j,n));    %tonf
    end
    plot(V(:,n),Pisos,'-o')
end
hold off
xlabel('V_{jn}')
ylabel('Pisos')
legend(strcat('Modo ',string(1:n_modos)))
grid on

% Luego, SSRS

%% P3 d) 
% Máximo desplazamiento lateral de cada piso
% SRSS
%% P3 e)
% Máxima razón de desplazamiento de piso