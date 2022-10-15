%% Examen Final Ingeniería Sísmica
% Cristobal Adasme - Alexis Contreras

%% Inicializar
clear variables
close all
clc

%% Propiedades
T1 = 1.5; % sec                                                             % Periodo buscado
FE = 1;                                                                     % Factor de escala = 1;
cant_pisos = 10;
Pisos = (cant_pisos:-1:1)';
Modos = flip(Pisos);
cant_modos = length(Modos);
L10 = 4.60; % m
L = 3.65; % m
lfact = L10/L; % L10 = lfact*L - Para escalar el primer piso, si primer piso es igual a tados lfact debe ser 1 (cambiar)
W_techo = 400.35; % kN
W_piso = 480.4; % kN
Wtotal = W_techo + W_piso*(cant_pisos-1); %kN
% g = 980; % cm/s2
g = 9.8; % m/s2

% Alturas
h = [L*ones(cant_pisos-1,1); 4.60];                                         % Altura de piso
h_base = zeros(cant_pisos,1);                                               % Altura desde la base al piso 
for i = 1:cant_pisos
    h_base(i) = sum(h(cant_pisos:-1:i,1));
end
% Masa
m1 = W_techo; % kN
m2 = W_piso;
m3 = W_piso;
m4 = W_piso;
m5 = W_piso;
m6 = W_piso;
m7 = W_piso;
m8 = W_piso;
m9 = W_piso;
m10 = W_piso;

M_tilde = diag([m1 m2 m3 m4 m5 m6 m7 m8 m9 m10]);
% Rigidez

k1 = 0.96; %kN/m   (Se escala el primer piso que tiene altura diferente)
k2 = 0.96;
k3 = 1.32;
k4 = 1.32;
k5 = 1.82;
k6 = 1.82;
k7 = 2.22;
k8 = 2.22;
k9 = 3.3;
k10 = 3.3/lfact^3;

K_tilde = [k1 -k1 0 0 0 0 0 0 0 0;
    -k1 k2+k1 -k2 0 0 0 0 0 0 0;
    0 -k2 k3+k2 -k3 0 0 0 0 0 0;
    0 0 -k3 k4+k3 -k4 0 0 0 0 0;
    0 0 0 -k4 k5+k4 -k5 0 0 0 0;
    0 0 0 0 -k5 k6+k5 -k6 0 0 0;
    0 0 0 0 0 -k6 k7+k6 -k7 0 0;
    0 0 0 0 0 0 -k7 k8+k7 -k8 0;
    0 0 0 0 0 0 0 -k8 k9+k8 -k9;
    0 0 0 0 0 0 0 0 -k9 k10+k9];

% Vector de influencia (?)
G = [-m1; -m2; -m3; -m4; -m5; -m6; -m7; -m8; -m9; -m10];
r = diag(ones(cant_pisos));

%% PREGUNTA 1 Obtención de EI

% Valores propios
[Phi_tilde,gamma_tilde] = eig(K_tilde,M_tilde);
gamma_tilde = diag(gamma_tilde);

% Evaluamos primer modo para encontrar rigidez que provoca un periodo del
% primer modo igual a 1.5
EI = (2*pi/T1)^2*(L^3)/(12*gamma_tilde(1)*g); %kN/m2

% Obteniendo los reales
K = 12*EI/L^3*K_tilde;                                                      % Rigidez
M = 1/g*M_tilde;                                                            % Masa

[Phi,lambda] = eig(K,M);
lambda = diag(lambda);                                                      % Valores propios

wn = lambda.^0.5;                                                           % Frecuencia angular de cada modo
wn2 = wn.^2;
Tn = 2*pi./wn;                                                              % Periodo de cada modo

% Matrices Modales
Mn = diag(Phi.'*M*Phi);                                                     % Masa modal
Kn = diag(Phi.'*K*Phi);                                                     % Rigidez modal

% Gamma_n
Gamma_n = Phi.'*M*r./Mn;

% Graficamos Gamma_n*phi_n para 5 modos Todos separados
cant_modos_dibujar = 5;                                                     % Cantidad de modos que queremos dibujar (cambiar si se quieren dibujar más)
figure
for i = 1:cant_modos_dibujar
    subplot(1,cant_modos_dibujar,i)
    plot([Gamma_n(i)*Phi(:,i); 0],[Pisos; 0], '-o')
    grid on
    xlabel(['\Gamma_', num2str(Modos(i)), '\{\phi_{', num2str(Modos(i)), '}\}'])
    xlim([-1.5 1.5])
    if i == 1
        ylabel('Piso')
        hold on
    end
end
hold off
sgtitle('\Gamma_n \{\phi_n\}')

% Graficamos Gamma_n*phi_n para 5 modos Todos juntos
figure
hold on
for i = 1:cant_modos_dibujar
    plot([Gamma_n(i)*Phi(:,i); 0],[Pisos; 0], '-o')
end
hold off
xlabel('\Gamma_n \{\Phi_n\}')
ylabel('Piso')
legend(strcat('Modo ',string(1:cant_modos_dibujar)))
grid on

% Comentarios Pregunta 1
disp('Pregunta 1')
disp('EI [kN/cm2]'), disp(EI*100^2)
disp('K'), disp(K)
disp('M'), disp(M)
tabla = table();
tabla.Modos = Modos;
tabla.ValPropios = lambda;
tabla.Tn = Tn;
tabla.Gamma_n = Gamma_n;
disp(tabla)
clear tabla

%% Pregunta 2
% PSa NCh433 of1996 Mod2009
T_init = 0;
T_step = 0.01;
T_final = 3;
T_rango = T_init:T_step:T_final;
Ro = 11;                                                                    % Marco de acer
I = 1;                                                                      % Categoría 2 (Edificio Habitacional)
Ao = 0.3;                                                                   % Zona sísmica 2 (Providencia)                                                               % Marco especial de hacero
To = 0.30;                                                                    % Tipo de suelo 2
p = 1.5;                                                                    % Tipo de suelo 2

T_ast = Tn(1);                                                              % Periodo con modo con mayor masa traslacional
R_ast = 1 + Tn(1)/(0.10*To + Tn(1)/Ro);                                     % R_ast
Sa = zeros(length(T_rango),1);                                              % Sa = I*A0*alfa/R_ast
Sa_red = zeros(length(T_rango),1);                                          % Sa/R_ast

for i = 1:length(T_rango)
    alfa = (1 + 4.5*(T_rango(i)/To)^p)/(1 + (T_rango(i)/To)^3);
    Sa(i) = I*Ao*alfa;                                                      % PSa sin reducir NCh433 of1996 Mod2009, NO TIENE factor 'S'                     

    Sa_red(i) = Sa(i)/R_ast;                                                % PSa reducido (R_ast con T* = T1 (del primer modo))
end

figure
hold on
plot(T_rango,Sa)
plot(T_rango,Sa_red)
hold off
xlabel('Periodo [sec]')
ylabel('PSa = Sa NCh433 [g]')
legend('PSa','PSa/R*')
title('Espectro de Diseño NC433 of1996 Mod.2012')
grid on

%% Pregunta 3
% Edificio del prediseño, analisis modal espectral, SRSS, xi = 5% para
% todos los modos, considerar solo 5 modos

% Mn_ast = zeros(cant_pisos,1);
% for n = 1:cant_pisos
%     for j = 1:cant_pisos
%         Mn_ast(n,1) = Mn_ast(n,1) + M(j,j)*Gamma_n(n)*Phi(j,n);
%     end
% end
% 
% PMass = Mn_ast/(Wtotal/g);                                                  % Porcentaje de Masa modal Equivalente de cada modo
% PMass_acum = zeros(cant_pisos,1);                                           % Porcentaje de Masa modal equivalente acumulado
% PMass_acum(1,1) = PMass(1,1);                                               
% counter = 0;                                                                % Contador para ver si ya encontró la solución
% for i = 2:cant_pisos
%     PMass_acum(i,1) = PMass_acum(i-1,1) + PMass(i,1);
%     if PMass_acum(i,1) > 0.9 && counter == 0
%         fprintf('Solo se requieren los primeros %i modos \n \n', Modos(i,1))
%         counter = 1;
%     end
% end

% Notar que Utilizando NCh433 (Sum(Mn_ast_i) > 90%W_total) se tiene que se 
% requieren solo 2 modos.

n_modos = 5;
xi_all = 0.05;
xi = xi_all*ones(cant_modos);

PSa = zeros(n_modos,1);

for i = 1:n_modos
    n = round(Tn(i)/T_step) + 1;
    PSa(i) = Sa_red(n)*g;
end
disp('PSa [m/s2] (Sa/R*)')
disp(PSa)

run Analisis_modal_espectral.m

% %% Análisis Modal Espectral
% %% Fuerzas Máximas de cada modo 
% fs_max = zeros(cant_pisos,n_modos);
% figure
% hold on
% for n = 1:n_modos
%     for j = 1:cant_pisos
%         fs_max(j,n) = Gamma_n(n)*Phi(j,n)*M(j,j)*PSa(n);
%     end
%     plot(fs_max(:,n),Pisos,'-o')
% end
% hold off
% xlabel('fs_{jn}^{max}')
% ylabel('Pisos')
% legend(strcat('Modo ',string(1:n_modos)))
% grid on
% 
% % SRSS
% fs_SRSSmax = zeros(cant_pisos,1);
% for j = 1:cant_pisos
%     fs_SRSSmax(j,1) = sqrt(sumsqr(fs_max(j,:)));
% end
% figure
% plot(fs_SRSSmax,Pisos,'-o')
% xlabel('fs_{SRSS}^{max} [kN]')
% ylabel('Pisos')
% grid on
% disp('fs_SRSSmax [kN]')
% disp(fs_SRSSmax)
% 
% %% Máximo desplazamiento lateral de cada piso SRSS
% 
% u_jnmax = zeros(cant_pisos,n_modos);
% for n = 1:n_modos
%     for j = 1:cant_pisos
%         u_jnmax(j,n) = Gamma_n(n)*Phi(j,n)*PSa(n)*(1/wn2(n));
%     end
% end
% 
% % Luego, SSRS
% u_SRSSmax = zeros(cant_pisos,1);
% for j = 1:cant_pisos
%     u_SRSSmax(j,1) = sqrt(sumsqr(u_jnmax(j,:)));
% end
% figure
% plot(u_SRSSmax*100,Pisos,'-o')
% xlabel('u_{SRSS}^{max} [cm]')
% ylabel('Pisos')
% grid on
% disp('u_SRSSmax [cm]')
% disp(u_SRSSmax*100)
% 
% %% Máxima razón de desplazamiento de piso
% 
% rdp = zeros(cant_pisos,n_modos);
% for j = 1:cant_pisos
%     for n = 1:n_modos
%         if j < cant_pisos
%             rdp(j,n) = (1/h(j))*Gamma_n(n)*(Phi(j,n)-Phi(j+1,n))*PSa(n)*(1/wn2(n));
%         elseif j == cant_pisos
%             rdp(j,n) = (1/h(j))*Gamma_n(n)*(Phi(j,n))*PSa(n)*(1/wn2(n));
%         end
%     end
% end
% 
% % Luego SRSS
% 
% rdp_SRSSmax = zeros(cant_pisos,1);
% for j = 1:cant_pisos
%     rdp_SRSSmax(j,1) = sqrt(sumsqr(rdp(j,:)));
% end
% figure
% plot(rdp_SRSSmax*100,Pisos,'-o')
% xlabel('rdp_{SRSS}^{max} [%]')
% ylabel('Pisos')
% grid on
% disp('rdp_SRSSmax [%]')
% disp(rdp_SRSSmax*100)
% 
% %% Máximo esfuerzo de corte de piso
% % Primero obtenemos el corte de piso para cada modo
% 
% V = zeros(cant_pisos,n_modos);
% figure 
% hold on
% for n = 1:n_modos
%     for j = 1:cant_pisos
%         V(j,n) = sum(fs_max(1:j,n));    %tonf
%     end
%     plot(V(:,n),Pisos,'-o')
% end
% hold off
% xlabel('V_{jn} [kN]')
% ylabel('Pisos')
% legend(strcat('Modo ',string(1:n_modos)))
% grid on
% 
% % Luego, SSRS
% V_SRSSmax = zeros(cant_pisos,1);
% for j = 1:cant_pisos
%     V_SRSSmax(j,1) = sqrt(sumsqr(V(j,:)));
% end
% figure
% plot(V_SRSSmax,Pisos,'-o')
% xlabel('V_{SRSS}^{max} [kN]')
% ylabel('Pisos')
% grid on
% disp('V_SRSSmax [kN]')
% disp(V_SRSSmax)
% 
% %% Corte basal máximo
% Vb_max = zeros(1,n_modos);
% for n = 1:n_modos
%     Vb_max(1,n) = sum(fs_max(:,n));
% end
% 
% % SRSS
% Vb_SRSSmax = sqrt(sumsqr(Vb_max(1,:)));
% disp('Vb_SRSSmax [kN]')
% disp(Vb_SRSSmax)
% 
% %% Momento basal máximo
% Mn_max = zeros(j,n_modos);
% for n = 1:n_modos
%     for j = 1:cant_pisos
%         Mn_max(j,n) = h(j)*fs_max(j,n);
%     end
% end
% 
% Mb_max = zeros(1,n_modos);
% for n = 1:n_modos
%     Mb_max(1,n) = sum(Mn_max(:,n));
% end
% 
% % SRSS
% Mb_SRSSmax = sqrt(sumsqr(Mb_max(1,:)));
% disp('Mb_SRSSmax [kN-m]')
% disp(Mb_SRSSmax)

%% PREGUNTA 4
% Comparar resultados de corte basal anterior con los requisitos de corte
% basal máximo y mínimo de la norma. Determinar si es necesario escalar o
% no la demanda sísmica e indique el factor de escala a utilizar (de ser
% necesario) si es necesario escalar recalcule los resultados de la parte
% anterior -> FUCK xd

% Corte basal norma
P = Wtotal;
R = 7;
T_prima = 0.35;
n = 1.33;
S = 1;

C = 2.75*Ao/R*(T_prima/T_ast)^n;
Qo = C*I*P;

% Máximo
Cmax = 0.35*S*Ao;
Qmax = Cmax*I*P;

% Mínimo
Cmin = Ao/6;
Qmin = Cmin*I*P;

disp('Qmin [kN]'), disp(Qmin)
disp('Qo [kN]'), disp(Qo)
disp('Qmax [kN]'), disp(Qmax)

if Vb_SRSSmax < Qmin
    fprintf('Se aprecia que Vb (A.Modal) con Qo (A.Estático)')
    fprintf('Se diferencian en %f [kN] \n',Qo-Vb_SRSSmax)
    fprintf('El corte basal es menor al mínimo requerido,\n por lo que se debe modificar \n')
    fprintf('El espectro de respuesta por Qmin/Qo \n \n')
    FE = Qmin/Vb_SRSSmax;                                                   % Factor de escala
    PSa = PSa*FE;                                                           % Aplicamos factor de escala
    run Analisis_modal_espectral.m                                          % Retorna gráficos y resultados en la consola
end

%% PREGUNTA 5
% Hacer figura graficar razon derivas piso en altura
%

delta = 0.002*ones(cant_pisos,1);
figure
plot(rdp_SRSSmax,Pisos,'-o')
hold on
plot(delta,Pisos,'--')
hold off
grid on
xlim([0 0.0025])
legend('rdp (\delta)','rdp_{max} (\delta_{max})')
xlabel('Razón de deriva de pisos')
ylabel('Piso')
title('Razón de deriva de pisos')

counter = 0;
for i = 1:cant_pisos
    if rdp_SRSSmax(i) > delta
        fprintf('Razón de deriva de pisos %.f > delta \n',rdp_SRSSmax(i))
        counter = counter + 1;
    end
end

if counter == 0
    fprintf('Como se observa en la figura, ningúna razón de deriva de pisos sobrepasa el valor límite de 0.002 \n')
end

