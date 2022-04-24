%% Inicializar
clear variables
close all
clc
%% P1
%% P1.1
% Graficar PGA promedio esperado para un terremoto tipo interplaca vs la
% distancia r, con las condiciones
% H = 37 km; r = 10:1:500; Mw = [7.0;8.0;8.8]
% Comparar resultados
C1_rock = 0.00;
C2_rock = 0.00;
C3_rock = -2.552;
C4_rock = 1.45; % * Standard deviation for magnitud greater than 8 is equal to the value for magnitud equal to 8
C5_rock = -0.1; % *

% Para suelo (soil)
C1_soil = 0.00;
C2_soil = 0.00;
C3_soil = -2.329;
C4_soil = 1.45;
C5_soil = -0.1;

Zt_var = 0;
H_var = 37;
Mw_vector = [7.0; 8.0; 8.8];
r = (10:1:500).';

for Mw = 1:3
    Mw_var = Mw_vector(Mw,1);
    for r_var = 1:length(r)
        puntos_rock(r_var,Mw) = exp(0.2418 + 1.414*Mw_var + C1_rock + C2_rock*(10-Mw_var)^3 + C3_rock*log(r(r_var,1) + 1.7818*exp(0.544*Mw_var)) + 0.00607*H_var + 0.3846*Zt_var);
        puntos_soil(r_var,Mw) = exp(-0.6687 + 1.438*Mw_var + C1_soil + C2_soil*(10-Mw_var)^3 + C3_soil*log(r(r_var,1) + 1.097*exp(0.617*Mw_var)) + 0.00648*H_var + 0.3643*Zt_var);
    end
    fprintf('Mw %f Roca \n',Mw_var)

    figure
    loglog(r,puntos_rock(:,Mw))
    xlabel('r [km]')
    ylabel('PGA roca[unidad]')
    title('PGA para roca, sismo magnitud', Mw_var)

    fprintf('Mw %f Suelo \n',Mw_var)
    figure
    loglog(r,puntos_soil(:,Mw))
    xlabel('r [km]')
    ylabel('PGA suelo[unidad]')
    title('PGA para suelo, sismo magnitud', Mw_var)
end


%% P1.2
% Para un periodo de $T_n = 2.0 [sec]$ y un amortiguamiento de $\xi = 0.05$

% Para roca
C1_rock = -3.328;
C2_rock = -0.0080;
C3_rock = -2.1070;
C4_rock = 1.55;
C5_rock = -0.1;

% Para suelo
C1_soil = -6.433;
C2_soil = -0.0164;
C3_soil = -1.290;
C4_soil = 1.55;
C5_soil = -0.1;

for Mw = 1:3
    Mw_var = Mw_vector(Mw,1);
    for r_var = 1:length(r)
        puntos_rock(r_var,Mw) = exp(0.2418 + 1.414*Mw_var + C1_rock + C2_rock*(10-Mw_var)^3 + C3_rock*log(r(r_var,1) + 1.7818*exp(0.544*Mw_var)) + 0.00607*H_var + 0.3846*Zt_var);
        puntos_soil(r_var,Mw) = exp(-0.6687 + 1.438*Mw_var + C1_soil + C2_soil*(10-Mw_var)^3 + C3_soil*log(r(r_var,1) + 1.097*exp(0.617*Mw_var)) + 0.00648*H_var + 0.3643*Zt_var);
    end
    fprintf('Mw %f Roca \n',Mw_var)

    figure
    loglog(r,puntos_rock(:,Mw))
    xlabel('r [km]')
    ylabel('PGA roca[unidad]')
    title('PGA para roca, sismo magnitud', Mw_var)

    fprintf('Mw %f Suelo \n',Mw_var)
    figure
    loglog(r,puntos_soil(:,Mw))
    xlabel('r [km]')
    ylabel('PGA suelo[unidad]')
    title('PGA para suelo, sismo magnitud', Mw_var)
end

%% P1.3
% No sé hacerla xd


%% P2
% Está en Excel completa


%% P3
%% P3.a
% Definir fRM(r)
r_min = 7.2; % km
L1 = 34.4; % km
L2 = 62.3; % km
L = 96.7; % km

% Definr fRM(r)
% La vamos a definir en función de 'l' para que quede con los subcasos
syms r l r1 r2
% fRM(r,l) = piecewise(and(l<L1,L1<L2),piecewise(r == r_min,l/(L-l),and(r<r1,r>r_min),2*(L1-l)/(L-l),and(r<r2,r>r1),(L-2*L1)/(L-l)),and(L1<l,l<L2),piecewise(r==r_min,L1/(L-l),and(r2>r,r>r_min),(L2-l)/(L-l)),and(l>L2,L2>L1),1);

% l<L1,L1<L2
fRM1 = piecewise(r == r_min,l/(L-l),and(r<=r1,r>r_min),r/((L-l)*sqrt(r^2-r_min^2)),and(r<=r2,r>r1),r/((L-l)*sqrt(r^2-r_min^2)));
r11 = (r_min^2+(L1-l)^2)^0.5;
r12 = (r_min^2+(L2-l)^2)^0.5;
% L1<l<L2
fRM2= piecewise(r==r_min,L1/(L-l),and(r2>=r,r>=r_min),r/((L-l)*sqrt(r^2-r_min^2)));
r22 = (r_min^2+(L2-l)^2)^0.5;
% L1<L2<l
fRM3 = piecewise(r==r_min,1,~(r==r_min),0);

%% P3.b
% En Excel
Mmax = 7.38;                                    % Magnitud de momento máxima (con L = 96.7) que la falla es capaz de producir

%% P3.c
% En Excel
l_vals = [43.72; 12.43; 6.62; 1.882]; %km
magn = [0.95*Mmax; 0.85*Mmax; 0.8*Mmax;0.7*Mmax];

%% P3.d

for i = 1:length(l_vals)
    l_val = l_vals(i,1);
    if and(l_val<=L1,l_val<L2)
        fRMj(i,2) = subs(r12,l,l_val);                                          % Coordenada 2 es el límite r2
        fRMj(i,3) = subs(r11,l,l_val);                                          % Coordenada 3 es r1
        fRMj(i,1) = subs(subs(subs(fRM1,l,l_val),r1,fRMj(i,3)),r2,fRMj(i,2));   % Coordenada 1 es fR(r), tiene r1 y r2
    elseif and(l_val>L1, l_val<=L2)
        fRMj(i,2) = subs(r22,l,l_val);                                          % Coordenada 2 es el límite r2  
        fRMj(i,1) = subs(subs(fRM2,l,l_val),r2,fRMj(i,2));                      % Coordenada 1 es fR(r), tiene r1, no tiene r2  
    elseif l_val>L2
        fRMj(i,1) = subs(fRM3,l,l_val);                                         % Coordenada 1 es fR(r), no tiene ni r1 ni r2
    end
end
clear l_val
close all
figure
hold on
grid on
fplot(fRMj(1,1),[r_min fRMj(1,2)])
fplot(fRMj(2,1),[r_min fRMj(2,2)])
fplot(fRMj(3,1),[r_min fRMj(3,2)])
fplot(fRMj(4,1),[r_min fRMj(4,2)])
xlabel('r [km]')
ylabel('fRM(r)')
% legend(['l = 43.72km; Magn = 7.01'; 'l = 12.42km; Magn = 6.27'; 'l = 6.62km; Magn = 5.90'; 'l = 1.88km; Magn = 5.16'])
hold off

%% P3.e

