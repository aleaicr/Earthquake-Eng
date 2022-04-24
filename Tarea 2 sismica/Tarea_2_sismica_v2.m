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
rmin = 7.2; % km
L1 = 34.4; % km
L2 = 62.3; % km
L = 96.7; % km

% Definr fRM(r)
% La vamos a definir en función de 'l' para que quede con los subcasos
syms r l
fRM(r,l) = piecewise(and(l<L1,L1<L2),piecewise(r == rmin,l/(L-l),and(r<r1,r>rmin),2*(L1-l)/(L-l),and(r<r2,r>r1),(L-2*L1)/(L-l)),and(L1<l,l<L2),piecewise(r==rmin,L1/(L-l),and(r2>r,r>rmin),(L2-l)/(L-l)),and(l>L2,L2>L1),1);

% l<L1,L1<L2
fRM1(r) = piecewise(r == rmin,l/(L-l),and(r<r1,r>rmin),2*(L1-l)/(L-l),and(r<r2,r>r1),(L-2*L1)/(L-l));
% L1<l<L2
fRM2(r)= piecewise(r==rmin,L1/(L-l),and(r2>r,r>rmin),(L2-l)/(L-l));
% L1<L2<l
fRM3(r) = piecewise(r==rmin,1,~(r==rmin),0);

Mmax = 7.38;
l_vals = [43.72; 12.43; 6.62; 1.882]; %km
magn = [0.95*Mmax; 0.85*Mmax; 0.8*Mmax;0.7*Mmax];

for i = 1:4
    l_val = l_vals(i,1);
    if and(l_val<L1,l<L2)
        funct1 = fRM1
    elseif and(l_val>L1, l<L2)
        
    elseif l_val>L2
    end


end



%% P3.b
% En Excel

%% P3.c
% En Excel

%% P3.d


%% P3.e

