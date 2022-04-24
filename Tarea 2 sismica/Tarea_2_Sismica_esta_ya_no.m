%% Inicializar
close all
clear variables
clc

%% P1

syms r_var C1_rock C2_rock C3_rock C4_rock C5_rock C1_soil C2_soil C3_soil C4_soil C5_soil H_var Zt_var Mw_var
% GMPE Para sitios en Roca
PGA_rock = exp(0.2418 + 1.414*Mw_var + C1_rock + C2_rock*(10-Mw_var)^3 + C3_rock*log(r_var + 1.7818*exp(0.544*Mw_var)) + 0.00607*H_var + 0.3846*Zt_var);
sigma_ln_PGA_rock = C4_rock + C5_rock*Mw_var;

% GMPE para sitios en suelo
PGA_soil = exp(-0.6687 + 1.438*Mw_var + C1_soil + C2_soil*(10-Mw_var)^3 + C3_soil*log(r_var + 1.097*exp(0.617*Mw_var)) + 0.00648*H_var + 0.3643*Zt_var);
sigma_ln_PGA_soil = C4_soil + C5_soil*Mw_var;

% Zt: 
% Mw: Moment Magnitude
% y: PGA: Peak ground acceleration
% r: Closest distance to rupture zone in km
% H : Focal depth
% Zt: source-type indicator
%       Zt = 0, for Interface earthquake
%       Zt = 1 for intraslab earthquake
%* Standard deviation

% Para calcular el PGA se utiliza el periodo T = 0, ya que es simular una
% vara rígida

r = (10:1:500).';
Mw_vals = [7.0; 8.0; 8.8];
H = 37; %km
Zt = 0; % Interface earthquake

% Utilizando las constantes de la Tabla 7.4 Villaverde(2009)
% Para roca (rock)
C1_rock_val = 0.00;
C2_rock_val = 0.00;
C3_rock_val = -2.552;
C4_rock_val = 1.45; % * Standard deviation for magnitud greater than 8 is equal to the value for magnitud equal to 8
C5_rock_val = -0.1; % *

% Para suelo (soil)
C1_soil_val = 0.00;
C2_soil_val = 0.00;
C3_soil_val = -2.329;
C4_soil_val = 1.45;
C5_soil_val = -0.1;

% Subs(GMPE)
PGA_rock = subs(PGA_rock, [C1_rock C2_rock C3_rock C4_rock C5_rock H_var Zt_var], [C1_rock_val C2_rock_val C3_rock_val C4_rock_val C5_rock_val H Zt]);
PGA_soil = subs(PGA_soil, [C1_soil C2_soil C3_soil C4_soil C5_soil H_var Zt_var], [C1_soil_val C2_soil_val C3_soil_val C4_soil_val C5_soil_val H Zt]);

for i = 1:3     
    PGA_rock = subs(PGA_rock,Mw_var,Mw_vals(i,1));
    PGA_soil = subs(PGA_soil,Mw_var,Mw_vals(i,1));
    for j = 1:length(r)
        puntos_rock(j,i) = subs(PGA_rock,r_var,r(j,1));
        puntos_soil(j,i) = subs(PGA_soil,r_var,r(j,1));
    end

    fprintf('Mw %f Roca \n',Mw_vals(i,1))
    figure
    loglog(r,puntos_rock(:,i))
    xlabel('r [km]')
    ylabel('PGA [unidad]')

    fprintf('Mw %f Suelo \n',Mw_vals(i,1))
    figure
    loglog(r,puntos_soil(:,i))
    xlabel('r [km]')
    ylabel('PGA [unidad]')
end

%% P2
