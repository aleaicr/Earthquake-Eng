clear all
close all
clc


%% P3
%% P3.a
% Definir fRM(r)
r_min = 7.2; % km
L1 = 34.4; % km
L2 = 62.3; % km
L = 96.7; % km

% Definr fRM(r)
% La vamos a definir en función de 'l' para que quede con los subcasos
syms r l
% fRM(r,l) = piecewise( ...
%     and(l<L1,L1<L2),piecewise(r == r_min,l/(L-l),and(r<r1,r>r_min),2*(L1-l)/(L-l),and(r<r2,r>r1),(L-2*L1)/(L-l)), ...
%     and(L1<l,l<L2),piecewise(r==r_min,L1/(L-l),and(r2>r,r>r_min),(L2-l)/(L-l)), ...
%     and(l>L2,L2>L1),1);

% % l<L1,L1<L2
% r11 = (r_min^2+(L1-l)^2)^0.5;
% r12 = (r_min^2+(L2-l)^2)^0.5;
% fRM1 = piecewise(r == r_min,l/(L-l),and(r<=r11,r>r_min),r/((L-l)*sqrt(r^2-r_min^2)),and(r<=r12,r>r11),r/((L-l)*sqrt(r^2-r_min^2)));
% % L1<l<L2
% r22 = (r_min^2+(L2-l)^2)^0.5;
% fRM2= piecewise(r==r_min,L1/(L-l),and(r22>=r,r>=r_min),r/((L-l)*sqrt(r^2-r_min^2)));
% % L1<L2<l
% fRM3 = piecewise(r==r_min,1,~(r==r_min),0);

%% P3.b
% En Excel
Mmax = 7.38;                                    % Magnitud de momento máxima (con L = 96.7) que la falla es capaz de producir

%% P3.c
% En Excel
l_vals = [43.72; 12.43; 6.62; 1.882]; %km
% magn = [0.95*Mmax; 0.85*Mmax; 0.8*Mmax;0.7*Mmax];

%% P3.d

for i = 1:length(l_vals)
    l = l_vals(i,1);
    l_val = l;
    if and(l_val<=L1,l_val<=L2)
        r1(i,1) = (r_min^2+(L1-l)^2)^0.5;
        r11 = r1(i,1);
        r2(i,1) = (r_min^2+(L2-l)^2)^0.5;
        r12 = r2(i,1);
        fRr(i,1) = piecewise(r == r_min,l/(L-l),and(r<=r11,r>r_min),r/((L-l)*sqrt(r^2-r_min^2)),and(r<=r12,r>r11),r/((L-l)*sqrt(r^2-r_min^2)));
    elseif and(l_val>L1, l_val<L2)
        r2(i,1) = (r_min^2+(L2-l)^2)^0.5;
        r22 = r2(i,1);
        fRr(i,1)= piecewise(r==r_min,L1/(L-l),and(r22>=r,r>=r_min),r/((L-l)*sqrt(r^2-r_min^2))); 
    elseif l_val>=L2
        fRr(i,1) = piecewise(r==r_min,1,~(r==r_min),0);                                        % Coordenada 1 es fR(r), no tiene ni r1 ni r2
    end
end

clear l_val r11 r12 r22
close all

figure
hold on
grid on
for i = 1:length(l_vals)
    fplot(fRr(i,1),[r_min r2(i,1)])
end
xlabel('r [km]')
ylabel('fRM(r)')
title('Funciones de densidad de probabilidad para las cuatro longitudes de ruptura')
legend('l = 43.72km; Magn = 7.01','l = 12.42km; Magn = 6.27','l = 6.62km; Magn = 5.90','l = 1.88km; Magn = 5.16')

%% P3.e

