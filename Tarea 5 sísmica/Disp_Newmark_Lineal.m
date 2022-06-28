function [u] = Disp_Newmark_Lineal(beta,xi,dt,ui,udi,uddg,Tn)
% Beta  -> Factor beta del método de Newmark (1/6 o 1/4)
% xi    -> Razón de amortiguamiento del modo
% dt    -> Rango de sampling del registro
% ui    -> Desplazamiento inicial
% udi   -> Velocidad Inicial
% uddg  -> Registro de aceleraciones del suelo
% Tn    -> Periodo a analizar (solo 1)

gamma = 0.5;
uddg_length = length(uddg);

u = zeros(uddg_length,1);
ud = zeros(uddg_length,1);
udd = zeros(uddg_length,1);

% Condiciones iniciales
u(1,1) = ui;
ud(1,1) = udi;

% 
wn = 2*pi/Tn;                                                               
udd(1,1) = uddg(1,1) - 2*xi*wn*udi - wn^2*ui;
a1 = 1/(beta*dt^2) + 2*xi*wn*gamma/(beta*dt);
a2 = 1/(beta*dt) + 2*xi*wn*(gamma/beta-1);
a3 = (1/(2*beta)-1) + 2*xi*wn*dt*(gamma/(2*beta)-1);
k_ton = a1 + wn^2;

for i = 2:uddg_length
    p_ton = -uddg(i,1) + a1*u(i-1,1) + a2*ud(i-1,1) + a3*udd(i-1,1);
    u(i,1) = p_ton/k_ton;
    ud(i,1) = gamma/(beta*dt)*(u(i,1)-u(i-1,1)) + (1-gamma/beta)*ud(i-1,1) + dt*(1-gamma/(2*beta))*udd(i-1,1);
    udd(i,1) = (u(i,1)-u(i-1,1))/(beta*dt^2) - ud(i-1,1)/(beta*dt) - (1/(2*beta)-1)*udd(i-1,1);
end
end