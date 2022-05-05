function [u,ud,udd,fs,fluencias,ductilidad] = Newmark_No_Lineal(beta,Tn,xi,dt,ui,udi,uddg,Fy,R,t)
% beta -> Factor de newmark-beta
% Tn -> Periodo 
% xi -> Razón de amortiguamiento
% dt -> Paso temporal del registro
% ui -> Condición Inicial de desplazamiento
% udi -> Condición inicial de velocidad
% uddg -> Registro de aceleraciones del suelo
% Fy -> Resistencia lateral
% R ->
% t ->
gamma = 0.5;
u = zeros(size(t,1),1);
ud = zeros(size(t,1),1);
udd = zeros(size(t,1),1);
fluencias = zeros(size(t,1),1);      % Fluencia
wn = 2*pi/Tn;
k = wn^2;
uy = Fy/k;
u(1,1) = ui;
ud(1,1) = udi;
if abs(ui) <= uy
    udd(1,1) = uddg(1,1) - 2*xi*wn*udi - wn^2*ui;
else
    udd(1,1)=uddg(1,1) - 2*xi*wn*udi - Fy*sign(ui);
end
fs(1,1) = k*u(1,1);
a1 = 1/(beta*dt^2) + 2*xi*wn*gamma/(beta*dt);
a2 = 1/(beta*dt)+2*xi*wn*(gamma/beta-1);
a3 = (1/(2*beta)-1)+2*xi*wn*dt*(gamma/(2*beta)-1);
kt = k;
for i = 1:size(t,1)-1
    u(i+1,1) = u(i,1);                                                       % Valor semilla de iteracion j.
    fs(i+1,1) = fs(i,1);                                                     % Estado de fs.
    p_ton = -uddg(i+1,1) + a1*u(i,1) + a2*ud(i,1) + a3*udd(i,1);
    R_ton = p_ton - fs(i+1,1) - a1*u(i+1,1);
    while abs(R_ton) > R                                                     % Revision de convergencia N-R.
        fluencias(i+1,1) = 0;
        kt_ton = kt + a1;
        Delta_u = R_ton/kt_ton;
        u(i+1,1) = u(i+1,1) + Delta_u;
        fs(i+1,1) = fs(i+1,1) + k*Delta_u;
        if abs(fs(i+1,1)) > Fy
            fs(i+1,1) = Fy*sign(fs(i+1,1));
            fluencias(i+1,1) = sign(fs(i+1,1));
            kt = k;
        else
            kt = k;
        end
        R_ton = p_ton - fs(i+1,1) - a1*u(i+1,1);
    end
        
    ud(i+1,1) = gamma/(beta*dt)*(u(i+1,1)-u(i,1)) + (1-gamma/beta)*ud(i,1) + dt*(1-gamma/(2*beta))*udd(i,1);
    udd(i+1,1) = (u(i+1,1)-u(i,1))/(beta*dt^2) - ud(i,1)/(beta*dt) - (1/(2*beta)-1)*udd(i,1);
end
ductilidad = max(abs(u))/uy;
end