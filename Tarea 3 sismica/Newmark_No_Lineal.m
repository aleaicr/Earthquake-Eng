function [u,ud,udd,fs,yield,ductility] = Newmark_No_Lineal(beta,Tn,xi,dt,ui,udi,uddg,Fy,R,t_vect)
% beta  -> Factor de newmark-beta
% Tn    -> Periodo 
% xi    -> Razón de amortiguamiento
% dt    -> Paso temporal del registro
% ui    -> Condición Inicial de desplazamiento
% udi   -> Condición inicial de velocidad
% uddg  -> Registro de aceleraciones del suelo
% Fy    -> Resistencia lateral
% R     -> Para criterio de convergencia
% t     -> Vector de tiempos

gamma = 0.5;
% Inicialización de desplazamientos, velocidades, aceleraciones y fluencias
l = length(t_vect);
u = zeros(l,1);
ud = zeros(l,1);
udd = zeros(l,1);
yield = zeros(l,1);                                                 % Fluencia

wn = 2*pi/Tn;                                                               % Frecuencia natural del oscilador 
k = wn^2;                                                                   % wn2
uy = Fy/k;

% Condiciones uniciales en el primer paso temporal
u(1,1) = ui;
ud(1,1) = udi;

if abs(ui) <= uy                                                            % Si el desplazamiento inicial es menor o igual al desplazamiento de fluencia
    udd(1,1) = uddg(1,1) - 2*xi*wn*udi - wn^2*ui;                           % La aceleración incial es udd = uddg - 2*xi*wn*udi - wn^2*ui
else
    udd(1,1) = uddg(1,1) - 2*xi*wn*udi - Fy*sign(ui);                       % Si es mayor, entonces es udd = uddg - 2*xi*wn*udi - Fy*signo(ui)
end

fs(1,1) = k*u(1,1);                                                         % Las fuerzas resistentes serán la rigidez por el desplazamiento

% Coeficientes a1, a2, a3 y kt
a1 = 1/(beta*dt^2) + 2*xi*wn*gamma/(beta*dt);                               
a2 = 1/(beta*dt) + 2*xi*wn*(gamma/beta-1);
a3 = (1/(2*beta)-1) + 2*xi*wn*dt*(gamma/(2*beta)-1);
kt = k;

for i = 1:l-1
    u(i+1,1) = u(i,1);                                                      % Valor semilla de iteracion j.
    fs(i+1,1) = fs(i,1);                                                    % Estado de fs.
    p_ton = -uddg(i+1,1) + a1*u(i,1) + a2*ud(i,1) + a3*udd(i,1);            % P tongo
    R_ton = p_ton - fs(i+1,1) - a1*u(i+1,1);                                % R tongo
    while abs(R_ton) > R                                                    % Revision de convergencia N-R.
        yield(i+1,1) = 0;                                               
        kt_ton = kt + a1;
        Du = R_ton/kt_ton;
        u(i+1,1) = u(i+1,1) + Du;
        fs(i+1,1) = fs(i+1,1) + k*Du;
        if abs(fs(i+1,1)) > Fy
            fs(i+1,1) = Fy*sign(fs(i+1,1));
            yield(i+1,1) = sign(fs(i+1,1));
            kt = k;
        else
            kt = k;
        end
        R_ton = p_ton - fs(i+1,1) - a1*u(i+1,1);
    end
        
    ud(i+1,1) = gamma/(beta*dt)*(u(i+1,1)-u(i,1)) + (1-gamma/beta)*ud(i,1) + dt*(1-gamma/(2*beta))*udd(i,1);
    udd(i+1,1) = (u(i+1,1)-u(i,1))/(beta*dt^2) - ud(i,1)/(beta*dt) - (1/(2*beta)-1)*udd(i,1);
end
ductility = max(abs(u))/uy;
end