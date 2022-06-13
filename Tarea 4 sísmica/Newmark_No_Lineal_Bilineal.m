function [u, up, upp, fs, yield, mu] = Newmark_No_Lineal_Bilineal(alpha, beta, Tn, xi, dt, u0, u0p, uppg, Fy, R, t)
% Método de Newmark para Modelo Bilineal (No Lineal)
% alfa: Pendiente de degradación o endurecimiento en modelo bilineal
% beta: Parámetro de Método de Newmark
% Tn, xi: Parámetros de dinámica de estructura
% dt: Paso temporal
% u0, u0p: Condiciones iniciales
% uppg: Registro de movimiento sísmico
% R: Criterio de Convergencia
% t: Vector de tiempo para el registro

% Entrega
% u: Historia de desplazamientos relativos
% up: Historia de velocidades relativas
% upp: Historia de aceleraciones absolutas
% fs: Historia de valores de fuerza resistente
% mu: Demanda de ductilidad

% Parámetros Método de Newmark
gamma = 0.5;                                                                % Gamma de Métodos de Newmark
t_length = length(t);                                                       % Largo del vector de tiempo

% Inicialización de vectores
u = zeros(t_length,1);
up = zeros(t_length,1);
upp = zeros(t_length,1);
yield = zeros(t_length,1);

%Propiedades
omega_n = 2*pi/Tn;                                                          % Frecuencia natural estructura
k = omega_n^2;                                                              % Rigidez Estructura (m = 1)
uy = Fy/k;                                                                  % Desplazamiento de fleuncia estructura

% Condiciones iniciales
u(1,1) = u0;
up(1,1) = u0p;

if abs(u0) <= uy
    upp(1,1) = uppg(1,1) - 2*xi*omega_n*u0p - omega_n^2*u0;
else
    upp(1,1) = uppg(1,1) - 2*xi*omega_n*u0p - Fy*sign(u0);
end

% Método de Newmark
fs(1,1) = k*u(1,1);                                                         % Resistencia Condición inicial
a1 = 1/(beta*dt^2) + 2*xi*omega_n*gamma/(beta*dt);                          % Factores a1, a2, a3 del método de Newmark
a2 = 1/(beta*dt) + 2*xi*omega_n*(gamma/beta-1);
a3 = (1/(2*beta)-1) + 2*xi*omega_n*dt*(gamma/(2*beta)-1);
kt = k;                                                                     % K_tongo

for i = 1:length(t) - 1
    j_counter = 0;                                                          % Contador para limitar cantidad de iteraciones para cada i (paso temporal), por si no converge
    u(i+1,1) = u(i,1);                                                      % Valor semilla de iteracion j.
    fs(i+1,1) = fs(i,1);                                                    % Estado de fs.
    p_tongo = -uppg(i+1,1) + a1*u(i,1) + a2*up(i,1) + a3*upp(i,1);          % p_tongo = -uppg + a1*u_i + u2*up_i
    R_tongo = p_tongo - fs(i+1,1) - a1*u(i+1,1);                            % R_tongo = p_tongo - fs_i+1 - a1*u_i+1
    while abs(R_tongo) > R                                                  % Revision de convergencia N-R.
        yield(i+1,1) = 0;
        kt_tongo = kt + a1;
        Delta_u = R_tongo/kt_tongo;
        u(i+1,1) = u(i+1,1) + Delta_u;
        fs(i+1,1) = fs(i+1,1) + k*Delta_u;
        if abs(fs(i+1,1) - u(i+1,1)*alpha*k) > Fy*(1-alpha)
            fs(i+1,1) = Fy*sign(fs(i+1,1))*(1-alpha) + u(i+1,1)*alpha*k;    % Parte EPP + Parte Elástica(Lineal)
            yield(i+1,1) = sign(fs(i+1,1));
            kt = alpha*k;
        else
            kt = k;
        end
        if j_counter >= 10000
            break                                                           % Break si no converge en el while
        end
        R_tongo = p_tongo - fs(i+1,1) - a1*u(i+1,1);                        % R_tongo = p_tongo - fs_i+1 - a1*u_i+1
        j_counter = j_counter + 1;
    end  
    up(i+1,1) = gamma/(beta*dt)*(u(i+1,1)-u(i,1)) + (1-gamma/beta)*up(i,1) + dt*(1-gamma/(2*beta))*upp(i,1);
    upp(i+1,1) = (u(i+1,1)-u(i,1))/(beta*dt^2) - up(i,1)/(beta*dt) - (1/(2*beta)-1)*upp(i,1);
end
mu = max(abs(u))/uy;
end