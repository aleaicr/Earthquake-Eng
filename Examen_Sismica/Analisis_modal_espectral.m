%% Fuerzas Máximas de cada modo 
fs_max = zeros(cant_pisos,n_modos);
figure
hold on
for n = 1:n_modos
    for j = 1:cant_pisos
        fs_max(j,n) = Gamma_n(n)*Phi(j,n)*M(j,j)*PSa(n);
    end
    plot(fs_max(:,n),Pisos,'-o')
end
hold off
xlabel('fs_{jn}^{max}')
ylabel('Pisos')
legend(strcat('Modo ',string(1:n_modos)))
grid on

% SRSS
fs_SRSSmax = zeros(cant_pisos,1);
for j = 1:cant_pisos
    fs_SRSSmax(j,1) = sqrt(sumsqr(fs_max(j,:)));
end
figure
plot(fs_SRSSmax,Pisos,'-o')
xlabel('fs_{SRSS}^{max} [kN]')
ylabel('Pisos')
grid on
disp('fs_SRSSmax [kN]')
disp(fs_SRSSmax)

%% Máximo desplazamiento lateral de cada piso SRSS

u_jnmax = zeros(cant_pisos,n_modos);
for n = 1:n_modos
    for j = 1:cant_pisos
        u_jnmax(j,n) = Gamma_n(n)*Phi(j,n)*PSa(n)*(1/wn2(n));
    end
end

% Luego, SSRS
u_SRSSmax = zeros(cant_pisos,1);
for j = 1:cant_pisos
    u_SRSSmax(j,1) = sqrt(sumsqr(u_jnmax(j,:)));
end
figure
plot(u_SRSSmax*100,Pisos,'-o')
xlabel('u_{SRSS}^{max} [cm]')
ylabel('Pisos')
grid on
disp('u_SRSSmax [cm]')
disp(u_SRSSmax*100)

%% Máxima razón de desplazamiento de piso

rdp = zeros(cant_pisos,n_modos);
for j = 1:cant_pisos
    for n = 1:n_modos
        if j < cant_pisos
            rdp(j,n) = (1/h(j))*Gamma_n(n)*(Phi(j,n)-Phi(j+1,n))*PSa(n)*(1/wn2(n));
        elseif j == cant_pisos
            rdp(j,n) = (1/h(j))*Gamma_n(n)*(Phi(j,n))*PSa(n)*(1/wn2(n));
        end
    end
end

% Luego SRSS

rdp_SRSSmax = zeros(cant_pisos,1);
for j = 1:cant_pisos
    rdp_SRSSmax(j,1) = sqrt(sumsqr(rdp(j,:)));
end
figure
plot(rdp_SRSSmax*100,Pisos,'-o')
xlabel('rdp_{SRSS}^{max} [%]')
ylabel('Pisos')
grid on
disp('rdp_SRSSmax [%]')
disp(rdp_SRSSmax*100)

%% Máximo esfuerzo de corte de piso
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
xlabel('V_{jn} [kN]')
ylabel('Pisos')
legend(strcat('Modo ',string(1:n_modos)))
grid on

% Luego, SSRS
V_SRSSmax = zeros(cant_pisos,1);
for j = 1:cant_pisos
    V_SRSSmax(j,1) = sqrt(sumsqr(V(j,:)));
end
figure
plot(V_SRSSmax,Pisos,'-o')
xlabel('V_{SRSS}^{max} [kN]')
ylabel('Pisos')
grid on
disp('V_SRSSmax [kN]')
disp(V_SRSSmax)

%% Corte basal máximo
Vb_max = zeros(1,n_modos);
for n = 1:n_modos
    Vb_max(1,n) = sum(fs_max(:,n));
end

% SRSS
Vb_SRSSmax = sqrt(sumsqr(Vb_max(1,:)));
disp('Vb_SRSSmax [kN]')
disp(Vb_SRSSmax)

%% Momento basal máximo
Mn_max = zeros(j,n_modos);
for n = 1:n_modos
    for j = 1:cant_pisos
        Mn_max(j,n) = h_base(j)*fs_max(j,n);
    end
end

%SRSS
Mb_max = zeros(1,n_modos);
for n = 1:n_modos
    Mb_max(1,n) = sum(Mn_max(:,n));
end

% SRSS
Mb_SRSSmax = sqrt(sumsqr(Mb_max(1,:)));
disp('Mb_SRSSmax [kN-m]')
disp(Mb_SRSSmax)