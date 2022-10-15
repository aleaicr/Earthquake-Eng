%% Fuerzas Máximas de cada modo 
fs_max = zeros(cant_pisos,n_modos);
for n = 1:n_modos
    for j = 1:cant_pisos
        fs_max(j,n) = Gamma_n(n)*Phi(j,n)*M(j,j)*PSa(n);
    end
end

% SRSS
fs_SRSSmax = zeros(cant_pisos,1);
for j = 1:cant_pisos
    fs_SRSSmax(j,1) = sqrt(sumsqr(fs_max(j,:)));
end


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

%% Máximo esfuerzo de corte de piso
% Primero obtenemos el corte de piso para cada modo

V = zeros(cant_pisos,n_modos);
for n = 1:n_modos
    for j = 1:cant_pisos
        V(j,n) = sum(fs_max(1:j,n));    %tonf
    end
end

% Luego, SSRS
V_SRSSmax = zeros(cant_pisos,1);
for j = 1:cant_pisos
    V_SRSSmax(j,1) = sqrt(sumsqr(V(j,:)));
end

%% Corte basal máximo
Vb_max = zeros(1,n_modos);
for n = 1:n_modos
    Vb_max(1,n) = sum(fs_max(:,n));
end

% SRSS
Vb_SRSSmax = sqrt(sumsqr(Vb_max(1,:)));

%% Momento basal máximo
Mn_max = zeros(j,n_modos);
for n = 1:n_modos
    for j = 1:cant_pisos
        Mn_max(j,n) = h_base(j)*fs_max(j,n);
    end
end

Mb_max = zeros(1,n_modos);
for n = 1:n_modos
    Mb_max(1,n) = sum(Mn_max(:,n));
end

% SRSS
Mb_SRSSmax = sqrt(sumsqr(Mb_max(1,:)));