function [Cy_min] = Cy_vs_mu(beta,T,xi,dt,u_i,ud_i,Datos,Error,t,mu_max)
% cont = 1;
k = 4*pi^2/T^2;
Cy = transpose(0.005:0.001:0.5);                                            % Barrido de Cy
W = 9.81;
mu = zeros(length(Cy),1);

for i=1:length(Cy)
    Fy = Cy(i,1)*W;
    uy = Fy/k;
    [u,ud,udd,fs,fluencia] = Newmark_No_Lineal(beta,T,xi,dt,u_i,ud_i,Datos,Fy,Error,t); % Para coef.sismico inelástico, es No Lineal
    mu(i,1) = max(abs(u))/uy;
end

plot(mu,Cy)

mu = flip(mu,1);
Cy = flip(Cy,1);
c = 0;

for i = 1:length(mu)
    if c == 1
        c = 1;
    elseif mu(i,1) > mu_max
        Cy_min = Cy(i-1,1);
        c = 1;
    elseif mu(i,1) == mu_max
        Cy_min = Cy(i,1);
        c = 1;
    end
end
end