clc
clear all
%% Datos
load('el_centro.mat');
load('resultados.mat');
Fs=fs;
alfa=-0.15;
xi=0.05;
Fy=0.08;
Tn=2;
R=10^(-10);
uddg=uppG;
T=transpose(0:dt:(size(uppG,1)-1)*dt);
beta=1/6;
u_i=0;
ud_i=0;
%% Parte i
[u,ud,udd,fs,fluencia,ductilidad]=Newmark_No_Lineal_Bilineal(alfa,beta,Tn,xi,dt,u_i,ud_i,uddg,Fy,R,T);
u_res=(uNL);
figure(1)
hold on
plot(T,u*980)
plot(t,u_res)
A=max(abs(u))*980
B=max(abs(u_res))
(A-B)/B
xlabel('Tiempo [s]','Fontsize',18)
ylabel('Desplazamiento [cm]','Fontsize',18)
legend("Aproximación Numérica","Real",'Fontsize',18)
hold off
%% Parte ii
figure(2)
hold on
plot(u*980,fs*980)
plot(u_res,Fs)
xlabel('Desplazamiento [cm]','Fontsize',18)
ylabel('Fuerza Restitutiva','Fontsize',18)
legend("Aproximación Numérica","Real",'Fontsize',13)
hold off