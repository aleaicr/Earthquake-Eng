clc
clear all
%% Datos
alfa=-0.15;
xi=0.05;
Tn=1;
R=10^(-10);
beta=1/6;
u_i=0;
ud_i=0;
u_colapso=100;
%% Concepción
Archivo='Concepcion2010-L.txt';
delimeterIn=' ';
headerlinesIn=1;
Registro=importdata(Archivo,delimeterIn,headerlinesIn);
Datos_C=Registro.data;
dt_C=0.005;
T_C=transpose(0:dt_C:(size(Datos_C,1)-1)*dt_C);
%% Santiago
Archivo='stgocentro2010-L.txt';
delimeterIn=' ';
headerlinesIn=1;
Registro=importdata(Archivo,delimeterIn,headerlinesIn);
Datos_S=Registro.data;
dt_S=0.005;
T_S=transpose(0:dt_S:(size(Datos_S,1)-1)*dt_S);
%% Talca
Archivo='Talca2010-L.txt';
delimeterIn=' ';
headerlinesIn=1;
Registro=importdata(Archivo,delimeterIn,headerlinesIn);
Datos_T=Registro.data;
dt_T=0.005;
T_T=transpose(0:dt_T:(size(Datos_T,1)-1)*dt_T);
%% Valparaíso
Archivo='valparaisoUTFSM2010-L.txt';
delimeterIn=' ';
headerlinesIn=1;
Registro=importdata(Archivo,delimeterIn,headerlinesIn);
Datos_V=Registro.data;
dt_V=0.005;
T_V=transpose(0:dt_V:(size(Datos_V,1)-1)*dt_V);
%% Parte i
Cy=0:0.01:0.7;
n=size(Cy,2);
for i=1:n
    [u_C,ud_C,udd_C,fs,fluencia,ductilidad]=Newmark_No_Lineal_Bilineal(alfa,beta,Tn,xi,dt_C,u_i,ud_i,Datos_C,Cy(1,n-i+1)*980,R,T_C);
    if max(abs(u_C))>u_colapso
        Cy_conce=Cy(1,n-i+2)
        break
    end
end
for i=1:n
    [u_S,ud_S,udd_S,fs,fluencia,ductilidad]=Newmark_No_Lineal_Bilineal(alfa,beta,Tn,xi,dt_S,u_i,ud_i,Datos_S,Cy(1,n-i+1)*980,R,T_S);
    if max(abs(u_S))>u_colapso
        Cy_stgo=Cy(1,n-i+2)
        break
    end
end
for i=1:n
    [u_T,ud_T,udd_T,fs,fluencia,ductilidad]=Newmark_No_Lineal_Bilineal(alfa,beta,Tn,xi,dt_T,u_i,ud_i,Datos_T,Cy(1,n-i+1)*980,R,T_T);
    if max(abs(u_T))>u_colapso
        Cy_talca=Cy(1,n-i+2)
        break
    end
end
for i=1:n
    [u_V,ud_V,udd_V,fs,fluencia,ductilidad]=Newmark_No_Lineal_Bilineal(alfa,beta,Tn,xi,dt_V,u_i,ud_i,Datos_V,Cy(1,n-i+1)*980,R,T_V);
    if max(abs(u_V))>u_colapso
        Cy_valpo=Cy(1,n-i+2)
        break
    end
end
%% Parte ii
[Sd,Sv,Sa_C,PSv,PSa]=Newmark_Lineal(beta,xi,dt_C,u_i,ud_i,Datos_C);
[Sd,Sv,Sa_S,PSv,PSa]=Newmark_Lineal(beta,xi,dt_S,u_i,ud_i,Datos_S);
[Sd,Sv,Sa_T,PSv,PSa]=Newmark_Lineal(beta,xi,dt_T,u_i,ud_i,Datos_T);
[Sd,Sv,Sa_V,PSv,PSa]=Newmark_Lineal(beta,xi,dt_V,u_i,ud_i,Datos_V);

for i=1:size(Sa_C,1)
    if Sa_C(i,1)==1
        Sa_conce=Sa_C(i,2);
        Sa_stgo=Sa_S(i,2);
        Sa_talca=Sa_T(i,2);
        Sa_valpo=Sa_V(i,2);
    end
end
Ce_C=Sa_conce/980
Rc_C=Ce_C/Cy_conce
Ce_S=Sa_stgo/980
Rc_S=Ce_S/Cy_stgo
Ce_T=Sa_talca/980
Rc_T=Ce_T/Cy_talca
Ce_V=Sa_valpo/980
Rc_V=Ce_V/Cy_valpo
%% Parte iii
a=0.26*(1-exp(-7.5*Tn));
b=0.89+0.04*Tn+0.15*log(Tn);
Rc=1+a*(-alfa)^(-b)