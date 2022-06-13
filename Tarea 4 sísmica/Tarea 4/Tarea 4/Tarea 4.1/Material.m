clc
clear all

Tn=2;
Fy=0.08*980;
alfa=-0.15;

omega_n=2*pi/Tn;
k=omega_n^2;
u_y=Fy/k;
u=0:0.0001:20;
fs=zeros(size(u,2),1);
for i=1:size(u,2)
    fs_EPP(i,1)=k*u(1,i)*(1-alfa);
    fs_lineal(i,1)=u(1,i)*alfa*k;
    fs(i,1)=k*u(1,i);
    if fs(i,1)>=Fy
        fs_EPP(i,1)=Fy*(1-alfa);
        fs(i,1)=fs_EPP(i,1)+fs_lineal(i,1);
    end
end
hold on
plot(u,fs)
plot(u,fs_EPP)
plot(u,fs_lineal)
xlabel('Desplazamiento [cm]','FontSize',18)
ylabel('Fuerza [N/kg]','FontSize',18)

txt = '\leftarrow k';
text(u_y/2,u_y/2*k,txt,'FontSize',18)
txt = '\leftarrow \alpha k';
text(u_y*1.75+0.5,Fy*(1-alfa)+u_y*1.75*alfa*k+1.5,txt,'FontSize',18)

txt = '(1-\alpha)k \rightarrow';
text(u_y/2-1.75,u_y/2*k*(1-alfa),txt,'FontSize',18)
txt = '\leftarrow \alpha k';
text(10,10*alfa*k+1.5,txt,'FontSize',18)

legend('Bilineal','EPP','Lineal','FontSize',18)