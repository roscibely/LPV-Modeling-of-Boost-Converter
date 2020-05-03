%% Simulation
Ts=1e-3;
t=0:Ts:500e-3;
x=[Pot(2)/Vg(2); Vg(2)]; y=0; u=0; v=0;
N=length(t);
alpha=(alf1(2)*rand(1,N))+alf1(1);
alpha2=((alf2(2))*rand(1,N))+alf2(1);
alpha3=((alf3(2))*rand(1,N))+alf3(1);
v=0; eh=0; ud=0;yd=0; u_d=0;y_d=0;
R=[Vo*ones(1,ceil(length(t)*0.5)) Vo*ones(1,ceil(length(t)*0.5))];
Ppot=[Pot(2) Pot(2)*ones(1,ceil((length(t)-1)*0.25)) Pot(2)*ones(1,ceil((length(t)-1)*0.25)) Pot(1)*ones(1,ceil((length(t)-1)*0.25)) Pot(1)*ones(1,ceil((length(t)-1)*0.25))];
Vvg=[Vg(1) Vg(1)*ones(1,ceil((length(t)-1)*0.25)) Vg(2)*ones(1,ceil((length(t)-1)*0.25)) Vg(2)*ones(1,ceil((length(t)-1)*0.25)) Vg(1)*ones(1,ceil((length(t)-1)*0.25))];
%
for i=1:N
    Kmpc = [0.000614287326110045 -0.00222411667432920 -0.000775919788174078];
    u(i+1)=(-Kmpc(1:end-1)*x(:,i)-Kmpc(end)*v)-ud;
    u_ant=u(i+1);
    if u(i+1)<0, u(i+1)=0;end
    if u(i+1)>umax, u(i+1)=umax;end
    u_dep=u(i+1);
    [A_, B_, C_]=boost_n_linear(Ppot(i),u(i),Vo,L,Co,Rco);
    % Saturação da corrente no indutor
    if x(1,i)<0, x(1,i)=0;end
    if x(1,i)>105, x(1,i)=105;end
    % Saturação da tensão no Capacitor
    if x(2,i)<0, x(2,i)=0;end
    if x(2,i)>63, x(2,i)=63;end
    % Modelo RK4
    k1=A_*x(:,i)+B_*Vvg(i);
    k2=A_*(x(:,i)+0.5*Ts*k1)+B_*Vvg(i);
    k3=A_*(x(:,i)+0.5*Ts*k2)+B_*Vvg(i);
    k4=A_*(x(:,i)+Ts*k3)+B_*Vvg(i);
    x(:,i+1)=x(:,i)+(Ts/6)*(k1+2*(k2+k3)+k4);
    % Equação da saida
    y(i)=C_*x(:,i);
    if y(i)<0, y(i)=0;end
    y_in=y(i)+yd;
    % Integrador
    v=g*v+h*(R(i)-y_in);
    Vo=y(i);
    y_d(i)=yd;
    u_d(i)=ud;
end
% %
figure(2)
plot(t,x(1,1:end-1),'k','linewidth',2), legend('i_L'), grid
axis([t(1) t(end) 0 105])
set(gca,'fontsize',16,'fontname','Times New Roman')
xlabel('Time (s)','fontsize',16,'fontname','Times New Roman','fontangle','normal'),
ylabel('Current I_L (A)','fontsize',16,'fontname','Times New Roman','fontangle','normal')

figure(3)
plot(t,y,'k',t,R(1:end-1),'r-.','linewidth',2), legend('y(k)','Ref'), grid
axis([t(1) t(end) 25 65])
set(gca,'fontsize',20,'fontname','Times New Roman')
xlabel('Time (s)','fontsize',25,'fontname','Times New Roman','fontangle','normal'),
ylabel('Output Voltage (V)','fontsize',25,'fontname','Times New Roman','fontangle','normal')

figure(4)
plot(t,u(1:end-1),'k','linewidth',2), legend('u(k)'), grid
axis([t(1) t(end) .0 0.55])
set(gca,'fontsize',20,'fontname','Times New Roman')
xlabel('Time (s)','fontsize',25,'fontname','Times New Roman','fontangle','normal'),
ylabel('Duty Cycle','fontsize',25,'fontname','Times New Roman','fontangle','normal')
figure(5)
plot(t,alpha,'r','linewidth',1.5), grid
set(gca,'fontsize',20,'fontname','Times New Roman')
xlabel('Time (s)','fontsize',25,'fontname','Times New Roman','fontangle','normal'),
ylabel('\alpha_1','fontsize',25,'fontname','Times New Roman','fontangle','normal')
figure(6)
plot(t,alpha2,'r','linewidth',1.5), grid
set(gca,'fontsize',20,'fontname','Times New Roman')
xlabel('Time (s)','fontsize',25,'fontname','Times New Roman','fontangle','normal'),
ylabel('\alpha_2','fontsize',25,'fontname','Times New Roman','fontangle','normal')
figure(7)
plot(t,alpha3,'r','linewidth',1.5), grid
set(gca,'fontsize',20,'fontname','Times New Roman')
xlabel('Time (s)','fontsize',25,'fontname','Times New Roman','fontangle','normal'),
ylabel('\alpha_3','fontsize',25,'fontname','Times New Roman','fontangle','normal')
figure(8)
plot(t,Ppot,'k','linewidth',2), legend('Pot'), grid
set(gca,'fontsize',16,'fontname','Times New Roman')
xlabel('Time (s)','fontsize',16,'fontname','Times New Roman','fontangle','normal'),
ylabel('Power (W)','fontsize',16,'fontname','Times New Roman','fontangle','normal')
figure(9)
plot(t,Vvg,'k','linewidth',2), legend('V_g'), grid
set(gca,'fontsize',16,'fontname','Times New Roman')
xlabel('Time (s)','fontsize',16,'fontname','Times New Roman','fontangle','normal'),
ylabel('Input voltage (V)','fontsize',16,'fontname','Times New Roman','fontangle','normal')
