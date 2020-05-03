%%
% LPV Modeling of Boost Converter and Gain Scheduling MPC Control
% Author: Rosana C B Rego
% 2019 IEEE 15th Brazilian Power Electronics Conference
%and 5th IEEE Southern Power Electronics Conference (COBEP/SPEC)
%%

%LMI Optimization 

clear; clc; close all
%% Model parameters
Pot=[1000 380]; Vo=48; Vg=[36 26]; L=36e-6; Co=4400e-6; Rco=26.7e-3;
Ts=1e-3;
Ro = (Vo^2)/Pot(2);
Romax = (Vo^2)/Pot(1);
Dcycle_min = 1 - Vg(2)/Vo; %duty cicle
Dcycle_max = 1 - Vg(1)/Vo; %duty cicle
%% Polytopics
alf1=[(1-Dcycle_min) (1-Dcycle_max)];
alf2=[((1 - Dcycle_min)^2)*Ro + Dcycle_min*(1 - Dcycle_min)*((Rco*Ro)/(Rco+Ro)) ((1 - Dcycle_max)^2)*Romax + Dcycle_max*(1 - Dcycle_max)*((Rco*Romax)/(Rco+Romax)) ].^(-1);
alf3=[((1 - Dcycle_min))*Ro + Dcycle_min*((Rco*Ro)/(Rco+Ro)) ((1 - Dcycle_max))*Romax + Dcycle_max*((Rco*Romax)/(Rco+Romax)) ].^(-1);
%Weighting matrix
Le=diag([1 1 1],0); R=1;
g=1;h=1;
%Constrain
umax = 0.5;
%Initial states
%% x_set
xset=[(Pot(1)/Vg(1)) Vo]';
N = 10;
%alpha=alf1(2)*rand(1,N) +alf1(1);
[A1, B1 , C1, D1] = BOOST(Vg(1),Rco,1-Vg(1)/Vo,L,(Vo^2)/Pot(1),Co); %Modelo
[Ad1,Bd1]=c2d(A1,B1,Ts);
for k=1:N
    xset(:,k+1)=Ad1*xset(:,k);
end
%% LPV variables
A_0=[];A_1=[];B_2=[];B_3=[];C_0=[];C_1=[]; A=[]; B=[]; C=[]; D=[];
alpha0=sdpvar(1);
alpha1=sdpvar(1);
alpha2=sdpvar(1);
alpha3=sdpvar(1);
A_0 = [-(((Rco*Ro)/(Rco+Ro)))/L -(Ro)/((Rco +Ro)*L);
    (Ro)/(Co*(Rco+Ro)) -1/(Co*(Rco+Ro))];
[A_0,B_0]=c2d(A_0, [0;0],Ts);
A_1 = [-(((Rco*Ro)/(Rco+Ro)))/L -(Ro)/((Rco +Ro)*L); (Ro)/(Co*(Rco+Ro)) -1/(Co*(Rco+Ro))];
B_2 = [(Ro/L)*((Rco)/(Ro+Rco)); -Ro/((Ro+Rco))];
[A_, B_2]=c2d([0 0; 0 0], B_2,Ts);
[A_1, B_]=c2d(A_1,[0;0],Ts);
B_3 = [(Ro/L)*((Ro)/(Ro+Rco)); 0];
[A_,B_3]=c2d([0 0; 0 0],B_3,Ts);
B_3=alpha3*B_3;
B_2=alpha2*B_2;
B=B_2+B_3;
B=Vg(1)*B;
A_1=alpha1*A_1;
A=A_0+A_1;
C_0=[0 Ro/(Rco+Ro)];
C_1=[alpha1*((Rco*Ro)/(Rco+Ro))  0];
C=C_0+C_1;
D_0=1;
D_2= -Vg(1)*((alpha2*Rco));
D=D_0+D_2;
% Extend model
g=1;h=1;
Bhat1=[B; -h*D];
n=size(Bhat1,1); m=size(Bhat1,2);
Ahat1=[A zeros(2,1); -h*C g];
A_0=[];A_1=[];B_2=[];B_3=[];C_0=[];C_1=[]; A=[]; B=[]; C=[]; D=[];
A_0 = [-(((Rco*Romax)/(Rco+Romax)))/L -(Romax)/((Rco +Romax)*L);
    (Romax)/(Co*(Rco+Romax)) -1/(Co*(Rco+Romax))];
[A_0,B_0]=c2d(A_0, [0;0],Ts);
A_1 = [-(((Rco*Romax)/(Rco+Romax)))/L -(Romax)/((Rco +Romax)*L); (Romax)/(Co*(Rco+Romax)) -1/(Co*(Rco+Romax))];
B_2 = [(Romax/L)*((Rco)/(Romax+Rco)); -Romax/((Romax+Rco))];
[A_, B_2]=c2d([0 0; 0 0], B_2,Ts);
[A_1, B_]=c2d(A_1,[0;0],Ts);
B_3 = [(Romax/L)*((Romax)/(Romax+Rco)); 0];
[A_,B_3]=c2d([0 0; 0 0],B_3,Ts);
B_3=alpha3*B_3;
B_2=alpha2*B_2;
B=B_2+B_3;
B=Vg(1)*B;
A_1=alpha1*A_1;
A=A_0+A_1;
C_0=[0 Romax/(Rco+Romax)];
C_1=[alpha1*((Rco*Romax)/(Rco+Romax))  0];
C=C_0+C_1;
D_0=1;
D_2= -Vg(1)*((alpha2*Rco));
D=D_0+D_2;;
% Extend model
g=1;h=1;
Bhat2=[B; -h*D];
n=size(Bhat1,1); m=size(Bhat1,2);
Ahat2=[A zeros(2,1); -h*C g];
%%
A_0=[];A_1=[];B_2=[];B_3=[];C_0=[];C_1=[]; A=[]; B=[]; C=[]; D=[];
A_0 = [-(((Rco*Romax)/(Rco+Romax)))/L -(Romax)/((Rco +Romax)*L);
    (Romax)/(Co*(Rco+Romax)) -1/(Co*(Rco+Romax))];
[A_0,B_0]=c2d(A_0, [0;0],Ts);
A_1 = [-(((Rco*Romax)/(Rco+Romax)))/L -(Romax)/((Rco +Romax)*L); (Romax)/(Co*(Rco+Romax)) -1/(Co*(Rco+Romax))];
B_2 = [(Romax/L)*((Rco)/(Romax+Rco)); -Romax/((Romax+Rco))];
[A_, B_2]=c2d([0 0; 0 0], B_2,Ts);
[A_1, B_]=c2d(A_1,[0;0],Ts);
B_3 = [(Romax/L)*((Romax)/(Romax+Rco)); 0];
[A_,B_3]=c2d([0 0; 0 0],B_3,Ts);
B_3=alpha3*B_3;
B_2=alpha2*B_2;
B=B_2+B_3;
B=Vg(2)*B;
A_1=alpha1*A_1;
A=A_0+A_1;
C_0=[0 Romax/(Rco+Romax)];
C_1=[alpha1*((Rco*Romax)/(Rco+Romax))  0];
C=C_0+C_1;
D_0=1;
D_2= -Vg(2)*((alpha2*Rco));
D=D_0+D_2;;
% Extend model
g=1;h=1;
Bhat3=[B; -h*D];
n=size(Bhat1,1); m=size(Bhat1,2);
Ahat3=[A zeros(2,1); -h*C g];
%%
A_0=[];A_1=[];B_2=[];B_3=[];C_0=[];C_1=[]; A=[]; B=[]; C=[]; D=[];
A_0 = [-(((Rco*Ro)/(Rco+Ro)))/L -(Ro)/((Rco +Ro)*L);
    (Ro)/(Co*(Rco+Ro)) -1/(Co*(Rco+Ro))];
[A_0,B_0]=c2d(A_0, [0;0],Ts);
A_1 = [-(((Rco*Ro)/(Rco+Ro)))/L -(Ro)/((Rco +Ro)*L); (Ro)/(Co*(Rco+Ro)) -1/(Co*(Rco+Ro))];
B_2 = [(Ro/L)*((Rco)/(Ro+Rco)); -Ro/((Ro+Rco))];
[A_, B_2]=c2d([0 0; 0 0], B_2,Ts);
[A_1, B_]=c2d(A_1,[0;0],Ts);
B_3 = [(Ro/L)*((Ro)/(Ro+Rco)); 0];
[A_,B_3]=c2d([0 0; 0 0],B_3,Ts);
B_3=alpha3*B_3;
B_2=alpha2*B_2;
B=B_2+B_3;
B=Vg(2)*B;
A_1=alpha1*A_1;
A=A_0+A_1;
C_0=[0 Ro/(Rco+Ro)];
C_1=[alpha1*((Rco*Ro)/(Rco+Ro))  0];
C=C_0+C_1;
D_0=1;
D_2= -Vg(2)*((alpha2*Rco));
D=D_0+D_2;
% Extend model
g=1;h=1;
Bhat4=[B; -h*D];
n=size(Bhat1,1); m=size(Bhat1,2);
Ahat4=[A zeros(2,1); -h*C g];
%%  OFF-line MPC
%LMI variables
Q = sdpvar(n,n, 'symmetric');
gamma = sdpvar(1,1);
Y0 = sdpvar(m,n, 'full');
Y1 = sdpvar(m,n, 'full');
Y2 = sdpvar(m,n, 'full');
Y3 = sdpvar(m,n, 'full');
Q0 = sdpvar(n,n, 'symmetric');
Q1 = sdpvar(n,n, 'symmetric');
Q2 = sdpvar(n,n, 'symmetric');
Q3 = sdpvar(n,n, 'symmetric');
X0 = sdpvar(m,m, 'full');
X1 = sdpvar(m,m, 'full');
X2 = sdpvar(m,m, 'full');
X3 = sdpvar(m,m, 'full');
xp = sdpvar(n,1);
% constraints and objective
objective = gamma;
%optimization object
ops = sdpsettings('solver','sedumi','sedumi.eps',1e-5);
% Optimization with LPV variables
LA_ = [Q0 Q0*Ahat1'+Y0'*Bhat1' Q0*sqrtm(Le) Y0'*sqrtm(R);
    Ahat1*Q0+Bhat1*Y0 Q0 zeros(n,n) zeros(n,m);
    sqrtm(Le)*Q0 zeros(n,n) gamma*eye(n) zeros(n,m);
    sqrtm(R)*Y0 zeros(m,n) zeros(m,n) gamma*eye(m)];
LA_2 = [Q1 Q1*Ahat2'+Y1'*Bhat2' Q1*sqrtm(Le) Y1'*sqrtm(R);
    Ahat2*Q1+Bhat2*Y1 Q1 zeros(n,n) zeros(n,m);
    sqrtm(Le)*Q1 zeros(n,n) gamma*eye(n) zeros(n,m);
    sqrtm(R)*Y1 zeros(m,n) zeros(m,n) gamma*eye(m)];
LA_3 = [Q2 Q2*Ahat3'+Y2'*Bhat3' Q2*sqrtm(Le) Y2'*sqrtm(R);
    Ahat3*Q2+Bhat3*Y2 Q2 zeros(n,n) zeros(n,m);
    sqrtm(Le)*Q2 zeros(n,n) gamma*eye(n) zeros(n,m);
    sqrtm(R)*Y2 zeros(m,n) zeros(m,n) gamma*eye(m)];
LA_4 = [Q3 Q3*Ahat4'+Y3'*Bhat4' Q3*sqrtm(Le) Y3'*sqrtm(R);
    Ahat4*Q3+Bhat4*Y3 Q3 zeros(n,n) zeros(n,m);
    sqrtm(Le)*Q3 zeros(n,n) gamma*eye(n) zeros(n,m);
    sqrtm(R)*Y3 zeros(m,n) zeros(m,n) gamma*eye(m)];
L3 = [[X0 Y0; Y0' Q0]>=0, X0<=umax.^2];
L31 = [[X1 Y1; Y1' Q1]>=0, X1<=umax.^2];
L32 = [[X2 Y2; Y2' Q2]>=0, X2<=umax.^2];
L33 = [[X3 Y3; Y3' Q3]>=0, X3<=umax.^2];
Qn_=1*eye(n);
Y=Y0+alpha1*Y1+alpha2*Y2+alpha3*Y3;
Q=Q0+alpha1*Q1+alpha2*Q2+alpha3*Q3;
for i = 1:1
    LMIs = [LA_ >= 0, LA_2 >= 0, LA_3 >= 0, LA_4 >= 0, L3,L31,L32,L33];
    LMIs = [LMIs, 0<=alpha0<=0.1, uncertain(alpha0)];
    LMIs = [LMIs, alf1(1)<=alpha1<=alf1(2), uncertain(alpha1)];
    LMIs = [LMIs, alf2(1)<=alpha2<=alf2(2), uncertain(alpha2)];
    LMIs = [LMIs, alf3(1)<=alpha3<=alf3(2), uncertain(alpha3)];
    LMIs = [LMIs, Q0>=0, Q1>=0,Q2>=0,Q3>=0];
    controller = optimizer(LMIs,objective,ops,xp,{Y0,Y1,Y2,Y3,gamma,Q0, Q1, Q2, Q3,alpha1, alpha2, alpha3, Y,Q,alpha0});
    x=[xset(1,k); xset(2,k);0];
    sol= controller{x};
    Yn0(:,:,i) = sol{1};
    Yn1(:,:,i) = sol{2};
    Yn2(:,:,i) = sol{3};
    Yn3(:,:,i) = sol{4};
    gammav(i) = sol{5};
    Qn(:,:,i) = sol{6};
    Qn1(:,:,i) = sol{7};
    Qn2(:,:,i) = sol{8};
    Qn3(:,:,i) = sol{9};
end

N_=1;
F0 = Yn0(:,:,N_)*inv(Qn(:,:,N_));
F1 = Yn1(:,:,N_)*inv(Qn1(:,:,N_));
F2 = Yn2(:,:,N_)*inv(Qn2(:,:,N_));
F3 = Yn3(:,:,N_)*inv(Qn3(:,:,N_));
