
%function dy=Matrix_def9(t,y)%for shaft with mass,2discs 2bearings 9elements
clc;
format short e ;
rho=7800;
E=2.1*10^11;
m=1.5;
m_b=0.005;
rho=7800;
D=0.5;
d=0.1;
w=100;
e=0.05;
%t=0:0.6;
t0=0.6;
l=1/9;
Id=(pi/64)*(D^4);
A=pi*(d^2)/4;
k_factor=103.11;
m_factor=1.46*10^(-3);
%-------------------------------
%control parameters
k_x=7.25*10^4;
k_p=3500;
k_d=8;
k_i=29;
%----------------------------------
% Defining Mass Matrix and deleting rows and columns for reduced form
M=m_factor*[156 22*l 54 -13*l 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
   22 4*l^2 13*l -3*l^2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
   54 13*l 312 0 54 -13*l 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
   -13*l -3*l^2 0 8*l^2 13*l -3*l^2 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
   0 0 54 13*l 312 0 54 -13*l 0 0 0 0 0 0 0 0 0 0 0 0;
   0 0 -13*l -3*l^2 0 8*l^2 13*l -3*l^2 0 0 0 0 0 0 0 0 0 0 0 0;
   0 0 0 0 54 13*l 312+m/m_factor 0 54 -13*l 0 0 0 0 0 0 0 0 0 0;
   0 0 0 0 -13*l -3*l^2 0 8*l^2+Id/m_factor 13*l -3*l^2 0 0 0 0 0 0 0 0 0 0;
   0 0 0 0 0 0 54 13*l 312 0 54 -13*l 0 0 0 0 0 0 0 0;
   0 0 0 0 0 0 -13*l -3*l^2 0 8*l^2 13*l -3*l^2 0 0 0 0 0 0 0 0;
   0 0 0 0 0 0 0 0 54 13*l 312 0 54 -13*l 0 0 0 0 0 0;
   0 0 0 0 0 0 0 0 -13*l -3*l^2 0 8*l^2 13*l -3*l^2 0 0 0 0 0 0;
   0 0 0 0 0 0 0 0 0 0 54 13*l 312+m/m_factor 0 54 -13*l 0 0 0 0;
   0 0 0 0 0 0 0 0 0 0 -13*l -3*l^2 0 8*l^2+Id/m_factor 13*l -3*l^2 0 0 0 0;
   0 0 0 0 0 0 0 0 0 0 0 0 54 13*l 312 0 54 -13*l 0 0;
   0 0 0 0 0 0 0 0 0 0 0 0 -13*l -3*l^2 0 8*l^2 13*l -3*l^2 0 0;
   0 0 0 0 0 0 0 0 0 0 0 0 0 0 54 13*l 312 0 54 -13*l;
   0 0 0 0 0 0 0 0 0 0 0 0 0 0 -13*l -3*l^2 0 8*l^2 13*l -3*l^2;
   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 54 13*l 156 -22*l ;
   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -13*l -3*l^2 -22*l 4*l^2];
  M(:,1)=[];
  M(:,19)=[];
  M(1,:)=[];
  M(19,:)=[];
  
 K=k_factor*[12 6*l -12 6*l 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     6*l 4*l^2 -6*l 2*l^2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     -12 -6*l 24 0 -12 6*l 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     6*l 2*l^2 0 8*l^2 -6*l 2*l^2 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 -12 -6*l 24 0 -12 6*l 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 6*l 2*l^2 0 8*l^2 -6*l 2*l^2 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 -12 -6*l 24 0 -12 6*l 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 6*l 2*l^2 0 8*l^2 -6*l 2*l^2 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 -12 -6*l 24 0 -12 6*l 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 6*l 2*l^2 0 8*l^2 -6*l 2*l^2 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 -12 -6*l 24 0 -12 6*l 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 6*l 2*l^2 0 8*l^2 -6*l 2*l^2 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 -12 -6*l 24 0 -12 6*l 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 6*l 2*l^2 0 8*l^2 -6*l 2*l^2 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 -12 -6*l 24 0 -12 6*l 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 6*l 2*l^2 0 8*l^2 -6*l 2*l^2 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 -12 -6*l 24 0 -12 6*l;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 6*l 2*l^2 0 8*l^2 -6*l 2*l^2;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -12 -6*l 12 -6*l;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6*l 2*l^2 -6*l 4*l^2];
  K(:,1)=[];
  K(:,19)=[];
  K(1,:)=[];
  K(19,:)=[];
F_uy=zeros(18,1);
F_uy(6,1)= -1i*m_b*e*w^2;
F_uy(12,1)=-1i*m_b*e*w^2;
% F=10^-4*[0.0617;0.8;6.67*10^-3;0.9;6.67*10^-3;0.96;-2.3*10^-3;0.833;-0.0111;0.667;-0.0111;0.56;2.2*10^-3;0.733;0.015;0.966;0.015;0.565]*w^2;
% F_u=F_uy+F;
K_del=zeros(18);% [F_c]=[K_del]{n}+[K_i]{n'} governing equation
K_del(6,6)=k_x;
K_del(12,12)=k_x;
%propotional constant matrix
K_p=zeros(18);
K_p(6,6)=k_p;
K_p(12,12)=k_p;
K_d=zeros(18);
K_d(6,6)=k_d;
K_d(12,12)=k_d;
K_i=zeros(18);
K_i(6,6)=k_i;
K_i(12,12)=k_i;
%-------------------------------
O=zeros(18);
I=eye(18);
E=inv(M);
A=[O I;-E*(K-K_del) O];
B=[O;E];
C=[O;E*K_i];
F=[K_p K_d];
X=A-C*F;
EIGEN=eig(X);
%dy=X*y+B*F_uy; % Equation With amb
%dy=A*y+B*F_uy; %Equation without amb original
% A1=[O I;-E*(K) O]; 
% dy=A1*y+B*F_uy; %Equation without amb modified
K1=K-K_del+K_i*K_p;
C_1=1j*w*K_i*K_d;

for w=1:1000
      for i=1:18
         
         eta=((K-w^2*M)\F_uy);
          kas(i,w)=eta(i,1);
       eta2=(((-w^2)*M+K1+C_1)\F_uy);
       kas2(i,w)=eta2(i,1);
     end
end
Phase=angle(kas);    %============calculating phase angle without AMB
Phase1=angle(kas2);  %============calculating phase angle with AMB
%==============================================================
%==========plots displacement vs w (I AMB)=============%
% semilogy(abs(kas(6,:)),'-') 
% hold on; grid off;
% semilogy(abs(kas2(6,:)),'-r')
%  xlabel('\omega (rad/sec)','fontsize',12)
% ylabel('Displacement Y in m')
% title('Plot of Disp vs \omega for 2 amb and 2 disc(I amb)[9 elements][distributed eccentricity]','FontSize',12)
% legend('Without AMB','With AMB')
%==========plots displacement vs w (II AMB)=============%
% semilogy(abs(kas(12,:)),'-') 
% hold on; grid off;
% semilogy(abs(kas2(12,:)),'-r')
%  xlabel('\omega (rad/sec)','fontsize',12)
% ylabel('Displacement Y in m')
% title('Plot of Disp vs \omega for 2 amb and 2 disc(II amb)[9 elements][distributed eccentricity]','FontSize',12)
% legend('Without AMB','With AMB')
% ==========plots Phase vs w (I AMB)====================%
plot(Phase(6,:),'-') 
hold on; grid off;
plot(Phase1(6,:),'-r')
 xlabel('\omega (rad/sec)','fontsize',14)
ylabel('Phase \phi')
title('Plot of Phase \phi vs \omega for 2 amb and 2 disc(I amb)[9 elements][distributed eccentricity]','FontSize',12)
legend('Without AMB','With AMB')
% ==========plots Phase vs w (II AMB)====================%
% plot(Phase(12,:),'-') 
% hold on; grid off;
% plot(Phase1(12,:),'-r')
%  xlabel('\omega (rad/sec)','fontsize',12)
% ylabel('Phase \phi')
% title('Plot of Phase \phi vs \omega for 2 amb and 2 disc(II amb)[9 elements][distributed eccentricity]','FontSize',12)
% legend('Without AMB','With AMB')
