
function dy=Matrix_def_withoutamb(t,y)%without amb one disc one amb!!!!
clc;
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
l=1;
Id=(pi/64)*(D^4);
A=pi*(d^2)/4;
k_factor=103.11;
m_factor=1.46*10^(-3);
% m_factor=rho*A*l/420;
% k_factor=E*Id;
%-------------------------------
%control parameters
k_x=7.25*10^4;
k_p=3500;
k_d=11;
k_i=29;
%----------------------------------
% Defining Mass Matrix and deleting rows and columns for reduced form
% M=diag([10^(-5) 10^(-5) 10^(-5) m Id 10^(-5)],0);

M=m_factor*[156 22*l 54 -3*l 0 0 0 0;
           22 4*l^2 13*l -3*l^2 0 0 0 0;
   54 13*l^2 312 0 54 -13*l^2 0 0;
   -13*l -3*l^2 0 8*l^2 13*l -3*l^2 0 0;
   0 0 54 13*l (312+m/m_factor) 0 54 -13*l;
   0 0 -13*l -3*l^2 0 (8*l^2+Id/m_factor) 13 -13*l;
   0 0 0 0 54 13*l^2 152 -22*l;
   0 0 0 0 -13*l -3*l^2 -22*l 4*l^2];
  M(:,1)=[];
  M(:,7)=[];
  M(1,:)=[];
  M(7,:)=[];
  
 K=k_factor*[ 12 6*l -12 6*l 0 0 0 0;
     6*l 4*l^2 -6*l 2*l^2 0 0 0 0;
     -12 -6*l 24 0 -12 -6*l 0 0;
     6*l 2*l^2 0 8*l^2 -6*l 2*l^2 0 0 ;
     0 0 -12 -6*l 24 0 -12 6*l;
     0 0 -6*l 2*l^2 0 8*l^2 -6*l 2*l^2;
     0 0 0 0 -12 -6*l 12 -6*l;
     0 0 0 0 6*l 2*l^2 -6*l 4*l^2];
  K(:,1)=[];
  K(:,7)=[];
  K(1,:)=[];
  K(7,:)=[];
  %=================Impulse Function @ t=0.6 sec=================
%    function y = dd1(t)
% % Our default value is 0
% y = 0; 
% t0=1;
% 
% % The function is 1 only if the input is 0
% if t == t0
%     y = 1;
% end
%    end 
%==================================================================
F_uy=zeros (6,1);
F_uy(4,1)= m_b*e*w^2*sin(w*t);
% F_uy(2,1)=50*dd1(t);
K_del=zeros(6);% [F_c]=[K_del]{n}+[K_i]{n'} governing equation
K_del(2,2)=k_x;
%propotional constant matrix
K_p=zeros(6);
K_p(2,2)=k_p;
K_d=zeros(6);
K_d(2,2)=k_d;
K_i=zeros(6);
K_i(2,2)=k_i;
%-------------------------------
O=zeros(6);
I=eye(6);
E=inv(M);
A=[O I;-E*(K+K_del) O];
A1=[O I;-E*(K) O];
B=[O;E];
C=[O;M\K_i];
F=[K_p K_d];
X=A-C*F;
dy=A*y+B*F_uy; %Equation without amb
   end 





