 function dy=Matrix_def3_withoutamb(t,y) % for shaft with mass, one disc 2 bearings

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
%=======================================
%control parameters
k_x=7.25*10^4;
k_p=3500;
k_d=11;
k_i=29;
%=======================================
% Defining Mass Matrix and deleting rows and columns for reduced form
M=m_factor*[4*l^2 13*l -3*l^2 0 0 0 0 0;13*l 312 0 54 -13*l 0 0 0;-3*l^2 0 8*l^2 13*l -3*l^2 0 0 0;0 54 13*l (312+m/m_factor) 0 54 -13*l 0;0 -13*l -3*l^2 0 (8*l^2+Id/m_factor) 13 -3*l^2 0;0 0 0 54 13*l 312 0 -13*l;0 0 0 -13*l -3*l^2 0 8*l^2 -3*l^2;0 0 0 0 0 -13*l -3*l^2 4*l^2];
  
 K=k_factor*[4*l^2 -6*l 2*l^2 0 0 0 0 0;-6*l 24 0 -12 -6*l 0 0 0;2*l^2 0 8*l^2 -6*l 2*l^2 0 0 0;0 -12 -6*l 24 0 -12 6*l 0;0 6*l 2*l^2 0 8*l^2 -6*l 2*l^2 0;0 0 0 -12 -6*l 24 0 6*l;0 0 0 6*l 2*l^2 0 8*l^2 2*l^2;0 0 0 0 0 6*l 2*l^2 4*l^2];
    %=================Impulse Function @ t=0.6 sec=================
%    function h = dd1(t)
% % Our default value is 0
% h = 0; 
% t0=0.6;
% 
% % The function is 1 only if the input is 0
% if t == t0
%     h = 1;
% end
%    end 
%==================================================================
F_uy=zeros (8,1);
F_uy(4,1)=10*dd1(t);
F_uy(4,1)=m_b*e*w^2*sin(w*t);
K_del=zeros(8);% [F_c]=[K_del]{n}+[K_i]{n'} governing equation
K_del(2,2)=k_x;
K_del(6,6)=k_x;
%propotional constant matrix
K_p=zeros(8);
K_p(2,2)=k_p;
K_p(6,6)=k_p;
K_d=zeros(8);
K_d(2,2)=k_d;
K_d(6,6)=k_d;
K_i=zeros(8);
K_i(2,2)=k_i;
K_i(6,6)=k_i;
%-------------------------------
O=zeros(8);
I=eye(8);
E=inv(M);
A=[O I;-E*(K-K_del) O];
B=[O;E];
C=[O;M\K_i];
F=[K_p K_d];
X=A-C*F;
%dy=X*y+B*F_uy; % Equation With amb
A1=[O I;-E*(K) O]; 
dy=A1*y+B*F_uy; %Equation without amb modified
 end
