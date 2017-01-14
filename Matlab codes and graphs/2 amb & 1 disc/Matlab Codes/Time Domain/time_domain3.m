%Solving the diffrential eqaution using ode23 %=============2 Amb 1 Disc
%options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5]);
y0=zeros(16,1);
[t,y] = ode23(@Matrix_def3,[0 5],y0);
[t1,y1] = ode23(@Matrix_def3_withoutamb,[0 5],y0);
%======================= I AMB===========%
% plot(t,y(:,2),'-',t1,y1(:,2),'-')
% xlabel('time (Sec)')
% ylabel('displacement Y (in m)')
% title('Plot of Displacement vs time (for 2 amb and 1 disc) 1st AMB','FontSize',10)
% legend('With AMB','Without AMB')
%============================ II AMB=====%

plot(t,y(:,6),'-',t1,y1(:,6),'-')
xlabel('time (Sec)')
ylabel('displacement Y (in m)')
title('Plot of Displacement vs time ( 2 amb and 1 disc) 2nd AMB','FontSize',10)
legend('With AMB','Without AMB')
