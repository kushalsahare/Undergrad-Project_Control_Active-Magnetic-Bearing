%Solving the diffrential eqaution using ode23
%options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5]);
y0=zeros(36,1);
[t,y] = ode23(@Matrix_def9,[0 5],y0);
[t1,y1] = ode23(@Matrix_def9_withoutamb,[0 5],y0);

%=========================I AMB ==================
% plot(t,y(:,6),'-',t1,y1(:,6),'-')
% xlabel('time (Sec)')
% ylabel('displacement Y (in m)')
% title('Plot of Displacement vs time for 2 amb and 2 disc (I AMB) with External impulse force','FontSize',12)
% legend('With AMB','Without AMB')
%=========================== II AMB===============
plot(t,y(:,12),'-',t1,y1(:,12),'-')
xlabel('time (Sec)')
ylabel('displacement Y (in m)')
title('Plot of Displacement vs time for 2 amb and 2 disc (II AMB)','FontSize',12)
legend('With AMB','Without AMB')