%Solving the diffrential eaution using ode23
%options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5]);
y0=zeros(12,1);
[t,y] = ode23(@Matrix_def5,[0 1],y0);
[t1,y1] = ode23(@Matrix_def5_withoutamb,[0 1],y0);
%=========================I AMB ==================
plot(t,y(:,2),'-',t1,y1(:,2),'-')
xlabel('time (Sec)')
ylabel('displacement Y')
title('Plot of Displacement vs time for 2 amb and 2 disc (I AMB)','FontSize',14)
legend('With AMB','Without AMB','FontSize',14)
%=========================== II AMB===============
% plot(t,y(:,4),'-',t1,y1(:,4),'-')
% xlabel('time (Sec)')
% ylabel('displacement Y')
% title('Plot of Displacement vs time for 2 amb and 2 disc (II AMB)','FontSize',10)
% legend('With AMB','Without AMB')