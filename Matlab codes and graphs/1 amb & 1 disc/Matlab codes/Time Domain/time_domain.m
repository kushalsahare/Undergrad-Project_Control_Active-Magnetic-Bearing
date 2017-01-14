%Solving the diffrential eaution using ode23
%options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5]);
y0=zeros(12,1);
[t,y] = ode23(@Matrix_def,[0 10],y0);
[t1,y1] = ode23(@Matrix_def_withoutamb,[0 10],y0);
plot(t,y(:,2),'-',t1,y1(:,2),'-')
xlabel('time (Sec)')
ylabel('displacement Y (in m)')
title('Plot of Displacement vs time for 1 amb and 1 disc','FontSize',10)
legend('With AMB','Without AMB')



