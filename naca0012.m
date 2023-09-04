clear
clc

L=0.4;
x_1=0:0.001:0.099;
x_2=0.1:0.01:L;
x=[x_1,x_2];
t_max=0.125*L;
y_2= + t_max/0.2.*(0.298222773*sqrt(x/L) - 0.127125232*x/L - 0.357907906*(x/L).^2 + 0.291984971*(x/L).^3 - 0.105174606*(x/L).^4);
y_1= - t_max/0.2.*(0.298222773*sqrt(x/L) - 0.127125232*x/L - 0.357907906*(x/L).^2 + 0.291984971*(x/L).^3 - 0.105174606*(x/L).^4);

plot(x,y_1,x,y_2)
xlim([-0.01,0.11])
M_X=[x';flipud(x')];
M_Y=[y_2';flipud(y_1')];
M=[M_X,M_Y];
axis equal
hold on
xlswrite('D:\1.fluent\[2]self-propelled fish\[1]fish_s_start\naca0012\naca0012.xls',M)
