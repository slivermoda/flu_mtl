clear; clc;
% x = (-180:0.1:180);
% y = normpdf(x, 0, 20);
% x1 = (-180:0.1:0);
% y1 = normpdf(x1, -180, 50);
% x2 = (0:0.1:180);
% y2 = normpdf(x2, 180, 50);
% figure, plot(x, y, '-.');
% hold on
% plot(x1, y1, '-.');
% hold on
% plot(x2, y2, '-.');
% grid on
%%
x = (-8:0.01:8);
y = normpdf(x, 3, 1);
y1 = normpdf(x, -3, 1);
figure, plot(x, y, '-.');
hold on
plot(x, y1, '-.');
grid on