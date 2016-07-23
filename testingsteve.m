x = -10:.0001:10;
x1 = normpdf(x,0,1);
x2 = normpdf(x,0,.5);
x3 = normpdf(x,0,.25)

figure;
hold on;
plot(x,x1,'r');
plot(x,x2,'b');
plot(x,x3,'g');
legend('\sigma = 1','\sigma = .5','\sigma = .25');
hold off;

figure;
hold on;
plot(x,x1-x2);
plot(x,x1-x3);
legend('1-.5','1-.25');
hold off;

sum(abs(x1-x2))
