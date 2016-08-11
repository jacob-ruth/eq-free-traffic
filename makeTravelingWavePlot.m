%% This file plots the start and two endpoint of a traveling wave after evolving the traffic system
% as well as the evolution of sigma

options = odeset('AbsTol',10^-8,'RelTol',10^-8); % ODE 45 options
h = 1.2;
v0 = 0.89;
numCars = 60;
len = 60;
mu = 0.1;
cars = zeros(2*numCars, 1);
for i = 1:(numCars)
    cars(i) = (i-1) * len/numCars + mu*sin(2*pi*i/numCars);
end
cars(numCars+1:end) = optimalVelocity(h, getHeadways(cars(1:numCars), len), v0);

[t,allTime] = ode45(@microsystem,[0 50000],cars, options, [v0 len h]);

 figure;
 hold on;
 plot(1:1:numCars, getHeadways(allTime(1,1:numCars)', len), 'k.', 'markersize', 12);
 plot(1:1:numCars, getHeadways(allTime(end-130,1:numCars)', len), 'r.', 'markersize', 12);
 plot(1:1:numCars, getHeadways(allTime(end-30,1:numCars)', len), 'b.', 'markersize', 12);
 xlabel('Car', 'fontsize', 20);
 ylabel('Headway', 'fontsize', 20);
 legend('t = 0', 't = 49,900', 't = 50,000');
 hold off;
 
 sig = std(getHeadways(allTime(:, 1:numCars)', len));
 figure;
 plot(t, sigma);
 xlabel('t', 'fontsize', 20);
 ylabel('\sigma', 'fontsize', 20);
 
 


