function displayTraffic(u, L)
[numVars,N] = size(u);
carPositions = u(1:numVars/2, :);
carPositions = rem(carPositions, L);
stepSize = 4;
specialCars = [1, numVars/8, numVars/4, 3*numVars/8];
close all; clc
figure
for i = 1:stepSize:N
    angle = carPositions(:, i)*(2*pi)/L;
    x = cos(angle);
    y = sin(angle);
    clf;
    hold on;
    plot(x,y,'or');
    plot(x(specialCars), y(specialCars), 'xb');
    shg;
    pause(.01)
end
end