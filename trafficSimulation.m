function trafficSimulation( )
    v0 = .91;
    h = 1.2;
    len = 60;
    numCars = 60;
    mu = .1;
    finalTime = 5 * 10^4;
    
    %% initialize car positions and velocities
    cars = zeros(2*numCars, 1);
    for i = 1:numCars
        cars(i) = (i-1) * len/numCars + mu*sin(2*pi*i/numCars);
        cars(i+numCars) = optimalVelocity(len/numCars);
    end
    
    allTime = rungeKutta(@microsystem, reshape(cars, 2*numCars,1));
%     fCarsEnd = circshift(allTime(1:numCars, end),-1);
%     hEnd = mod(fCarsEnd - allTime(1:numCars,end),len);
%     fCarsStart = circshift(allTime(1:numCars, 1),-1);
%     hStart = mod(fCarsStart - allTime(1:numCars,1),len)
    hEnd = getHeadways(allTime(1:numCars,end));
    hStart = getHeadways(allTime(1:numCars,1));
    figure;
    hold on;
    plot(1:1:numCars,hEnd);
    plot(1:1:numCars,hStart, '--or');
    hold off;
    
    
    
    
    function v = optimalVelocity(headway)
        v = v0 * (tanh(headway - h) + tanh(h));
    end

    function hways = getHeadways(v)
        futureCars = circshift(v,-1);
        hways = mod(futureCars - v, len);
    end

    %% ODE that governs individual cars
    function u = microsystem(colCars)
        invT = 1.7;
        %colCars = reshape(cars, 2*numCars,1);
%         futureCars = zeros(numCars,1);
        colCars(1:numCars,1) = mod(colCars(1:numCars,1),len);
        futureCars = circshift(colCars(1:numCars,1),-1,1);
        headways = mod(futureCars - colCars(1:numCars,1),len);
        
        u = zeros(2*numCars,1);
        u(1:numCars,1) = colCars(numCars+1:2*numCars,1);
        u(numCars+1:2*numCars,1) = invT*(optimalVelocity(headways) - colCars(numCars+1:2*numCars,1));
    end

    %% Runge-Kutta
    function xr = rungeKutta(f, u)
        N = 500000;
        xr = zeros(2*numCars,N+1);
        xr(:,1) = u;
        x = u;
        ts = .1;
        
        for j=1:N
            k1 = f(x);
            k2 = f(x+ts/2*k1);
            k3 = f(x+ts/2*k2);
            k4 = f(x+ts*k3);
            x = x+ts/6*(k1+2*k2+2*k3+k4);
            xr(:,j+1) = x;
        end
    end
end