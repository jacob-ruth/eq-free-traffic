function trafficSimulation()
    h = 1.2;
    len = 60;
    numCars = 60;
    mu = .1;
    finalTime = 200000;
    tskip = 300;
    delta = 2000;
    stepSize = .001;
    delSigma = 0.001;
    delv0 = 0.001;
    
    %% initialize car positions and velocities
    cars_1 = zeros(2*numCars, 1);
    cars_2 = zeros(2*numCars, 1);
    v0_base1 = 0.89;
    v0_base2 = 0.8898;
    
    for i = 1:numCars
        cars_1(i) = (i-1) * len/numCars + mu*sin(2*pi*i/numCars);
        cars_1(i+numCars) = optimalVelocity(len/numCars, v0_base1);
        
        cars_2(i) = (i-1) * len/numCars + mu*sin(2*pi*i/numCars);
        cars_2(i+numCars) = optimalVelocity(len/numCars, v0_base2);
    end
    
    options = odeset('AbsTol',10^-8,'RelTol',10^-8);
% %     tic;
%      [t1,allTime_1] = ode45(@microsystem,[0 finalTime],[cars_1; v0_base1],options);
%      ref_1 = allTime_1(end,:)';
%      [t2,allTime_2] = ode45(@microsystem,[0 finalTime],[cars_2; v0_base2],options);
%      ref_2 = allTime_2(end,:)';
% %     toc;
%     
%      save('refStats.mat','t1','allTime_1','ref_1','t2','allTime_2','ref_2');
     load('refStats.mat','ref_1','ref_2');
    
    %% plot the results
%     hEnd = getHeadways(allTime_2(end,1:numCars)');
%     hStart = getHeadways(allTime_2(1,1:numCars)');
%     figure;
%     hold on;
%     plot(1:1:numCars,hEnd);
%     plot(1:1:numCars,hStart, '--or');
%     hold off;
%     
%     s = std(getHeadways(allTime_2(:,1:numCars)'));
%     
%     figure;
%     plot(t, s);
    
    %% initialize secant continuation
    sigma_1 = std(getHeadways(ref_1(1:numCars)));
    sigma_2 = std(getHeadways(ref_2(1:numCars)));
    w = [sigma_2 - sigma_1 ; v0_base2 - v0_base1];
    newGuess = [sigma_2; v0_base2] + stepSize *(w/norm(w));
    
    %% Newton and that other guy's method
    u = newGuess;
    for k=1:10
        fprintf('starting iteration %f \n', k);
        f = F(ref_2, u(2),u(1));
        fprintf('f is %f \n', f)
        Df = jacobian(ref_2, u(1), u(2), w);
        u = u - Df^(-1)*[f;w(1)*(u(1)-newGuess(1)) + w(2)*(u(2) - newGuess(2))]
    end
    
    %% lift and evolve to relatively steady state
    init = getHeadways(ref_1(1:numCars)); %init is headways
    new_s = 100000;
    olds = 0;
    while(abs(new_s - olds) > .01)
        l = lift(.15, 1, init, v0_base1);
        olds = std(init);
        [t,evolved] = ode45(@microsystem,[0 tskip+delta],[l;v0_base1],options);
        evolved = evolved(end,:)';
        init = getHeadways(evolved(1:numCars));
        new_s = std(init);
    end
    
    figure;
%     hold on
%     plot(1:1:numCars,getHeadways(ref_1(1:numCars)),'or')
%     plot(1:1:numCars,init,'xb')
%     plot(1:1:numCars,getHeadways(l(1:numCars)),'*g')
    hold off
    
    function out = ler(sigma,ref,t,p,v0)
        lifted = lift(sigma, p, getHeadways(ref(1:numCars)),v0);
        [~,evo] = ode45(@microsystem,[0 t],[lifted;v0],options);
        evo = evo(end, 1:numCars)';
        out = std(getHeadways(evo));
    end
    
    function dif = F(ref, sigma,v0)
        r0 = ler(sigma, ref, tskip, 1,v0);
        r1 = ler(sigma, ref, tskip+delta, 1,v0);
        dif = (r1-r0)/delta;
    end
    
    function J = jacobian(ref, sigma, v0,w)
        J = zeros(2);
        unchanged = F(ref, sigma, v0);
        J(1,1) = (F(ref, sigma + delSigma, v0) - unchanged)/delSigma;
        J(1,2) = (F(ref, sigma,v0 + delv0) - unchanged)/delv0;
        J(2,:) = w';
    end

    %% lifting operator
    function new = lift(s, p, hways, v0)
        del = p * s/std(hways) * (hways - mean(hways)) + mean(hways);
        new = zeros(length(hways)*2,1);
        for i=2:length(hways)
            new(i) = mod(sum(del(1:i-1)),len);
        end
        hways = getHeadways(new(1:length(hways)));
        for i=1:length(hways)
            new(i+length(hways)) = optimalVelocity(hways(i), v0);
        end
    end
    
    %% optimal velocity function given in paper
    function v = optimalVelocity(headway,v0)
        v = v0 * (tanh(headway - h) + tanh(h));
    end

    function hways = getHeadways(v)
        futureCars = circshift(v,[-1,0]);
        hways = mod(futureCars - v, len);
    end

    %% ODE that governs individual cars
    function u = microsystem(~,params)
        invT = 1.7;
        colCars = params(1:(end - 1));
        v0 = params(end);
        colCars(1:numCars,1) = mod(colCars(1:numCars,1),len);
        futureCars = circshift(colCars(1:numCars,1),-1,1);
        headways = mod(futureCars - colCars(1:numCars,1),len);
        
        u = zeros(2*numCars,1);
        u(1:numCars,1) = colCars(numCars+1:2*numCars,1);
        u(numCars+1:2*numCars,1) = invT*(optimalVelocity(headways,v0) - colCars(numCars+1:2*numCars,1));
        u(2*numCars+1,1) = 0;
    end

    %% Runge-Kutta
    function xr = rungeKutta(f, u)
        ts = .1;
        N = finalTime/ts;
        xr = zeros(2*numCars,N+1);
        xr(:,1) = u;
        x = u;

        for j=1:N
            k1 = f(0,x);
            k2 = f(0,x+ts/2*k1);
            k3 = f(0,x+ts/2*k2);
            k4 = f(0,x+ts*k3);
            x = x+ts/6*(k1+2*k2+2*k3+k4);
            xr(:,j+1) = x;
        end
    end
end