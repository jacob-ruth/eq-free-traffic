function trafficSimulation()
    h = 1.2;
    len = 60;
    numCars = 60;
    mu = .1;
    finishedEvolution = 200000;
    tskip = 300;
    delta = 2000;
    stepSize = .001;
    delSigma = 0.00001;
    delv0 = 0.00001;
    tolerance = 10^(-12);
    %% initialize car positions and velocities
    cars_1 = zeros(2*numCars, 1);
    cars_2 = zeros(2*numCars, 1);

    origv01 = 0.91;
    origv02 = 0.9;
        options = odeset('AbsTol',10^-8,'RelTol',10^-8);
   % c = approximateWaveSpeed(origv02);
    v0_base1 = origv01;
    v0_base2 = origv02;
    
    for i = 1:numCars
        cars_1(i) = (i-1) * len/numCars + mu*sin(2*pi*i/numCars);
        cars_1(i+numCars) = optimalVelocity(len/numCars, v0_base1);
        
        cars_2(i) = (i-1) * len/numCars + mu*sin(2*pi*i/numCars);
        cars_2(i+numCars) = optimalVelocity(len/numCars, v0_base2);
    end
    
%      tic;
%      [t1,allTime_1] = ode45(@microsystem,[0 finalTime],cars_1, options, v0_base1);
%      ref_1 = allTime_1(end,:)';
%      [t2,allTime_2] = ode45(@microsystem,[0 finalTime],cars_2, options, v0_base2);
%      ref_2 = allTime_2(end,:)';
%      toc; 
%     save('refStats.mat','ref_1','ref_2');
      load('refStats919.mat','ref_1','ref_2');
      hways1 = getHeadways(ref_1(1:numCars));
      hways2 = getHeadways(ref_2(1:numCars));
      [f,~] = microFJ([hways1; -0.8581 ; 0], [hways2 ; -0.8581; 0], numCars, len, v0_base2, 1/1.7)
      
     % std(getHeadways(ref_1(1:numCars)))
     % std(getHeadways(ref_2(1:numCars)))
    %load('things.mat', 'allOfTheThings');
    %std(getHeadways(allOfTheThings(1:numCars, 100)))
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
    steps = 250;
    bif = zeros(2,steps);
    
    %allOfTheThings = zeros(2*numCars, 500);
    %thingCounter = 1;
    %thingLabel = zeros(1,500);
    
    sigma_1 = std(getHeadways(ref_1(1:numCars)));
    sigma_2 = std(getHeadways(ref_2(1:numCars)));
    
    for iEq=1:steps
        iEq
        w = [sigma_2 - sigma_1 ; v0_base2 - v0_base1];
        newGuess = [sigma_2; v0_base2] + stepSize *(w/norm(w));

        %% Newton and that other guy's method
        u = newGuess;
        firstGuess = u;
        f = F(ref_2, u(1),u(2));
        neww = w(1)*(u(1)-newGuess(1)) + w(2)*(u(2) - newGuess(2));
        k=1;
        
        while((abs(f)>tolerance || abs(neww)>tolerance) && k < 20)
%             fprintf('starting iteration %f \n', k)
            f = F(ref_2, u(1),u(2));
%             fprintf('f is %d \n', f);
            Df = jacobian(ref_2, u(1), u(2), w);
            neww = w(1)*(u(1)-newGuess(1)) + w(2)*(u(2) - newGuess(2));

            u = u - Df^(-1)*[f;neww];
            k = k + 1;
        end        
        u
        
%         u = fsolve(@(u)FW(u,ref_2,w,newGuess), newGuess);

        bif(:,iEq) = u;

        ref_1 = ref_2;
        sigma_1 = sigma_2;
        v0_base1 = v0_base2;
        v0_base2 = u(2);
        [sigma_2,ref_2] = ler(u(1),ref_1,tskip+delta,1,u(2));
    end
    

    %save('normDeltathings.mat','allOfTheThings','thingLabel');

    save('paramsts300t10-6.mat','bif','tolerance','h','len','numCars','mu','finalTime','tskip','delta','stepSize','delSigma','delv0','origv01','origv02');

    
    figure;
    scatter(bif(2,:),bif(1,:),'*');
    
    function fw = FW(u,ref,W,newGuess)
        fw = zeros(2,1);
        fw(1) = F(ref,u(1),u(2));
        fw(2) = W(1)*(u(1)-newGuess(1)) + W(2)*(u(2) - newGuess(2));
    end
    

    function c = approximateWaveSpeed(v0)
            waveCars = zeros(2*numCars, 1);
            for waveI = 1:numCars
                waveCars(waveI) = (waveI - 1) * len/numCars + mu*sin(2*pi*waveI/numCars);
                waveCars(waveI+numCars) = optimalVelocity(len/numCars, v0);
            end
            [t,carsEvolved] = ode45(@microsystem,[0 50000],waveCars, options, v0);
            waveCarsEvolved  = carsEvolved(:,1:numCars)';
                        save('waveEvolved91.mat', 't', 'carsEvolved')

            waveCarsHeadways = getHeadways(waveCarsEvolved);
            [headWay, nextCar] = max(getHeadways(waveCarsEvolved(1:numCars)));
            function [diff, isTerminal, direction] = trigger(~,y,~)
                isTerminal = 1;
                direction = 0;
                diff = prevCarMaxHeadway(y,headWay, mod((nextCar - 1),60));
            end
            stoppingOptions = odeset('AbsTol',10^-8,'RelTol',10^-8, 'Events', @trigger);
            [~, carsEnd,te,ye,ie] = ode45(@microsystem,[0,10000],waveCarsEvolved, stoppingOptions, v0);
            c = -1/te;
    end
    function outDiff = prevCarMaxHeadway(y, headway, carIndex)
        headways = getHeadways(y(1:numCars));
        outDiff = headways(carIndex) - headway;
    end

    %% lift and evolve to relatively steady state
%     init = getHeadways(ref_1(1:numCars)); %init is headways
%     new_s = 100000;
%     olds = 0;
%     while(abs(new_s - olds) > .01)
%         l = lift(.15, 1, init, v0_base1);
%         olds = std(init);
%         [t,evolved] = ode45(@microsystem,[0 tskip+delta],l,options,v0_base1);
%         evolved = evolved(end,:)';
%         init = getHeadways(evolved(1:numCars));
%         new_s = std(init);
%     end
%     
%     figure;
%     hold on
%     plot(1:1:numCars,getHeadways(ref_1(1:numCars)),'or')
%     plot(1:1:numCars,init,'xb')
%     plot(1:1:numCars,getHeadways(l(1:numCars)),'*g')
%     hold off
    
    %% lift, evolve, restrict
    % sigma - the current value of the std, used to seed the lifting
    % ref   - a reference state, used to seed the lifting
    % t     - the duration to evaluate the lifted parameters
    % p     -  parameter which will degrade the accuracy of the lifting
    %               operator.  Should be kept at p = 1 by default.
    % v0    - the optimal velocity parameter for this state
    % RETURNS: 
    % sigma     - the std. of the headways after restricting has occured
    % new_state - the final state of the evolution, which can be used as a
    %               future reference state
    function [sigma,new_state, sigma2, new_state2] = ler(sigma,ref,t,p,v0, tReference)
        lifted = lift(sigma, p, getHeadways(ref(1:numCars)),v0);
         [~,evo] = ode45(@microsystem,[0 t],lifted, options,v0);
         if (nargin > 5)
             [~,evo2] = ode45(@microsystem,[0 tReference],evo(end,1:2*numCars)',options,v0);
             sigma2 = std(getHeadways(evo2(end,1:numCars)'));
         end
         evoCars = evo(end, 1:numCars)';
         sigma = std(getHeadways(evoCars));
         if(nargout >= 2)
             new_state = evo(end,1:2*numCars)';
         end
         if(nargout == 4)
             new_state2 = evo2(end,1:2*numCars)';
         end
    end
    %% Finite Difference Quotient
    %  ref - reference state to base the lifting 
    %  sigma - the current value of the std. at which point to approximate
    %       the time derivative
    %  v0 - the velocity parameter for this state
    %  RETURNS:
    %  dif - the difference which approximates the time derivative
    function dif = F(ref, sigma,v0)
        [r0, ~, r1] = ler(sigma, ref, tskip, 1, v0, delta);
        dif = (r1-r0)/delta;
    end
    %% Jacobian for newton's method
    % ref - The previous reference state used to compute F.
    % sigma - the current value of sigma 
    % v0 - the velocity parameter for this state
    % w - the secant direction
    % J- The Jacobian, which will be given by
    % | F_sigma    F_v0  |
    % | w_sigma    w_vo  |
    function J = jacobian(ref, sigma, v0,w)
        J = zeros(2);
        unchanged = F(ref, sigma, v0);
        J(1,1) = (F(ref, sigma + delSigma, v0) - unchanged)/delSigma;
        J(1,2) = (F(ref, sigma,v0 + delv0) - unchanged)/delv0;
        J(2,:) = w';
    end
    
    %% lifting operator
    % s 	-   the target std. for this state
    % p     -   parameter which will degrade the accuracy of the lifting
    %           operator.  Should be kept at p = 1 by default.
    % hways -   headways of the cars at a reference state
    % v0    -   the goal velocity at this state
    % RETURNS: 
    % new   -   the lifted vector of cars and their positions. The cars
    %           will be reset to new positions with the first car centered 
    %            at 0, and with headways that will give the desired std. of
    %            s.  The cars' new velocities will be given directly by the
    %            optimal velocity function.
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
    %    h         -  a parameter to represent the optimal velocity of the
    %               	car 
    %    headway   - the distance between this car and the car ahead of it
    %    v0        - ideal goal speed of each driver
    %    returns:
    %    v         - the optimal velocity of this car, who will speed or slow to try
    %                   to meet it
    function v = optimalVelocity(headway,v0)
        v = v0 * (tanh(headway - h) + tanh(h));
    end

    %% function to calculate headways of car vector
    %  v   - column vector of the cars' positions
    %  returns:
    %  hways - column vector of cars' headways
    function hways = getHeadways(v)
        futureCars = circshift(v,[-1,0]);
        hways = mod(futureCars - v, len);
    end

    %% ODE that governs individual cars
    % Governs the movement of the individual cars (microvariables)
    % ~         - dummy parameter for time, to allow use in ode45
    % params    - a column vector of the distribution of the cars and their
    %               velocities, of size 2*numCars.  The position of car i
    %               and its velocity are given at num params(i), 
    %               params(i + numcars)
    %             
    function u = microsystem(~,colCars, v0)
        invT = 1.7;
        headways = getHeadways(colCars(1:numCars));
        
        u = zeros(2*numCars,1);
        u(1:numCars,1) = colCars(numCars+1:2*numCars,1);
        u(numCars+1:2*numCars,1) = invT*(optimalVelocity(headways,v0) - colCars(numCars+1:2*numCars,1));
    end

    %% Runge-Kutta
    function xr = rungeKutta(f, u)
        ts = .1;
        N = finishedEvolution/ts;
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