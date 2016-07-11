function runDiffMap()
theta = linspace(0,4*pi, 100)';
f  = [theta.*cos(theta) theta.*sin(theta)];
figure;
plot(f(:,1), f(:,2),'.');
eps = 0.5;


D = zeros(length(theta));
for r = 1:length(theta)
    for c = 1:length(theta)
        D(r,c) = norm(f(r,:)-f(c,:));
    end
end

D;

stepSize = .0001;
maxEps = 50;

[vec, val] = diffusionMap(eps, D,5);
t = 5;


figure;
plot3(f(:,1), f(:,2),(val(1,1)^t)*vec(:,1));

figure;
plot3(f(:,1), f(:,2),(val(2,2)^t)*vec(:,2));

figure;
plot3(f(:,1), f(:,2),(val(3,3)^t)*vec(:,3));

figure;
plot3(f(:,1), f(:,2),(val(4,4)^t)*vec(:,4));

return;
t = 1;
figure;
plot((val(1,1)^t)*vec(:,1));
title('Time = 1');
figure;
plot((val(2,2)^t)*vec(:,2));
title('Time = 1');


t = 10;
figure;
plot((val(1,1)^t)*vec(:,1), (val(2,2)^t)*vec(:,2));
title('Time = 10');


t = 100;
figure;
plot((val(1,1)^t)*vec(:,1), (val(2,2)^t)*vec(:,2));
title('Time = 100');

%load('trajectory90.mat','allTime');

%allTime = allTime';
%allTime = allTime(:,end-5000:10:end);

eps = .06;
return;

allTime(1:60,:) = getHeadways(allTime(1:60,:));

[minVal1, ~] = min(allTime(1:60,:),[],1);
[maxVal1, max1] = max(allTime(1:60,:),[],1);
[minVal2, ~] = min(allTime(61:end,:),[],1);
[maxVal2, ~] = max(allTime(61:end,:),[],1);
for c = 1:length(max1)
    allTime(1:60,c) = (allTime(1:60,c) - minVal1(c))./(maxVal1(c)-minVal1(c));
    allTime(1:60,c) = circshift(allTime(1:60,c), [-max1(c)+1,0]);
    
    allTime(61:end,c) = (allTime(61:end,c)-minVal2(c))./(maxVal2(c)-minVal2(c));
    allTime(61:end,c) = circshift(allTime(61:end,c), [-max1(c)+1,0]);
end

D = zeros(length(allTime));
for r = 1:length(allTime)
    for c = 1:length(allTime)
        D(r,c) = norm(allTime(:,r)-allTime(:,c));
    end
end

diffusionMap(eps,D);


    function plotEps()
        epsilon = 0:stepSize:maxEps;
        L = zeros(size(epsilon));
        for iEps = 1:length(epsilon)
            curEps = epsilon(iEps);
            L(iEps) = sum(sum(D<curEps,1),2);
        end
        
        figure;
        loglog(epsilon,L);
        xlabel('\epsilon','FontSize',24);
        ylabel('L(\epsilon)','FontSize',20);
    end

    function hways = getHeadways(v)
        futureCars = circshift(v,[-1,0]);
        hways = mod(futureCars - v, 60);
    end
end