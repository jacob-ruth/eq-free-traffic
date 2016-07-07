load('waveEvolved91.mat');
carPositions = carsEvolved(:,1:60)';
futureCars = circshift(carPositions,[-1,0]);

carHeadways = mod(futureCars - carPositions, 60);
%finds the max headways, and which index they occur
[maxHway,maxHwayIndex] = max(carHeadways);
figure;
plot(t, maxHway');
changePoints = mod([1 maxHwayIndex(1:end - 1)] - maxHwayIndex,60);
changePoints(1) = 0;
L = t(changePoints~=0);
dT = L(2:end) - L(1:end-1);
figure;
scatter(L(1:end-1), dT, '.');

figure;
plot(t, std(carHeadways));
avDt = mean(dT(30000:end));
c = -1/avDt
