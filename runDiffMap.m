function runDiffMap()
<<<<<<< HEAD
load('tracyOutput.mat','output9');
load('rebeccaOutput.mat','outputRebecca');

allTime = output9;
=======
load('tracyOutput','output9');
load('rebeccaOutput.mat','outputRebecca');

allTime = output9;

allTime = allTime(:,end-5000:10:end);
>>>>>>> origin/DiffusionMaps

allTime = allTime(:,1:5:end);

eps = .7;

stepSize = .0001;
maxEps = 1;

allTime = getHeadways(allTime(1:60,:));


[minVal1, ~] = min(allTime(1:60,:),[],1);
[maxVal1, max1] = max(allTime(1:60,:),[],1);
[minVal2, ~] = min(allTime(61:end,:),[],1);
[maxVal2, ~] = max(allTime(61:end,:),[],1);
for c = 1:length(max1)
    allTime(1:60,c) = (allTime(1:60,c) - minVal1(c))./(maxVal1(c)-minVal1(c));
    allTime(1:60,c) = circshift(allTime(1:60,c), [-max1(c)+1,0]);
    
%     allTime(61:end,c) = (allTime(61:end,c)-minVal2(c))./(maxVal2(c)-minVal2(c));
%     allTime(61:end,c) = circshift(allTime(61:end,c), [-max1(c)+1,0]);
end

% figure;
% for iPlot = 1:5:size(allTime,2)
%     plot(allTime(1:60,iPlot));
%     pause;
% end


D = zeros(length(allTime));
for r = 1:length(allTime)
    for c = 1:length(allTime)
        D(r,c) = norm(allTime(:,r)-allTime(:,c));
    end
end

k=2;

t = 3;
[vec,val] = diffusionMap(eps,D,k);

val = val^t;

figure;
scatter3(1:1:size(vec,2),vec(1,:),vec(2,:),'.');
xlabel('point');
ylabel('vec1');
zlabel('vec2');





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