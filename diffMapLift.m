function l = diffMapLift(newVal, evec, eval, eps, oldData)

fprintf('lifting to %d \n',newVal);

closeidxs = dsearchn(evec,newVal');
closeidx = closeidxs(1);
close = oldData(:,closeidx);

fprintf('initial guess: %d \n', evec(closeidx,:));

soptions = saoptimset('ObjectiveLimit',1e-5,'AnnealingFcn',@annealingboltz,'InitialTemperature',1e-3);
soptions.ReannealInterval = 50;

[l,fval] = simulannealbnd(@(x)toMin(x,newVal,evec,eval,eps,oldData),close,...
    zeros(length(close),1),60*ones(length(close),1),...
    soptions);

    function f = toMin(x,target,evec,eval,eps,old)
        lambda = 10;
        f = lambda * norm(diffMapRestrict(x,eval,evec,old,eps)-target)^2;
        f = f + abs(sum(x) - 60);
        LD = fourdif(length(x),1)*2*pi/length(x);
        phase = (LD * x)'*(x - close);
        f = f + abs(phase);
    end
end