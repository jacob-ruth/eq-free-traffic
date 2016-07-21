function l = diffMapLift(newVal, evec, eval, eps, oldData)

fprintf('lifting to %d \n',newVal);

closeidxs = dsearchn(evec,newVal');
closeidx = closeidxs(1);
close = oldData(:,closeidx);

fprintf('initial guess: %d \n', evec(closeidx,:));

soptions = saoptimset('ObjectiveLimit',1e-14);


[l,~] = simulannealbnd(@(x)toMin(x,newVal,evec,eval,eps,oldData),close,...
    zeros(length(close),1),60*ones(length(close),1),...
    soptions);

norm(l - close)

    function f = toMin(x,target,evec,eval,eps,old)
        lambda = 100;
        f = lambda * norm(diffMapRestrict(x,eval,evec,old,eps)-target)^2;
        f = f + abs(sum(x) - 60);
    end
end