function l = diffMapLift(newVal, evec, eval, eps, oldData)
[~,closeidx] = min(abs(evec-newVal));
close = oldData(:,closeidx);

soptions = saoptimset('ObjectiveLimit',1e-7);

tic;
[l,fval] = simulannealbnd(@(x)toMin(x,newVal,evec,eval,eps,oldData),close,...
    zeros(length(close),1),60*ones(length(close),1),...
    soptions);
fval
toc;


    function f = toMin(x,target,evec,eval,eps,old)
        lambda = 1;
        f = lambda * (diffMapRestrict(x,eval,evec,old,eps)-target)^2;

        f = f + abs(sum(x) - 60);
    end
end