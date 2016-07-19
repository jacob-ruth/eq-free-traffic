function l = diffMapLift(newVal, evec, eval, eps, oldData)

fprintf('lifting to %d \n',newVal);

[~,closeidx] = min(abs(evec-newVal));
close = oldData(:,closeidx);

soptions = saoptimset('ObjectiveLimit',1e-7);

if(newVal < min(evec))
    newVal
end


[l,fval] = simulannealbnd(@(x)toMin(x,newVal,evec,eval,eps,oldData),close,...
    zeros(length(close),1),60*ones(length(close),1),...
    soptions);
% fprintf('f value = %d \n',fval);


    function f = toMin(x,target,evec,eval,eps,old)
        lambda = 1;
        f = lambda * (diffMapRestrict(x,eval,evec,old,eps)-target)^2;

        f = f + abs(sum(x) - 60);
    end
end