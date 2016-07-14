function l = diffMapLift(newVal, evec, oldData)
%sort eigenvalues and eigenvectors
[sevec,ind] = sort(evec,'descend');
sdata = oldData(:,ind);

[val,idx] = min(abs(sevec - newVal));
if(val == newVal)
    l = sdata(:,idx);
else
    if(val < newVal)
        lowidx = idx;
    else
        lowidx = idx-1;
    end
    
    lowdist = newVal - sevec(lowidx);
    highdist = sevec(lowidx + 1) - newVal;
    
    l = (lowdist/(lowdist + highdist))*sdata(:,lowidx)+...
        (highdist/(lowdist + highdist))*sdata(:,lowidx+1);
end
end