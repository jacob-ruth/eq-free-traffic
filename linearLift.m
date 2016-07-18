function l = diffMapLift(newVal, evec, oldData)
%sort eigenvalues and eigenvectors
[sevec,ind] = sort(evec,'ascend');
sdata = oldData(:,ind);
idx = 1;

while(sevec(idx) < newVal)
    idx = idx + 1;
end
lowidx = idx - 1;
val = sevec(idx);

if(val == newVal)
    l = sdata(:,idx);
else    
    lowdist = newVal - sevec(lowidx);
    highdist = sevec(lowidx + 1) - newVal;
    
    l = (lowdist/(lowdist + highdist))*sdata(:,lowidx)+...
        (highdist/(lowdist + highdist))*sdata(:,lowidx+1);
end
end