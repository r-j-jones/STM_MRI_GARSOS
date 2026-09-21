function switch_gpu(index1, index2)
% index1 index2 chooses from 1 to gpuDeviceCount

p = gcp('nocreate'); % check the current numWorker
        
if isempty(p) 
    p = parpool(2);
elseif p.NumWorkers ~= 2
    delete(p);
    p = parpool(2);
end     
        
spmd(2)
    a = labindex;
    if (a == 1) 
        a = index1;
    else 
        a = index2;
    end
    gpuDevice(a);
end
end