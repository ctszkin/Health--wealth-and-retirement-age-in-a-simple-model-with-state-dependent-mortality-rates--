function s = stateGen(par,n)
    lengthOfHeathly = geornd(par.muh+par.theta,[n 1]) + 1;
    lengthOfMorbid = ( geornd(par.mum,[n 1]) + 1 ) .* binornd(1,par.theta./(par.muh+par.theta),n,1);

    s = [lengthOfHeathly lengthOfMorbid];
end
