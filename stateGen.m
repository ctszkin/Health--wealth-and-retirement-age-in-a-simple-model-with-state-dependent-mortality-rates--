function s = stateGen(par,n)
    % parSet(rho,lambdam,lambdah,r,w,muh,mum,theta,sigma,a);
    % par = parSet(0.04,0.3,0.2,0.05,0.15,0.001,0.02,0.08,1,(0:0.1:20)');

    % method 1: using a geometric random variable for the simulation
    
    % Reason for +2 :  +1 for the first period; +1 for the first time is sucess (turn to morbid or death)   
    lengthOfHeathly = geornd(par.muh+par.theta,[n 1]) + 1;
    lengthOfMorbid = ( geornd(par.mum,[n 1]) + 1 ) .* binornd(1,par.theta./(par.muh+par.theta),n,1);

    s = [lengthOfHeathly lengthOfMorbid];
end
