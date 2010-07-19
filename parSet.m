function par = parSet(rho,lambdam,lambdah,r,wage,max_edu,muh,mum,theta,sigma,a,edu)
    par.rho=rho;
    par.beta=1/(1+rho);
    par.lambdam=lambdam;
    par.lambdah=lambdah;
    par.r=r;
    par.wage=wage;
    par.max_edu=max_edu;
    par.muh=muh;
    par.mum=mum;
    par.theta=theta;
    par.sigma=sigma;
    par.a=a;
    par.edu=edu;
end