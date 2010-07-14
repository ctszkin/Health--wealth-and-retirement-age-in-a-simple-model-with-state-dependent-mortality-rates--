function value = cont_sol(par,vs,a)
    if ( nargin<2)
        error('par and vs must be specified');
    end
    if ( nargin<3)
        a=par.a;
    end
    muh=par.muh;
    mum=par.mum;
    theta=par.theta;
    rho=par.rho;
    r=par.r;
    lambda=par.lambda;
    w=par.w;
    
    if strcmp(vs,'hr')
        alpha = 1./(rho+mum) ;
        beta = 1./(rho+mum).* ( log(rho+mum) + (r-rho) / (rho +mum) );

        aa= (1 + theta .* alpha) ./ (rho + muh + theta) ;
        bb= ( aa .* (r + muh) -1 - log(aa) + theta.*beta ) ./ (rho + muh + theta) ;

        value = log(a) .* aa  +  bb;
        return;
    elseif strcmp(vs,'mr')
        value = log(a)./(rho+mum)  +  1./(rho+mum).* ( log(rho+mum) + (r-rho) / (rho +mum) );
        return;
    elseif strcmp(vs,'mw')
        value = log(a + w/(r+mum) )./(rho+mum)  +  1./(rho+mum).* ( log(rho+mum)-lambda + (r-rho) / (rho +mum) );
        return;
    else
        error('vs can only be ''hr'' ''mr'' ''mw''');
    end
end
