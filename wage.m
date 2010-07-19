function w = wage(e,par)

    e(e>=par.max_edu) = par.max_edu;

    w = exp(par.wage*e);
    
end