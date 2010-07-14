function values = dis_sol(par)
%     beta=par.beta/(1+par.mum);
%     r=(par.r).*(1+par.mum);
%     beta=par.beta;
%     r=par.r;
% %     values = (1./(1-beta)^2) .*( log(r) + beta.*log(beta)) + 1./(1-beta).*log(1-beta)  + 1./(1-beta).*log(a)    ;
%     values = (1./(1-beta)^2) .*( log((1+r)*(1-beta*r)) + beta.* log( beta / (1-beta) ) ) + 1./(1-beta).*log(a) ;
    b=1./(1-par.beta.*(1-par.mum));
    
    a= ( log((1+par.r)./(1-par.mum).*(1-par.beta.*(1-par.mum))) + par.beta.*(1-par.mum).*b .* log((1+par.r).*par.beta ) ) .*b;
    
    values = a + b.*log(par.a);
    

end
    