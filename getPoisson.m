function [Poisson] = getPoisson(F, C)

    Poisson = [] ;
   
    Theta = exp(2*pi*sqrt(-1)*[0:length(F)-1]/length(F)) ;
    tau = 2*pi/length(F) ;

    for qq = 1: length(C)

        r = abs(C(qq)) ; 

			% Poisson integral \frac{1}{2pi} * tau * \sum f.*\frac{1-r^2}{1-2*r*cc+r^2}
		if r

			z = C(qq) ./ r ;
			cc = real(z*conj(Theta)) ;
        	Pll = (1-r.^2) ./ (1 - 2*r*cc + r.^2) ; 
        	Poisson(qq) = tau * F * Pll' ./ (2*pi) ;

		else

			Poisson(qq) = mean(F)/2/pi ;

		end

    end
