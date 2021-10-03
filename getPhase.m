function [phase_B, der_phase_B] = getPhase(f, N, over, IFmethod)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	f is a 1-D analytic function (only pos. freq.) 
%%%
%	over: oversampling rate
%%%	
%	D = degre for approximation polynomial
%%%
%	IFmethod: the method chosen to evaluate phase and IF
%
	[l,ll] = size(f);
	ll2    = floor(ll / 2) ;
	unit   = [1:ll];
	unit_m = [1:over * ll];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% oversampling by over
	f_fft = fft(f);
	f_fft_m = zeros(l, ll * over);
	f_fft_m(1:ll2) = f_fft(1:ll2) * over;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: compute derivative of phase of b with different algorithms
	B_ana_m = ifft(f_fft_m) ;

	if IFmethod == 1
			% method 1: directly with B phase
		thetaB     = unwrap(angle(B_ana_m));
		der_thetaB = gradient(thetaB);

	elseif IFmethod == 2
			% method 2: using Grad(B) / B 	
		der_thetaB = imag(gradient(B_ana_m) ./ B_ana_m) .* (abs(B_ana_m) > eps);
	
		thetaB = zeros(over*ll,1);
		for k = 1:(over*ll-1),
			thetaB(k+1) = thetaB(k) + der_thetaB(k+1);
		end


	elseif IFmethod == 3
			% method 3: same way but using Fourier decomposition
		der_thetaB =  der_phase_fft(B_ana_m, over) .* (abs(B_ana_m) > eps);
		thetaB = zeros(over*ll,1);

		coef           = (1 - 4/N).^(abs(unit_m - 1 - (over * ll)/2) / over);
		der_thetab_per = zeros(over * ll, 1);
		thetab_per     = zeros(over * ll, 1);
		
		der_thetab_per            = der_thetaB;
		der_thetab_per(ll * over) = 0;
		der_thetab_per            = ifft(fft(der_thetab_per) .* fftshift(coef));
	
		for k = 1:(ll * over - 1),
			thetab_per(k+1) = thetab_per(k) + der_thetab_per(k);
		end	
	
		conv_phase    = thetab_per';
		b_conv_phase  = exp(i * conv_phase);

		der_thetaB = der_phase_fft(b_conv_phase, over) .* (abs(B_ana_m) > eps);
        thetaB = zeros(over*ll,1);
        for k = 1:(over*ll-1),
            thetaB(k+1) = thetaB(k) + der_thetaB(k+1);
        end

	end
	
		



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 5: from signal oversample creates signal with original size

	phase_B  = zeros(l,ll);
	der_phase_B = zeros(1,ll) ;

	phase_B(1:ll) = thetaB(1:over:ll*over);	
	der_phase_B(1:ll) = over * der_thetaB(1:over:over*ll) ;




%%%%%%%%%%%%%%%%%%%%---END---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
