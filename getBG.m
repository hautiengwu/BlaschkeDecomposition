function [B_ana, G_ana] = getBG(f, N, over, eps);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	f is a 1-D analytic function (only pos. freq.) 
%	Decompose in F = F0 + B*G
%	Filter low frequencies only 1 coeff.
%%%
%	B and G are filtered in the fourier space using N
% 	phase_b is the phase of b, using the oversample b!
%%%
%	over: oversampling rate
%%%	
%	eps: value for threshold
%%%
	[l,ll] = size(f);
	ll2    = floor(ll / 2) ;
	unit   = [1:ll];
	unit_m = [1:over * ll];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Plus oversampling by over

	mask_N = zeros(l,ll*over);
	mask_N = (unit_m < (N+1));	
		
	f_fft = fft(f);
	f_fft_m = zeros(l, ll * over);
	f_fft_m(1:ll2) = f_fft(1:ll2) * over;
	f_ana_m = ifft(f_fft_m) ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: evaluate G

		% compute log(abs(Hi)) of the analytic function Hi
	f_abs_m    = abs(f_ana_m);
	if min(f_abs_m) < eps * max(f_abs_m)
		fprintf('close to zero! add eps\n') ;
		eps2 = (eps * max(f_abs_m))^2;
		log_abs_f  = 0.5 * log(f_abs_m.^2 + eps2);
	else
		eps2 = 0 ;
		log_abs_f  = log(f_abs_m) ;
	end


		% Ana_log_abs_Hi = P_+(ln|Hi|)
	m       = mean(log_abs_f);
	Ana_log_abs_F = m + 2 * ifft(fft(log_abs_f - m) .* mask_N);

%filter the exp_plus
	G_ana_m = ifft(fft(exp(Ana_log_abs_F)) .* mask_N);

	norm_errg = 100 * (1 - norm(fft(exp(Ana_log_abs_F)) .* mask_N) / norm(fft(exp(Ana_log_abs_F))));
	


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: evaluate B
	B_ana_m = f_ana_m ./ G_ana_m;
	
		% (Optional) convoles B, suggested by Nahon
	if 0
		B_ana_m = ifft( fft(B_ana_m) .* mask_N .* (1 - 4/N).^(unit_m/over));
	end

	




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 5: from signal oversample creates signal with original size

	B_ana  = zeros(l,ll);
	G_ana  = zeros(l,ll);

	B_ana(1:ll)  = B_ana_m(1:over:ll*over);
	G_ana(1:ll)  = G_ana_m(1:over:ll*over);
	
	
%%%%%%%%%%%%%%%%%%%%---END---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
