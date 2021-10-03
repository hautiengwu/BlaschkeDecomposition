function [High_ana, Low_ana] = getHL(f, N, over, D);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	f is a 1-D analytic function (only pos. freq.) 
%	Filter low frequencies only 1 coeff.
%%%
%	over: oversampling rate
%%%	
%	D = degre for approximation polynomial
%%%
%
	[l,ll] = size(f);
	ll2    = floor(ll / 2) ;
	unit   = [1:ll];
	unit_m = [1:over * ll];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: mask for hilbert transform with degre D, for low frequency part. 
% we suppose that f is already analytic (support freq pos)

		% design the filter
	mask_N = zeros(l,ll*over);
	maskLow = zeros(l,ll*over);
	maskHigh  = zeros(l,ll*over);

	mask_N = (unit_m < (N+1));	
	maskLow  = (unit_m < (D+1));
    maskHigh   = (unit_m > D) .* mask_N;
		
		% Filter in the fourier space with the 2 masks 
		% Plus oversampling by over
	f_fft = fft(f);
	f_fft_m = zeros(l, ll * over);
	f_fft_m(1:ll2) = f_fft(1:ll2) * over;

		% Low_ana is low freq
	Low_ana_m = ifft(f_fft_m .* maskLow);
		% High_ana is high freq
	High_ana_m  = ifft(f_fft_m .* maskHigh); 



	High_ana  = zeros(l,ll);
	Low_ana = zeros(l,ll);

	High_ana(1:ll)  = High_ana_m(1:over:ll*over);
	Low_ana(1:ll) = Low_ana_m(1:over:ll*over);
	
%%%%%%%%%%%%%%%%%%%%---END---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
