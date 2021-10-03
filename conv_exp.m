function [y, y_over] = conv_exp(f, over, eps) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	[y, y_over] = conv_exp(f, over, eps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Filter fft(f) that corresponds to convole f 
%	with an exponential, using the ovrsampling over
% 	and the parameter 'eps' for exp(-eps.|x|)
%
%	by default: eps  = 10;	
%				over = 1;
%
%	Function used in BG.m
%
%	Created 20/10/99 Michel Nahon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	[N,l]  = size(f) ;
	ll     = round(N/ over) ;
	unit_m = [1:N];
	over;
	if (nargin < 2)
		eps  = 10;
		over = 1;
	elseif (nargin < 3)
		eps = 10;
	end
	
% 	computes FFT of B
	f_fft      = fft(f');
	unit_shift = mod(unit_m + (N / 2) - 1, N) - N / 2;
	convoled   = (1 - eps / ll) .^ (abs(unit_shift) / over);
	f_conv     = ifft(f_fft .* convoled)';
	y_over     = f_conv(1:N);
	y          = f_conv(1:over:N);
%%%%%%%%%%%%%%%%%%%%%---END---%%%%%%%%%%%%%%%%%%%%%
