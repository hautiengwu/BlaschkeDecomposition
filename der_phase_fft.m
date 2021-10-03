function y = der_phase_fft(b, over) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	y = der_phaseb_fft(b, over)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Computes first derivative of phase of B using
%	the FFT of B to obtain gradient(B) and then
%	computes gradient(B) / B
%
%
%
%	function used in BG.m
%
%	Created 23/9/99 Michel Nahon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	[l,N] = size(b);
	ll = N / over;
	unit_m = [1:N];

	if (nargin < 2)
		over = 1;
	end
	
% 	computes FFT of B
	b_fft      = fft(b) .* ((unit_m <= ll) | (unit_m > N -ll));
	unit_shift = mod(unit_m + (N / 2) - 1, N) - N / 2;
	convol_exp = 1.;
%    convol_exp = (1 - 10 / ll) .^ (abs(unit_shift) / over);
	der_b      = ifft(b_fft .* unit_shift .* convol_exp) * 2 * pi / N;
	
	y = real(der_b ./ b);
	y = real(conv_exp(y',1,10))';
%%%%%%%%%%%%%%%%%%%%%---END---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
