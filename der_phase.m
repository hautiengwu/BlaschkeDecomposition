function y = der_phase(theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   y = der_phase(theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Compute the first derivate of theta, using a 
% 	convolution afterword.
%
% 	y = Dx(theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ex:
%		y = der_phase(f)
%
% Created 15/9/99 Michel Nahon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	filter = [1 2 3 2 1];
	len    = length(filter);
	
%	compute filter
	m      = mean(filter);
	filter = filter / (m * len);
	len2   = (len + 1) / 2;

	[l,ll] = size(theta);

% compute derivative
	der_theta  = gradient(theta);

% convole derivative
	conv_der   = conv(der_theta, filter);
	y          = conv_der(len2:(l + len2 - 1));	
	
%%%%%%%%%%%%%%%%%%%%---END---%%%%%%%%%%%%%%%%%%%%
