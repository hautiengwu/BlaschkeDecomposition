function [High, Low, B, G, B_phase, B_phase_der, B_prod] = BK_decomp(f, pass, N, over, D, eps, IFmethod)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	The main code for Blaschke decomposition. 
%	It is originally prepared by M. Nahon in 2000
% 	and later updated/modified by Hau-tieng Wu in 2016
%%%
%	f: the input analytic signal. Assume it is analytic.
%%%
%	number of passes = pass (get pass-1 components)
%%%
%	N corresponds to the degre of th polynom 
%	in Fourrier, it is used to truncated the 
%	functions B & G.
%%%
%	over: oversampling rate
%%%
%	D: the order of the low pass filter 
%%%
%	eps : value of the threshold to cut f
%%%
%	IFmethod: which method to evaluate the phase
% 		method 1: directly with B phase
% 		method 2: using Grad(B) / B 
% 		method 3: same as method 2 but using Fourier decomposition
%
%

	% prepare for the decomposition
[l,ll] = size(f);

if ll == 1 ; error('The signal must be saved as a row') ; end
if min(l,ll) > 1 ; error('The code only support one channel signal right now') ; end

High      = zeros(pass,ll);
Low       = zeros(pass,ll);
G      	  = zeros(pass,ll);
B      	  = zeros(pass,ll);
B_prod 	  = zeros(pass,ll);
B_phase   = zeros(pass,ll);
B_phase_B = zeros(pass,ll);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first decomposition

[High_ana, Low_ana] = getHL(f, N, over, D) ;
[B_ana, G_ana] = getBG(High_ana, N, over, eps) ;
[phase_B, der_phase_B] = getPhase(B_ana, N, over, IFmethod);

%[B_ana, G_ana, High_ana, Low_ana, phase_B, der_phase_B] = getBG(f, N, over, D, eps, IFmethod);


High(1,:)     	 = High_ana;
Low(1,:)      	 = Low_ana ;
G(1,:)        	 = G_ana;
B(1,:)        	 = B_ana;
B_phase(1,:)  	 = phase_B;
B_phase_der(1,:) = der_phase_B;
B_prod(1,:)   	 = B_ana;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We do 'pass' iterations

for iter = 2: pass

	[High_ana, Low_ana] = getHL(G(iter-1, :), N, over, D) ;
	[B_ana, G_ana] = getBG(High_ana, N, over, eps) ;
	[phase_B, der_phase_B] = getPhase(B_ana, N, over, IFmethod);

	%[B_ana, G_ana, High_ana, Low_ana, phase_B, der_phase_B] = getBG(G(iter-1, :), N, over, D, eps, IFmethod);

	High(iter,:)      	= High_ana;
	Low(iter,:)       	= Low_ana;
	B(iter,:)        	= B_ana;
	B_phase(iter,:)  	= phase_B;
	B_phase_der(iter,:) = der_phase_B;
	B_prod(iter,:)   	= B_prod(iter-1,:) .* B_ana ;
	G(iter,:)        	= G_ana;

end

