
t = 2*pi*[1:1024]/1024;
	% root 1/2, i/3, -sqrt(2)/3+i*sqrt(2)/3
%x = (exp(sqrt(-1)*t)-0.5).*(exp(sqrt(-1)*t)-sqrt(-1)/3) ; %.*(exp(sqrt(-1)*5*t)-(sqrt(2)/3+sqrt(2)*sqrt(-1)/3)^5);
%x = (exp(sqrt(-1)*t)-0.5).^20 ;
%x = (exp(sqrt(-1)*t)-0.5).*(exp(sqrt(-1)*t)-sqrt(-1)*.6).*(exp(sqrt(-1)*t)-(sqrt(2)/3+sqrt(2)*sqrt(-1)/3));

zC = exp(sqrt(-1)*t) ;
x = (zC.^10 - 0.4.^10) .* (zC - 0.7*i) .* (zC + 0.7) .* (zC + 0.2 - 0.1*i) ; 
%x = (zC.^4 - 0.5.^4) .* (zC + 0.1 + 0.2*i) ; 



if 0 

	Hz = 128 ; T = 8; t = [1/Hz:1/Hz:T] ;
	N = length(t) ;
	NoiseID = 1 ; multiplicative = 0 ; iter = 1; 

	tau = 6 ; snrdb = 200 ; Ratio1 = 0 ; Ratio2 = 0 ;
	am1a = 3 ; am2a = pi*3/5; am1d = 0 ; am2d = 0 ;
	if1a = pi/2 ; if2a = 3 ; if1d = 1 ; if2d = 1 ;

	loadExamples ;
	x = transpose(y) ;

end




%C = [] ;
%r = log2(linspace(2^(0.01),2^(0.99),160)) ;
%for qq = 1: length(r)
%    z = exp(sqrt(-1)*2*pi*[0:16*(length(r)-qq+1)-1]/(16*(length(r)-qq+1))) ;
%    C = [C r(qq)*z] ;
%end

m = 32 ; 
z = cplxgrid(m) ;
C = z(:) ;
Poisson = getPoisson(x, C) ;
Poisson = reshape(Poisson, m+1, 2*m+1) ;
Poisson(end, :) = Poisson(end-1, :) ;
subplot(221);
cplxmap(z, Poisson) ; shading interp
set(gca,'fontsize',20) ; colorbar ; axis tight ; %title('Poisson integral') ;


rr = 0.99 - [1:98]/100 ;
PP = zeros(length(rr), length(t)) ;
Phi = zeros(length(rr), length(t)) ;
Phider = zeros(length(rr), length(t)-1) ;
BB = zeros(length(rr), length(t)) ;
True = zeros(length(rr), length(t)) ;

for ii = 1:length(rr)
	C = rr(ii) * exp(sqrt(-1)*t) ;
	P = getPoisson(x, C) ;
	PP(ii, :) = P ;
	[B_ana, G_ana] = getBG(P, length(t), 16, 1e-4) ;
	tmp = phase(B_ana) ;
	tmp2 = tmp(2:end)-tmp(1:end-1) ;
	Phi(ii, :) = phase(B_ana) ;
	BB(ii, :) = B_ana ;
	Phider(ii, :) = tmp2 ;
end

subplot(222); imagesc(2*[1:1024]/1024, rr, log(1+abs(PP))) ; set(gca,'fontsize', 20); xlabel('time (unit: pi)'); ylabel('r') ; colorbar
%title('real part of |real(P_r)|^{0.5}sign(real(P_r))')

subplot(223); imagesc(2*[1:1024]/1024, rr, real(BB)) ; set(gca,'fontsize', 20); xlabel('time (unit: pi)'); ylabel('r') ; colorbar %title('real(B)') ;

subplot(224); imagesc(2*[1:1024]/1024, rr, Phi) ; set(gca,'fontsize', 20); xlabel('time (unit: pi)'); ylabel('r') ; colorbar %title('phase of B') ;
colormap(jet)
