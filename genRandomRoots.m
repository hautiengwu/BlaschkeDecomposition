close all;  clear;
addpath('/Users/hautiengwu/Dropbox/code/EMD') ;


%========================================
	% some global parameters
scrsz = get(0,'ScreenSize');
PlotRecon = 0 ;



%========================================
	% BK analysis parameters
    % nb passes for reconstruction
pass = 10 ;

    % Ratio for oversampling (over=1 no oversampling)
over = 16 ;

    % Degree of polynomial for approximation (D=1 => cste)
D = 1 ;

    % Value of epsilon as a threshold for f.
eps = 1e-4 ;

    % four methods to evaluate phase and IF 
    % method 1: directly with B phase
    % method 2: using Grad(B) / B 
    % method 3: same way but using Fourier decomposition
IFmethod = 3 ;


%========================================
	% parameters for SST, if used
Alpha = 5e-5 ;

%========================================
	% Signal parameters
Hz = 4096 ;
T = 1;
t = 2*pi*[1/Hz:1/Hz:T] ;
N = length(t) ;

RN = 30 ; sigma = 1 ;
ROOT = (randn(RN, 1) + i* randn(RN, 1))*sigma ;
y = ones(1, N) ;
for jj = 1:RN 
%    y = y.*(exp(i*t) - ROOT(jj)./norm(ROOT(jj))) ; %./(1-exp(i*t)*conj(ROOT(jj))) ;
    y = y.*(exp(i*t) - ROOT(jj)) ;
    if max(abs(y)) > 1e9 ; y = 1e9*y./max(abs(y)) ; end
end
subplot(121) ; plot(real(y)) ; 
subplot(122) ; hold off; plot(cos(t), sin(t), 'color', [.7 .7 .7]) ; hold on ; plot(real(ROOT), imag(ROOT), 'x') ; axis image




[High, Low, B, G, B_phase, B_phase_der, B_prod] = BKdecomp(y, pass, length(y), over, D, eps, IFmethod) ;




% pass iterations will give pass-1 components
figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)/3]) ;
plot(t/pi, real(y), 'k', 'linewidth', 2) ; set(gca,'fontsize', 20) ; hold on;
recon = zeros(size(y)) ;
for pp = 1:pass-1
    recon = recon + B_prod(pp,:) .* Low(pp+1, :) ;
    % plot each component
    % the i-th component is L_{i+1}B_1B_2...B_i
    if pp == 1 | pp == pass-1 ; plot(t/pi, real(recon), 'color', [1 pp/pass pp/pass], 'linewidth', 2) ; hold on ; end
    
end
axis tight ;

figure;
plot(cos(t), sin(t), 'color', [.7 .7 .7], 'linewidth', 2) ; hold on; axis equal
plot(ROOT, 'x') ;
