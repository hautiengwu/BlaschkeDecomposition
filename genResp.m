close all; 
clear;
addpath('/Users/hautiengwu/Dropbox/code/EMD') ;


%========================================
	% some global parameters
initstate(1)
scrsz = get(0,'ScreenSize');
Example = 2 ;
PlotSST = 1 ;
PlotRecon = 0 ;
PlotBG = 0 ;



%========================================
	% BK analysis parameters
    % nb passes for reconstruction
pass = 3 ;

    % Ratio for oversampling (over=1 no oversampling)
over = 16 ;

    % Degree of polynomial for approximation (D=1 => cste)
D = 5 ;

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
load flow ;
t = [1:length(flow)]/25 ;

y = hilbert(flow'-mean(flow)) ;
y = y .* exp(sqrt(-1)*2*pi*5*t) ;


Hz = 25 ;
N = length(t) ;


	% flip to remove the boundary effect
y = [y fliplr(conj(y))];

	% apply P_+ to get the holomorphic signal
%f = hilbert(y) ;
mask = ones(1, length(y)) ; 
mask(length(y)/2+1:end) = 0 ;
f = ifft(fft(y).*mask) ;


[High, Low, B, G, B_phase, B_phase_der, B_prod] = BKdecomp(f, pass, length(f), over, D, eps, IFmethod) ;

f = f(1:end/2) .* exp(-sqrt(-1)*2*pi*5*t);
High = High(:, 1:end/2) ;
Low = Low(:, 1:end/2) ;
B = B(:, 1:end/2)  .* repmat(exp(-sqrt(-1)*2*pi*5*t),3,1) ;
G = G(:, 1:end/2) ;
B_phase = B_phase(:, 1:end/2)  .* repmat(exp(-sqrt(-1)*2*pi*5*t),3,1);
B_phase_der = B_phase_der(:, 1:end/2) ;
B_prod = B_prod(:, 1:end/2) .* repmat(exp(-sqrt(-1)*2*pi*5*t),3,1) ;




[tfr0, tfrtic, tfrsq0, ConceFT0, tfrsqtic0] = ConceFT_STFT(transpose(f), 0, 0.03, Alpha, 1, Hz*5*8+1, 1, 6, 1, 0, 0, 0) ;
[tfr1, tfrtic, tfrsq1, ConceFT1, tfrsqtic] = ConceFT_STFT(transpose(Low(2,:).*B_prod(1,:)), 0, 0.03, Alpha, 1, Hz*5*8+1, 1, 6, 1, 0, 0, 0) ;
[tfr2, tfrtic, tfrsq2, ConceFT2, tfrsqtic] = ConceFT_STFT(transpose(Low(3,:).*B_prod(2,:)), 0, 0.03, Alpha, 1, Hz*5*8+1, 1, 6, 1, 0, 0, 0) ;


			% plot Figure 1
	figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
	subplot(3,1,1) ;
	plot(t, flow, 'k', 'linewidth', 2) ; axis tight ; set(gca,'fontsize',20) ; ylabel('(a.u.)') ;
	
	subplot(3,1,2) ; plot(t, real(B_prod(1,:) .* Low(2, :)), 'color', 'k', 'linewidth', 2) ; axis tight ; set(gca,'fontsize',20) ; ylabel('(a.u.)') ;

	subplot(3,1,3) ;  plot(t, real(B_prod(2,:) .* Low(3, :))-30, 'color', 'k', 'linewidth', 2) ; axis tight ; 
	set(gca,'fontsize',20) ; ylabel('(a.u.)') ; 

	figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)*2/3])
	subplot(2,1,1) ; imageSQ(t, tfrsqtic*Hz, abs(tfrsq0), 0.999) ;
	hold on; colormap(1-gray) ; axis([-inf inf 0 0.8]) 
	ylabel('Freq (Hz)') ;

	subplot(2,1,2) ; imageSQ(t, tfrsqtic*Hz, sqrt(abs(tfrsq1).^2 + abs(tfrsq2).^2), 0.999) ; hold on; colormap(1-gray) ; xlabel('Time (sec)') ; axis([-inf inf 0 0.8]) ;
	ylabel('Freq (Hz)') ;







%=======


