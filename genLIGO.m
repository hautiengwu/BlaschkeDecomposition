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
load HanfordH1.txt ;
t = HanfordH1(:,1) ; t = [t; 2*t(end)-t(end-1)] ; t = t' ;

y = zeros(1, size(HanfordH1,1)+1) ;
y(1:size(HanfordH1,1)) = hilbert(HanfordH1(:,2)'-mean(HanfordH1(:,2))) ;
y(end) = y(end-1) ; y = y .* exp(sqrt(-1)*2*pi*20*t) ;

load LivingstonL1.txt ;
yL = zeros(1, size(LivingstonL1,1)+1) ;
yL(1:size(LivingstonL1,1)) = hilbert(LivingstonL1(:,2)'-mean(LivingstonL1(:,2))) ;
yL(end) = yL(end-1) ; yL = yL  .* exp(sqrt(-1)*2*pi*20*t) ;


Hz = 1./(t(2)-t(1)) ;
N = length(t) ;


	% flip to remove the boundary effect
y = [y fliplr(conj(y))];
yL = [yL fliplr(conj(yL))];

	% apply P_+ to get the holomorphic signal
%f = hilbert(y) ;
mask = ones(1, length(y)) ; 
mask(length(y)/2+1:end) = 0 ;
f = ifft(fft(y).*mask) ;
fL = ifft(fft(yL).*mask) ;


[High, Low, B, G, B_phase, B_phase_der, B_prod] = BKdecomp(f, pass, length(f), over, D, eps, IFmethod) ;
[HighL, LowL, BL, GL, B_phaseL, B_phase_derL, B_prodL] = BKdecomp(fL, pass, length(f), over, D, eps, IFmethod) ;

f = f(1:end/2)  .* exp(-sqrt(-1)*2*pi*20*t);
High = High(:, 1:end/2) ;
Low = Low(:, 1:end/2) ;
B = B(:, 1:end/2)  .* repmat(exp(-sqrt(-1)*2*pi*20*t),3,1) ;
G = G(:, 1:end/2) ;
B_phase = B_phase(:, 1:end/2)  .* repmat(exp(-sqrt(-1)*2*pi*20*t),3,1);
B_phase_der = B_phase_der(:, 1:end/2) ;
B_prod = B_prod(:, 1:end/2) .* repmat(exp(-sqrt(-1)*2*pi*20*t),3,1) ;


fL = fL(1:end/2) .* exp(-sqrt(-1)*2*pi*20*t) ;
HighL = HighL(:, 1:end/2) ;
LowL = LowL(:, 1:end/2) ;
BL = BL(:, 1:end/2) .* repmat(exp(-sqrt(-1)*2*pi*20*t),3,1) ;
GL = GL(:, 1:end/2) ;
B_phaseL = B_phaseL(:, 1:end/2) .* repmat(exp(-sqrt(-1)*2*pi*20*t),3,1) ;
B_phase_derL = B_phase_derL(:, 1:end/2) ;
B_prodL = B_prodL(:, 1:end/2) .* repmat(exp(-sqrt(-1)*2*pi*20*t),3,1) ;



[tfr0, tfrtic, tfrsq0, ConceFT0, tfrsqtic0] = ConceFT_STFT(transpose(f), 0, 0.03, Alpha, 1, 801, 1, 6, 1, 0, 0, 0) ;
[tfr1, tfrtic, tfrsq1, ConceFT1, tfrsqtic] = ConceFT_STFT(transpose(Low(2,:).*B_prod(1,:)), 0, 0.03, Alpha, 1, 801, 1, 6, 1, 0, 0, 0) ;
[tfr2, tfrtic, tfrsq2, ConceFT2, tfrsqtic] = ConceFT_STFT(transpose(Low(3,:).*B_prod(2,:)), 0, 0.03, Alpha, 1, 801, 1, 6, 1, 0, 0, 0) ;

[tfr0L, tfrtic, tfrsq0L, ConceFT0L, tfrsqtic0] = ConceFT_STFT(transpose(fL), 0, 0.03, Alpha, 1, 801, 1, 6, 1, 0, 0, 0) ;
[tfr1L, tfrtic, tfrsq1L, ConceFT1L, tfrsqtic] = ConceFT_STFT(transpose(LowL(2,:).*B_prodL(1,:)), 0, 0.03, Alpha, 1, 801, 1, 6, 1, 0, 0, 0) ;
[tfr2L, tfrtic, tfrsq2L, ConceFT2L, tfrsqtic] = ConceFT_STFT(transpose(LowL(3,:).*B_prodL(2,:)), 0, 0.03, Alpha, 1, 801, 1, 6, 1, 0, 0, 0) ;


			% plot Figure 1
	figure ;
	subplot(2,3,1) ;
	plot(HanfordH1(:,1), HanfordH1(:,2), 'k', 'linewidth', 2) ; axis tight ; set(gca,'fontsize',20) ; xlabel('Time (sec)') ; ylabel('Strain (10^{-21})') ;
	
	subplot(2,3,2) ; plot(t, real(B_prod(1,:) .* Low(2, :)), 'color', 'k', 'linewidth', 2) ; axis tight ; set(gca,'fontsize',20) ; xlabel('Time (sec)') ; ylabel('Strain (10^{-21})') ;

	subplot(2,3,3) ; plot(t, real(B_prod(2,:) .* Low(3, :)), 'color', 'k', 'linewidth', 2) ; axis tight ; set(gca,'fontsize',20) ; xlabel('Time (sec)') ; ylabel('Strain (10^{-21})') ;

    subplot(2,3,4) ;
    plot(LivingstonL1(:,1), LivingstonL1(:,2), 'k', 'linewidth', 2) ; axis tight ; set(gca,'fontsize',20) ; xlabel('Time (sec)') ; ylabel('Strain (10^{-21})') ;

    subplot(2,3,5) ; plot(t, real(B_prodL(1,:) .* LowL(2, :)), 'color', 'k', 'linewidth', 2) ; axis tight ; set(gca,'fontsize',20) ; xlabel('Time (sec)') ; ylabel('Strain (10^{-21})') ;

    subplot(2,3,6) ; plot(t, real(B_prodL(2,:) .* LowL(3, :)), 'color', 'k', 'linewidth', 2) ; axis tight ; set(gca,'fontsize',20) ; xlabel('Time (sec)') ; ylabel('Strain (10^{-21})') ;




	figure 
	subplot(2,3,1) ; imageSQ(t, tfrsqtic*Hz, abs(tfrsq0), 0.999) ;
	hold on; colormap(1-gray) ; xlabel('Time (sec)') ; axis([-inf inf 0 300]) ;
	ylabel('Freq (Hz)') ;

	subplot(2,3,2) ; imageSQ(t, tfrsqtic*Hz, sqrt(abs(tfrsq1).^2 + abs(tfrsq2).^2), 0.999) ; hold on; colormap(1-gray) ; xlabel('Time (sec)') ; axis([-inf inf 0 300]) ;
	ylabel('Freq (Hz)') ;

	subplot(2,3,3) ; imageSQ(t, tfrsqtic*Hz, sqrt(abs(tfrsq1).^2 + abs(tfrsq2).^2), 0.999) ; hold on; colormap(1-gray) ; xlabel('Time (sec)') ; axis([0.35 0.45 20 180]) ;
	ylabel('Freq (Hz)') ;

    subplot(2,3,4) ; imageSQ(t, tfrsqtic*Hz, abs(tfrsq0L), 0.999) ;
    hold on; colormap(1-gray) ; xlabel('Time (sec)') ; axis([-inf inf 0 300]) ;
	ylabel('Freq (Hz)') ;

    subplot(2,3,5) ; imageSQ(t, tfrsqtic*Hz, sqrt(abs(tfrsq1L).^2 + abs(tfrsq2L).^2), 0.999) ; hold on; colormap(1-gray) ; xlabel('Time (sec)') ; axis([-inf inf 0 300]) ;
	ylabel('Freq (Hz)') ;

    subplot(2,3,6) ; imageSQ(t, tfrsqtic*Hz, sqrt(abs(tfrsq1L).^2 + abs(tfrsq2L).^2), 0.999) ; hold on; colormap(1-gray) ; xlabel('Time (sec)') ; axis([0.35 0.45 20 180]) ;
	ylabel('Freq (Hz)') ;







%=======


