%close all; 
clear;

% EMD toolbox can be downloaded online
addpath('/Users/hautiengwu/Dropbox/code/EMD') ;


%========================================
	% some global parameters
initstate(1)
scrsz = get(0,'ScreenSize');
Example = 1 ;
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
Hz = 512 ;
T = 10;
t = [1/Hz:1/Hz:T]' ;
N = length(t) ;
NoiseID = 1 ;
multiplicative = 1 ;

Allprofiles = [...
	% Example 1: close IF and carrier
6 200   3 pi*3/5 0 0     pi/2 3 1 1	0 0 ; ...
	% Example 2: noise
6 10 	3 pi*3/5 0 0 	 pi/2 3 1 1	0 0 ; ...
	% dynamics
6 20000 3 2      0 0	 pi   6 0 0 0.3 0.4 ; ...
] ;


tau = Allprofiles(Example,1) ;
snrdb = Allprofiles(Example,2) ;
am1a = Allprofiles(Example,3) ;
am2a = Allprofiles(Example,4) ;
am1d = Allprofiles(Example,5) ;
am2d = Allprofiles(Example,6) ;
if1a = Allprofiles(Example,7) ;
if2a = Allprofiles(Example,8) ;
if1d = Allprofiles(Example,9) ;
if2d = Allprofiles(Example,10) ;
Ratio1 = Allprofiles(Example,11) ;
Ratio2 = Allprofiles(Example,12) ; 

loadExamples

	% save the ground truth
am0 = zeros(2, N) ; if0 = zeros(2, N) ; y0 = zeros(2, N) ; phi0 = zeros(2, N) ;
am0(1,:) = am1; am0(2,:) = am2 ;
if0(1,:) = if1; if0(2,:) = if2 ;
phi0(1,:) = phi1; phi0(2,:) = phi2 ;
y0(1,:) = y1; y0(2,:) = y2 ;



%======================================
	% start to run BK analysis

Carrier = 20 ; 
yC = y.*exp(sqrt(-1)*2*pi*Carrier*t) ;

allmode = eemd(real(y),0,1);

	% flip to remove the boundary effect
y = [transpose(y) fliplr(y')];
yC = [transpose(yC) fliplr(yC')];

	% apply P_+ to get the holomorphic signal
%f = hilbert(y) ;
mask = ones(1, length(y)) ; 
mask(length(y)/2+1:end) = 0 ;
f = ifft(fft(y).*mask) ;
fC = ifft(fft(yC).*mask) ;


[High, Low, B, G, B_phase, B_phase_der, B_prod] = BKdecomp(f, pass, length(f), over, D, eps, IFmethod) ;
[HighC, LowC, BC, GC, B_phaseC, B_phase_derC, B_prodC] = BKdecomp(fC, pass, length(f), over, D, eps, IFmethod) ;

f = f(1:end/2) ;
High = High(:, 1:end/2) ;
Low = Low(:, 1:end/2) ;
B = B(:, 1:end/2) ;
G = G(:, 1:end/2) ;
B_phase = B_phase(:, 1:end/2) ;
B_phase_der = B_phase_der(:, 1:end/2) ;
B_prod = B_prod(:, 1:end/2) ;

fC = fC(1:end/2) ;
HighC = HighC(:, 1:end/2) ;
LowC = LowC(:, 1:end/2) ;
BC = BC(:, 1:end/2) .*repmat(exp(sqrt(-1)*2*pi*(-Carrier)*t'), pass,1) ;
GC = GC(:, 1:end/2) ;
B_phaseC = B_phaseC(:, 1:end/2) .*repmat(exp(sqrt(-1)*2*pi*(-Carrier)*t'), pass,1) ;
B_phase_derC = B_phase_derC(:, 1:end/2) - Carrier ;
B_prodC = B_prodC(:, 1:end/2) .*repmat(exp(sqrt(-1)*2*pi*(-Carrier)*t'), pass,1) ;



	[tfr0, tfrtic, tfrsq0, ConceFT0, tfrsqtic0] = ConceFT_STFT(transpose(f), 0, 0.5/20, Alpha, 1, Hz*6+1, 1, 6, 1, 0, 0, 0) ;
	[tfr1, tfrtic, tfrsq1, ConceFT1, tfrsqtic] = ConceFT_STFT(transpose(Low(2,:).*B_prod(1,:)), 0, 0.5/20, Alpha, 1, Hz*2.5+1, 1, 6, 1, 0, 0, 0) ;
	[tfr2, tfrtic, tfrsq2, ConceFT2, tfrsqtic] = ConceFT_STFT(transpose(Low(3,:).*B_prod(2,:)), 0, 0.5/20, Alpha, 1, Hz*2.5+1, 1, 6, 1, 0, 0, 0) ;

    [tfr1C, tfrtic, tfrsq1C, ConceFT1C, tfrsqtic] = ConceFT_STFT(transpose(LowC(2,:).*B_prodC(1,:)), 0, 0.5/20, Alpha, 1, Hz*2.5+1, 1, 6, 1, 0, 0, 0) ;
    [tfr2C, tfrtic, tfrsq2C, ConceFT2C, tfrsqtic] = ConceFT_STFT(transpose(LowC(3,:).*B_prodC(2,:)), 0, 0.5/20, Alpha, 1, Hz*2.5+1, 1, 6, 1, 0, 0, 0) ;

		% get IF from SST
		% use the ground truth to speed up the calculation 
		% if not using the ground truth, will be slow but the results
		% are the same
	IFSST = zeros(2, length(f)) ;
	IFSSTC = zeros(2, length(f)) ;
	Q = zeros(31, length(f)) ; cc = zeros(length(f),1) ;
	QC = zeros(31, length(f)) ; 
	for cidx = 1:length(f) ; 
		cc(cidx) = round(if0(1,cidx) ./ (Hz*Alpha)) ;
		Q(:, cidx) = tfrsq1(cc(cidx)-15:cc(cidx)+15, cidx) ;
		QC(:, cidx) = tfrsq1C(cc(cidx)-15:cc(cidx)+15, cidx) ;
	end
	[c] = CurveExt_M(abs(Q'), 1) ; IFSST(1,:) = tfrsqtic(cc'-16+c')*Hz ; 
	[c] = CurveExt_M(abs(QC'), 1) ; IFSSTC(1,:) = tfrsqtic(cc'-16+c')*Hz ; 

    for cidx = 1:length(f) ; 
        cc(cidx) = round(if0(2,cidx) ./ (Hz*Alpha)) ;
        Q(:, cidx) = tfrsq2(cc(cidx)-15:cc(cidx)+15, cidx) ;
        QC(:, cidx) = tfrsq2C(cc(cidx)-15:cc(cidx)+15, cidx) ;
    end
	[c] = CurveExt_M(abs(Q'), 1) ; IFSST(2,:) = tfrsqtic(cc'-16+c')*Hz ; 
	[c] = CurveExt_M(abs(QC'), 1) ; IFSSTC(2,:) = tfrsqtic(cc'-16+c')*Hz ; 

    fprintf(['SER of 1st, without C = ',num2str(norm(real(B_prod(1,:) .* Low(2, :))-real(y0(1,:)))/norm(real(y0(1,:)))),'\n']) ;
    fprintf(['SER of 2nd, without C = ',num2str(norm(real(B_prod(2,:) .* Low(3, :))-real(y0(2,:)))/norm(real(y0(2,:)))),'\n']) ;
    fprintf(['SER of 1st, with C = ',num2str(norm(real(B_prodC(1,:) .* LowC(2, :))-real(y0(1,:)))/norm(real(y0(1,:)))),'\n']) ;
    fprintf(['SER of 2nd, with C = ',num2str(norm(real(B_prodC(2,:) .* LowC(3, :))-real(y0(2,:)))/norm(real(y0(2,:)))),'\n']) ;



    fprintf(['SER of 1st IF, without C = ',num2str(norm(IFSST(1,:)-if0(1,:))/norm(if0(1,:))),'\n']) ;
    fprintf(['SER of 2nd IF, without C = ',num2str(norm(IFSST(2,:)-if0(2,:))/norm(if0(2,:))),'\n']) ;
    fprintf(['SER of 1st IF, with C = ',num2str(norm(IFSSTC(1,:)-if0(1,:))/norm(if0(1,:))),'\n']) ;
    fprintf(['SER of 2nd IF, with C = ',num2str(norm(IFSSTC(2,:)-if0(2,:))/norm(if0(2,:))),'\n']) ;



	if PlotSST 	% apply SST to get IF
			% plot Figure 1
        figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
        subplot(2,3,1) ;
		plot(t, real(f), 'k', 'linewidth', 2) ; axis tight ; set(gca,'fontsize',20) ; xlabel('Time (sec)') ;
		
		subplot(2,3,2) ; 
        plot(t, real(B_prod(1,:) .* Low(2, :)), 'color', [1 0 0], 'linewidth', 2) ; hold on ; xlabel('Time (sec)') ;
        plot(t, real(y0(1,:)), 'k', 'linewidth', 2) ; set(gca,'fontsize', 20) ;

        subplot(2,3,3) ;
        plot(t, real(B_prod(2,:) .* Low(3, :)), 'color', [1 0 0], 'linewidth', 2) ; hold on ;xlabel('Time (sec)') ;
        plot(t, real(y0(2,:)), 'k', 'linewidth', 2) ; set(gca,'fontsize', 20) ;

		subplot(2,3,4) ;
		plot(t, real(B_prod(2,:) .* Low(3, :))-real(y0(2,:))-3, 'color', [1 0 0], 'linewidth', 2) ; hold on
		plot(t, real(B_prodC(2,:) .* LowC(3, :))-real(y0(2,:))-3, 'color', 'k', 'linewidth', 2) ; hold on
		plot(t, real(B_prod(1,:) .* Low(2, :))-real(y0(1,:)), 'color', [1 0 0], 'linewidth', 2) ; hold on
		plot(t, real(B_prodC(1,:) .* LowC(2, :))-real(y0(1,:)), 'color', 'k', 'linewidth', 2) ; hold on
		axis tight ; set(gca,'fontsize', 20) ;xlabel('Time (sec)') ;
		set(gca,'ytick',[-4:1]);
		set(gca,'yticklabel',[-1 0 1 -1 0 1]);

		subplot(2,3,5) ; 
        plot(t, real(B_prodC(1,:) .* LowC(2, :)), 'color', [1 0 0], 'linewidth', 2) ; hold on ; xlabel('Time (sec)') ;
        plot(t, real(y0(1,:)), 'k', 'linewidth', 2) ; set(gca,'fontsize', 20) ;

        subplot(2,3,6) ;
        plot(t, real(B_prodC(2,:) .* LowC(3, :)), 'color', [1 0 0], 'linewidth', 2) ; hold on ; xlabel('Time (sec)') ;
        plot(t, real(y0(2,:)), 'k', 'linewidth', 2) ; set(gca,'fontsize', 20) ;


			% plot Figure 2
        figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
        subplot(2,2,2) ;
		plot(t, if0(1,:) , 'k', 'linewidth', 2) ; 
		hold on ; plot(t, if0(2,:) , 'k', 'linewidth', 2) ; 
		plot(t, IFSST(1,:) , 'r', 'linewidth', 2) ;
		plot(t, IFSST(2,:) , 'r', 'linewidth', 2) ;
        colormap(1-gray) ; axis([0.5 9.5 0 1.1*max(if0(end,:))]) ;
		set(gca,'fontsize', 20) ; xlabel('Time (sec)') ; ylabel('Frequency (Hz)') ;

        subplot(2,2,4) ;
        plot(t, if0(1,:) , 'k', 'linewidth', 2) ; 
        hold on ; plot(t, if0(2,:) , 'k', 'linewidth', 2) ;
        plot(t, IFSSTC(1,:) , 'r', 'linewidth', 2) ;
        plot(t, IFSSTC(2,:) , 'r', 'linewidth', 2) ;
        colormap(1-gray) ; axis([0.5 9.5 0 1.1*max(if0(end,:))]) ; 
		set(gca,'fontsize', 20) ; xlabel('Time (sec)') ; ylabel('Frequency (Hz)') ;

        subplot(2,2,1) ;
        imageSQ(t, tfrsqtic*Hz, abs(tfrsq1) + abs(tfrsq2), 0.999) ;
        colormap(1-gray) ; xlabel('Time (sec)') ; ylabel('Frequency (Hz)') ;
		axis([0.5 9.5 0 1.1*max(if0(end,:))]) ; 

        subplot(2,2,3) ;
        imageSQ(t, tfrsqtic*Hz, abs(tfrsq1C) + abs(tfrsq2C), 0.999) ;
		xlabel('Time (sec)') ; ylabel('Frequency (Hz)') ;
        colormap(1-gray) ; axis([0.5 9.5 0 1.1*max(if0(end,:))]) ; 

	end

        % pass iterations will give pass-1 components
	if PlotRecon	% plot decomposed phase and am
    	for pp = 1:2

            % plot each component
            % the i-th component is L_{i+1}B_1B_2...B_i
        	figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)/3])
        	subplot(1,3,1) ;
        	plot(t, real(B_prod(pp,:) .* Low(pp+1, :)), 'color', [1 0 0], 'linewidth', 2) ; hold on
        	plot(t, real(y0(pp,:)), 'k', 'linewidth', 2) ; set(gca,'fontsize', 20) ;

            % plot AM
			subplot(1,3,2) ;
        	plot(t, real(Low(pp+1,:)), 'color', [1 0 0],  'linewidth', 2) ; hold on ;
        	plot(t, am0(pp,:), 'k', 'linewidth', 2) ;
        	axis([-inf inf min([real(Low(pp+1,:)) am0(pp,:)])-0.3 max([real(Low(pp+1,:)) am0(pp,:)])+0.3])
        	set(gca,'fontsize', 20) ;

            % plot IF
        	subplot(1,3,3) ;
			estif = sum(B_phase_der(1:pp,:)*Hz/2/pi, 1) ;
        	plot(t, estif, 'color', [0 0 1], 'linewidth', 2) ; hold on ;
        	plot(t, IFSST(pp,:), 'color', [1 0 0], 'linewidth', 2) ; hold on ;
        	plot(t, if0(pp,:), 'k', 'linewidth', 2) ;
        	axis([0 inf 0 1.1*max(if0(end,:))]) ; set(gca,'fontsize', 20) ;


            % plot phi
        	%subplot(4,1,4) ;
        	%thetab = phihat ; %unwrap(angle(B_prod(pp, :)))/2/pi ;
        	%plot(t, thetab, 'r', 'linewidth', 2) ; hold on ;
        	%plot(t, phi0(pp,:)-phi0(pp,1), 'k', 'linewidth', 2) ;
        	%axis tight ; set(gca, 'fontsize', 20) ;
        	%legend(['estimated \phi_',num2str(pp)], ['\phi_',num2str(pp)]) ;
        	%xlabel('Time (sec)') ;
		end

            % plot B and G
        if PlotBG
            figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)/2])
            gx = real(G(pp,:)); gy=imag(G(pp,:));
            subplot(1, 2, 1),
            plot(gx, gy, 'b'), hold on ; scatter(gx, gy, 50,t,'fill') ;
            axis([-inf inf -inf max(gy)*1.2])  ;
            grid on ; set(gca, 'fontsize', 20) ; colorbar
            legend(['G_',num2str(pp),' in the complex domain']) ;

            bx = real(B(pp,:)); by=imag(B(pp,:));
            subplot(1, 2, 2); plot(bx,by,'b'), hold on ;
            scatter(bx,by,50,t,'fill') ;
            axis([-1.1 1.1 -1.1 1.3]) ; colorbar ; %colormap(1-gray) ;
            set(gca, 'fontsize', 20) ;
            legend(['B_',num2str(pp),' in the complex domain']) ;
        end

    end








%=======




%{
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)/2]) ;

plot(t, if1,'b','linewidth',2)
hold on;
plot(t, if2,'r','linewidth',2)
set(gca,'fontsize',24);
xlabel('Time (sec)');
text(5, if1(end/2)-0.2,'\phi''_1(t)','color','b','fontsize', 24);
text(5, if2(end/2)+0.2,'\phi''_2(t)','color','r','fontsize', 24);
export_fig(['signalIF_Over',num2str(over),'D',num2str(D)],'-transparent');

hold off
plot(t, real(y1)+8, 'color', [.5 .5 .5]) ;
hold on;
plot(t, real(y2)+16, 'k') ;
plot(t, real(f), 'b','linewidth',2) ; set(gca, 'fontsize', 24) ; hold;
set(gca,'fontsize',24);
xlabel('Time (sec)'); axis tight
export_fig(['signal_Over',num2str(over),'D',num2str(D)],'-transparent');

hold off
am2(end)=am2(end-1);
plot(t, am1, 'b','linewidth',2)
hold on;
plot(t, am2,'r','linewidth',2)
text(5, am1(end/2)-0.2,'A_1(t)','color','b','fontsize', 24);
text(5, am2(end/2)+0.2,'A_2(t)','color','r','fontsize', 24);
set(gca,'fontsize',24); xlabel('Time (sec)');
export_fig(['signalAM_Over',num2str(over),'D',num2str(D)],'-transparent');


subplot(2,1,1) ;
plot(t, real(B_prod(1,:) .* Low(2, :)), 'r', 'linewidth', 2) ; hold on
plot(t, real(y0(1,:)), 'k', 'linewidth', 2) ; set(gca,'fontsize', 24) ;

subplot(2,1,2) ;
plot(t, real(B_prod(2,:) .* Low(3, :)), 'r', 'linewidth', 2) ; hold on
plot(t, real(y0(2,:)), 'k', 'linewidth', 2) ; set(gca,'fontsize', 24) ;
export_fig(['signalDecomp_Over',num2str(over),'D',num2str(D)],'-transparent');


subplot(1,3,1) ; imageSQ(t, tfrsqtic*Hz, abs(tfrsq0), 0.999) ;
hold on; colormap(1-gray) ; axis([1 9 2 4]) ; ylabel('SST') ;

subplot(1,3,2) ; imageSQ(t, tfrsqtic*Hz, abs(tfrsq1) + abs(tfrsq2), 0.999) ;
hold on; colormap(1-gray) ; 
axis([1 9 2 4]) ; ylabel('BG-SST') ;

subplot(1,3,3) ; imageSQ(t, tfrsqtic*Hz, abs(tfrsq1) + abs(tfrsq2), 0.999) ;
hold on; colormap(1-gray) ; plot(t, if0(1,:) , 'r', 'linewidth', 2) ;
plot(t, if0(2,:) , 'r--') ; axis([1 9 2 4]) ; ylabel('BG-SST') ;
export_fig(['signalSST_Over',num2str(over),'D',num2str(D)],'-transparent');
%}

figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)/2])
subplot(3,3,1) ; plot(t, allmode(:,2), 'k','linewidth',2) ; axis tight ; set(gca,'fontsize', 20) ;
subplot(3,3,4) ; plot(t, allmode(:,3), 'k','linewidth',2) ; axis tight ; set(gca,'fontsize', 20) ;
subplot(3,3,7) ; plot(t, allmode(:,4), 'k','linewidth',2) ; axis tight ; set(gca,'fontsize', 20) ;
subplot(3,3,[2 5 8]) ; imageSQ(t, tfrsqtic0*Hz, abs(tfrsq0), 0.999) ;
colormap(1-gray) ; xlabel('Time (sec)') ; ylabel('Frequency (Hz)') ;
axis([0.5 9.5 0 1.1*max(if0(end,:))]) ;
subplot(3,3,[3 6 9]) ; imageSQ(t, tfrsqtic0*Hz, abs(tfrsq0), 0.999) ;
hold on; plot(t, if0(1,:) , 'r', 'linewidth', 2) ;
plot(t, if0(2,:) , 'r', 'linewidth', 2) ;
colormap(1-gray) ; xlabel('Time (sec)') ; ylabel('Frequency (Hz)') ;
axis([0.5 9.5 0 1.1*max(if0(end,:))]) ;



