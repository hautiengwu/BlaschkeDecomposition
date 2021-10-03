%close all;
clear;



%========================================
% some global parameters
initstate(1)
scrsz = get(0,'ScreenSize');
Example = 1 ;



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
Hz = 256 ;
T = 8;
t = [1/Hz:1/Hz:T]' ;
N = length(t) ;



% save the ground truth
am0 = zeros(2, N) ; if0 = zeros(2, N) ; 
y0 = zeros(2, N) ; phi0 = zeros(2, N) ;

am0(1,:) = ones(1,N); 
am0(2,:) = pi*ones(1,N) ;
if0(1,:) = ones(1,N); 
if0(2,:) = pi*ones(1,N) ;
phi0(1,:) = t ; 
phi0(2,:) = pi*t ;
y0(1,:) = am0(1,:).*exp(i*2*pi*phi0(1,:)) ; 
y0(2,:) = am0(2,:).*exp(i*2*pi*phi0(2,:)) ;
y = y0(1,:) + y0(2,:) ;



%=========

figure 
m = 32 ; 
z = cplxgrid(m) ;
z(m+1, :) = z(m+1, :) * .999;
C = z(:) ;
Poisson = getPoisson(y, C) ;
Poisson = reshape(Poisson, m+1, 2*m+1) ;
%Poisson(end, :) = Poisson(end-1, :) ;
subplot(221) 
cplxmap(z, Poisson) ; shading interp
set(gca,'fontsize',20) ; colorbar ; axis tight ; %title('Poisson integral') ;


rr = 0.999 - [0:1:98]/150 ;
PP = zeros(length(rr), length(t)) ;
Phi = zeros(length(rr), length(t)) ;
Phider = zeros(length(rr), length(t)-1) ;
BB = zeros(length(rr), length(t)) ;
True = zeros(length(rr), length(t)) ;

for ii = 1:length(rr)
        % get Poinsson at radius r
	C = rr(ii) * exp(sqrt(-1)*2*pi*t/T) ;
	P = getPoisson(y, C) ;
	PP(ii, :) = P ;
    
        % get BKD on radius r
	[B_ana, G_ana] = getBG(P, length(t), 16, 1e-4) ;
	%tmp = phase(B_ana) ;
	%tmp2 = tmp(2:end)-tmp(1:end-1) ;
	Phi(ii, :) = phase(B_ana) ;
	BB(ii, :) = B_ana ;
	%Phider(ii, :) = tmp2 ;
end

subplot(222); imagesc(2*[1:1024]/1024, rr, log(1+abs(PP))) ; set(gca,'fontsize', 20); 
xlabel('time (unit: pi)'); ylabel('r') ; colorbar

subplot(223); imagesc(2*[1:1024]/1024, rr, real(BB)) ; set(gca,'fontsize', 20); 
xlabel('time (unit: pi)'); ylabel('r') ; colorbar; 

subplot(224); imagesc(2*[1:1024]/1024, rr, Phi) ; set(gca,'fontsize', 20); 
xlabel('time (unit: pi)'); ylabel('r') ; colorbar; 
colormap(jet)

pause

%======================================
% start to run BK analysis

y = [y fliplr(conj(y))];
mask = ones(1, length(y)) ;
mask(length(y)/2+1:end) = 0 ;
f = ifft(fft(y).*mask) ;


[High, Low, B, G, B_phase, B_phase_der, B_prod] = BKdecomp(f, pass, length(f), over, D, eps, IFmethod) ;

f = f(1:end/2) ;
High = High(:, 1:end/2) ;
Low = Low(:, 1:end/2) ;
B = B(:, 1:end/2) ;
G = G(:, 1:end/2) ;
B_phase = B_phase(:, 1:end/2) ;
B_phase_der = B_phase_der(:, 1:end/2) ;
B_prod = B_prod(:, 1:end/2) ;





%% plot results

% pass iterations will give pass-1 components
for pp = 1:2
    
    % plot each component
    % the i-th component is L_{i+1}B_1B_2...B_i
    figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)/3])
    subplot(2,2,1) ;
    plot(t, real(B_prod(pp,:) .* Low(pp+1, :)), 'color', [1 0 0], 'linewidth', 2) ; hold on
    plot(t, real(y0(pp,:)), 'k', 'linewidth', 2) ; set(gca,'fontsize', 20) ;
    
    % plot AM
    subplot(2,2,2) ;
    plot(t, real(Low(pp+1,:)), 'color', [1 0 0],  'linewidth', 2) ; hold on ;
    plot(t, am0(pp,:), 'k--', 'linewidth', 2) ;
    axis([-inf inf min([real(Low(pp+1,:)) am0(pp,:)])-0.3 max([real(Low(pp+1,:)) am0(pp,:)])+0.3])
    set(gca,'fontsize', 20) ;
    
    % plot IF
    subplot(2,2,3) ;
    estif = sum(B_phase_der(1:pp,:)*Hz/2/pi, 1) ;
    plot(t, estif, 'r', 'linewidth', 2) ; hold on ;
    plot(t, if0(pp,:), 'k--', 'linewidth', 2) ;
    axis([0 inf 0 4*max(if0(end,:))]) ; set(gca,'fontsize', 20) ;

    subplot(2,2,4) ;
    estphase = sum(B_phase(1:pp,:)/2/pi, 1) ;
    plot(t, estphase, 'r', 'linewidth', 2) ; hold on ;
    plot(t, phi0(pp,:), 'k--', 'linewidth', 2) ;
    axis tight ; set(gca,'fontsize', 20) ;
end

% plot B and G
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







