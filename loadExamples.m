
initstate(1)
TN1 = round(N*Ratio1);
TN2 = round(N*Ratio2);

%% the amplitude modulation of the simulated signal
am1 = smooth(cumsum(randn(N,1)) ./ Hz, Hz*7, 'loess') ;
am1 = am1a + am1d * am1 ./ max(abs(am1)) ;
am2 = smooth(cumsum(randn(N,1)) ./ Hz, Hz*6, 'loess') ;
am2 = am2a + am2d * am2 ./ max(abs(am2)) ;
%am3 = am3a ;
if TN1>0 ; am1(1:TN1) = 0 ; end
if TN2>0 ; am2(end-TN2+1:end) = 0 ; end


%% the instantaneous frequency of the simulated signal
		% Example1 = close IFs
Xi1 = cumsum(randn(N,1)) ./ Hz ;
if1 = smooth(Xi1, Hz*5, 'loess') ;
%tmp = smooth(cumsum(randn(N,1)) ./ Hz, Hz*6, 'loess') ;
%if2 = if2a + if1d * if1 ./ max(abs(if1)) + tmp/10 ;
if1 = if1a + if1d * if1 ./ max(abs(if1)) ;

Xi2 = cumsum(randn(N,1)) ./ Hz ;
Xi2(1:end/4-1) = if2a + if2d * Xi2(1:end/4-1) ./ max(abs(Xi2)) ;
Xi2(end*3/4+1:end) = if2a + if2d * Xi2(end*3/4+1:end) ./ max(abs(Xi2)) ;
Xi2(end/4:end*3/4) = if1a + if1d*Xi1(end/4:end*3/4)./max(abs(if1)) + 0.3 ;
if2 = smooth(Xi2, Hz*5, 'loess') ;





%%
phi1 = cumsum(if1) / Hz ;
phi2 = cumsum(if2) / Hz ;
if TN1>0 ; if1(1:TN1) = 0 ; end
if TN2>0 ; if2(end-TN2+1:end) = 0 ; end



	%% the simulated signal.
y1 = am1 .* exp(sqrt(-1)*2*pi*phi1) ;
y2 = am2 .* exp(sqrt(-1)*2*pi*phi2) ;
clean = y1 + y2  ; 

sigma = sqrt( var(clean)*10.^( -snrdb /10 ) );

if Example == 2 ; initstate(iter) ; end
switch NoiseID
    
    case 1
        %% add noise (Gaussian white noise)
        noise = randn(N, 1) ;

    case 2
        %% ARMA
        dis = random('t',4,N, 1) ;
        e = armaxfilter_simulate(dis, .5, 1, .5, 1, -.5) ;
        noise = e ./ std(e) ;

    case 3
        %% poisson
        noise = random('poiss',1,N,1);

end


	%% Add noise
if snrdb < 100 & multiplicative == 0

    y = clean + sigma * noise ;

elseif snrdb < 100 & multiplicative == 1
	
	disp('multiplicative noise') ;
    y = clean .* exp(noise/2) ;

else

	y = clean ;

end
