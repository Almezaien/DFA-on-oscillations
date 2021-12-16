load dfa_data.mat

% try with narrowband amplitude time series
xfilt = abs(hilbert(filterFGx(x,srate,10,5)));

% create data with DFA=.5
N = length(x);
randnoise = randn(N,1);


% setup parameters
nScales = 20;
ranges  = round(N*[.01 .2]);
scales  = ceil(logspace(log10(ranges(1)),log10(ranges(2)),nScales));
rmses   = zeros(2,nScales);

% plot the two signals
figure(11), clf
subplot(221)
plot(timevec,randnoise)
title('Signal 1: white noise')
xlabel('Time (seconds)')

subplot(222)
plot(timevec,xfilt)
title('Signal 2: real data')
xlabel('Time (seconds)')



% integrate and mean-center the signals
randnoise = cumsum(randnoise(:)-mean(randnoise));
x4dfa = cumsum(xfilt(:)-mean(xfilt));

% and show those time series for comparison
subplot(223)
plot(timevec,randnoise)
title('Integrated noise')
 
subplot(224)
plot(timevec,x4dfa)
title('Integrated signal')

%%

% compute RMS over different time scales
for scalei = 1:nScales
    
    % number of epochs for this scale
    n = floor(N/scales(scalei)); 
    
    % compute RMS for the random noise
    epochs  = reshape( randnoise(1:n*scales(scalei)) ,scales(scalei),n);
    depochs = detrend(epochs);
    % here is the root mean square computation
    rmses(1,scalei) = mean( sqrt( mean(depochs.^2,1) ) );
    
    % repeat for the signal
    epochs  = reshape( x4dfa(1:n*scales(scalei)) ,scales(scalei),n);
    depochs = detrend(epochs);
    rmses(2,scalei) = mean( sqrt( mean(depochs.^2,2) ) );
end


% fit a linear model to quantify scaling exponent
A = [ ones(nScales,1) log10(scales)' ];  % linear model
dfa1 = (A'*A) \ (A'*log10(rmses(1,:))'); % fit to noise
dfa2 = (A'*A) \ (A'*log10(rmses(2,:))'); % fit to signal


% plot the 'linear' fit (in log-log space)
figure(12), clf, hold on

% plot results for white noise
plot(log10(scales),log10(rmses(1,:)),'rs','linew',2,'markerfacecolor','w','markersize',10)
plot(log10(scales),dfa1(1)+dfa1(2)*log10(scales),'r--','linew',2)

% plot results for the real signal
plot(log10(scales),log10(rmses(2,:)),'bs','linew',2,'markerfacecolor','w','markersize',10)
plot(log10(scales),dfa2(1)+dfa2(2)*log10(scales),'b--','linew',2)

legend({'Data (noise)';[ 'Fit (DFA=' num2str(round(dfa1(2),3)) ')' ]; ...
        'Data (signal)';[ 'Fit (DFA=' num2str(round(dfa2(2),3)) ')' ] })
xlabel('Data scale (log)'), ylabel('RMS (log)')
title('Comparison of Hurst exponent for different noises')
axis square

%% DFA scanning through frequencies

frex = linspace(1,40,80);
dfas = zeros(size(frex));
rms  = zeros(1,nScales);

for fi=1:length(frex)
    
    % get power time series
    x4dfa = abs(hilbert(filterFGx(x,srate,frex(fi),5)));
    x4dfa = cumsum(x4dfa-mean(x4dfa));
    
    % compute RMS over different time scales
    for scalei = 1:nScales
        
        % number of epochs for this scale
        n = floor(N/scales(scalei));
        
        % repeat for signal1
        epochs  = reshape( x4dfa(1:n*scales(scalei)) ,scales(scalei),n);
        depochs = detrend(epochs);
        rms(scalei) = mean( sqrt( mean(depochs.^2,1) ) );
    end
    
    dfa = (A'*A) \ (A'*log10(rms)');
    dfas(fi) = dfa(2);
end

figure(14), clf
plot(frex,dfas,'ko-','markerfacecolor','w')
xlabel('Frequency (Hz)')
ylabel('Hurst exponent')