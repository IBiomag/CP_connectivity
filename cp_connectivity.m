c=clock;
for k0 = 1:500
  k0
  perf = zeros(6,6,5,3); % freqs, measure, SNR, perfmeasure
  snr  = [0.5 1 1.5 2 2.5];
  for s1 = 1:5
    
    % simulation of signals
    fs    = 100;
    nsamp = 6000;
    n1    = randn(1,nsamp);
    [b,a] = butter(4,2*[3 6]/fs);
    sig0  = filtfilt(b,a,n1);
    x     = sig0 + snr(s1)*randn(1,nsamp);
    y     = sig0 + snr(s1)*randn(1,nsamp);
    
    % spectral estimation
    opt      = [];
    opt.maxf = 10;
    opt.tap  = 0;
    freq1    = CP_spectral(x, y, fs, 1, opt);
    opt.tap  = 1;
    freq2    = CP_spectral(x, y, fs, 1, opt);
    opt.tap  = 2;
    freq3    = CP_spectral(x, y, fs, 1, opt);
    freq4    = CP_spectral(x, y, fs, 2);
    freq5    = CP_spectral(x, y, fs, 3);
    freq6    = CP_spectral(x, y, fs, 4);
    
    freq     = {freq1 freq2 freq3 freq4 freq5 freq6};
    
    % connectivity estimation
    k0
    for k = 1:6
      for m = 1:6
        res = CP_measure(freq{m}, k);
        [perf(m,k,s1,1), perf(m,k,s1,2)] = plot_shaded(res,0,0);
        perf(m,k,s1,3) = res.etime;
      end % 6 spectral estimates
    end % 6 connectivity metrics
  end % 5 SNR levels
  
  r{k0} = perf;
keyboard
end % iterations
