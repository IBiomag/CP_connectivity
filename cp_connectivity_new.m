c = clock;

nrepeat = 500;
snr     = [0.5 1 1.5 2 2.5];
r       = cell(1,nrepeat); % holds the results per repeat

for k0 = 1:nrepeat
  perf = zeros(6,6,5,3); % freqs, measure, SNR, perfmeasure
  for s1 = 1:5
    freq = cell(1,6);
    
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
    freq{1}  = cp_spectral(x, y, fs, 1, opt);
    opt.tap  = 1;
    freq{2}  = cp_spectral(x, y, fs, 1, opt);
    opt.tap  = 2;
    freq{3}  = cp_spectral(x, y, fs, 1, opt);
    freq{4}  = cp_spectral(x, y, fs, 2);
    freq{5}  = cp_spectral(x, y, fs, 3);
    freq{6}  = cp_spectral(x, y, fs, 4);
    
    % connectivity estimation
    for k = 1:6
      for m = 1:6
        res = cp_measure(freq{m}, k);
        [perf(m,k,s1,1), perf(m,k,s1,2)] = plot_shaded(res,0,0);
        perf(m,k,s1,3) = res.etime;
      end % 6 spectral estimates
    end % 6 connectivity metrics
  end % 5 SNR levels
  
  r{k0} = perf;
end % iterations
