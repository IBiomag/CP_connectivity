fs    = 100;
nsamp = 6000;
n1    = randn(1,nsamp);
[b,a] = butter(4,2*[3 6]/fs);
sig0  = filtfilt(b,a,n1);
x     = sig0 + 1*randn(1,nsamp);
y     = sig0 + 1*randn(1,nsamp);

% spectral_new estimation
opt      = [];
opt.maxf = 10;
opt.tap  = 0;
freq{1}  = cp_spectral_new(x, y, fs, 1, opt);
opt.tap  = 1;
freq{2}  = cp_spectral_new(x, y, fs, 1, opt);
opt.tap  = 2;
freq{3}  = cp_spectral_new(x, y, fs, 1, opt);
freq{4}  = cp_spectral_new(x, y, fs, 2);
freq{5}  = cp_spectral_new(x, y, fs, 3);
freq{6}  = cp_spectral_new(x, y, fs, 4);

opt      = [];
opt.maxf = 10;
opt.tap  = 0;
freq2{1}  = CP_spectral(x, y, fs, 1, opt);
opt.tap  = 1;
freq2{2}  = CP_spectral(x, y, fs, 1, opt);
opt.tap  = 2;
freq2{3}  = CP_spectral(x, y, fs, 1, opt);
freq2{4}  = CP_spectral(x, y, fs, 2);
freq2{5}  = CP_spectral(x, y, fs, 3);
freq2{6}  = CP_spectral(x, y, fs, 4);
