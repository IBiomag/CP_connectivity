c=clock;
%parfor k0=1:500,
for k0=1:500
  k0
  perf=zeros(6,6,5,3); % freqs, measure, SNR, perfmeasure
  snr=[0.5 1 1.5 2 2.5];
  for s1=1:5
    
    % simulation of signals and spectral estimation
    fs=100;
    nsamp=6000;
    n1=randn(1,nsamp);
    [b,a]=butter(4,2*[3 6]/fs);
    sig0=filtfilt(b,a,n1);
    x=sig0+snr(s1)*randn(1,nsamp);
    y=sig0+snr(s1)*randn(1,nsamp);
    opt=[];
    opt.maxf=10;
    opt.tap=0;
    freq1=CP_spectral(x,y,fs,1,opt);
    opt.tap=1;
    freq2=CP_spectral(x,y,fs,1,opt);
    opt.tap=2;
    freq3=CP_spectral(x,y,fs,1,opt);
    freq4=CP_spectral(x,y,fs,2,opt);
    freq5=CP_spectral(x,y,fs,3,opt);
    freq6=CP_spectral(x,y,fs,4,opt);
    
    k0
    for k=1:6
      res=CP_measure(freq1,k);k2=1;
      [perf(k2,k,s1,1), perf(k2,k,s1,2)]=plot_shaded(res,0,0);
      perf(k2,k,s1,3)=res.etime;
      res=CP_measure(freq2,k);k2=2;
      [perf(k2,k,s1,1), perf(k2,k,s1,2)]=plot_shaded(res,0,0);
      perf(k2,k,s1,3)=res.etime;
      res=CP_measure(freq3,k);k2=3;
      [perf(k2,k,s1,1), perf(k2,k,s1,2)]=plot_shaded(res,0,0);
      perf(k2,k,s1,3)=res.etime;
      res=CP_measure(freq4,k);k2=4;
      [perf(k2,k,s1,1), perf(k2,k,s1,2)]=plot_shaded(res,0,0);
      perf(k2,k,s1,3)=res.etime;
      res=CP_measure(freq5,k);k2=5;
      [perf(k2,k,s1,1), perf(k2,k,s1,2)]=plot_shaded(res,0,0);
      perf(k2,k,s1,3)=res.etime;
      res=CP_measure(freq6,k);k2=6;
      [perf(k2,k,s1,1), perf(k2,k,s1,2)]=plot_shaded(res,0,0);
      perf(k2,k,s1,3)=res.etime;
    end % 6 connectivity metrics
  end % 5 SNR levels
  r{k0}=perf;
end % iterations
