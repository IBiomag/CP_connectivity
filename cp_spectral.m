function freq = cp_spectral(x,y,fs,meth,opt)

%function to compute cortico-peripheral connectivity
%opt.maxf highest frequency to compute

if nargin<5
  opt = struct([]);
end
opt.maxf    = ft_getopt(opt, 'maxf', 10);
opt.length  = ft_getopt(opt, 'length', 2);
opt.overlap = ft_getopt(opt, 'overlap', 0.5);

switch meth
    case 1
        % create fieldtrip style data structure
        data          = [];
        data.label    = {'Env';'MEG'};
        data.fsample  = fs;
        data.trial{1} = [x(:)';y(:)'];
        data.time{1}  = (1:length(x))/fs;

        % multitaper analysis; needs opt.tap (tapsmofrq, if 0 use 'hanning')
        opt.taper = ft_getopt(opt, 'taper', 0);
        
        % redefine into shorter overlapping segments
        cfg  = keepfields(opt, {'length' 'overlap'});
        data = ft_redefinetrial(cfg, data);

        ft_warning off
        ft_info off
        
        % spectral estimation 
        cfg            = [];
        cfg.channel    = {'all'};
        cfg.method     = 'mtmfft';
        cfg.feedback   = 'no';
        cfg.output     = 'fourier'; %fourier
        cfg.foi        = (1:0.5:opt.maxf);
        % cfg.pad        = 'nextpow2'; % I don't understand this, because
        % it leads to incomparable frequency bins, e.g. in comparison to
        % bpfilter+hilbert
        if opt.tap == 0
            cfg.taper     = 'hanning';
        else
            cfg.taper     = 'dpss';
            cfg.tapsmofrq = opt.tap;
        end
        freq           = ft_freqanalysis(cfg, data);
        freq.meth      = ['FT' num2str(opt.tap)];
               
    case 2
        fb = cwtfilterbank('SignalLength', numel(x), 'SamplingFrequency', fs, 'FrequencyLimits', [1 opt.maxf], 'Wavelet', 'amor'); %3 60

        [Envc, fw] = wt(fb,x);
        [megc, fw] = wt(fb,y);
        
        % create fieldtrip style freq structure
        freq               = [];
        freq.label         = {'Env'; 'MEG'};
        freq.dimord        = 'rpttap_chan_freq';
        freq.freq          = flipud(fw);
        freq.fourierspctrm = zeros(size(Envc,2),2,size(Envc,1));
        freq.fourierspctrm(:,1,:) = flip(Envc,1).';
        freq.fourierspctrm(:,2,:) = flip(megc,1).';
        freq.cumsumcnt     = ones(size(Envc,2),1);
        freq.cumtapcnt     = ones(size(Envc,2),1);
        freq.meth          = 'CWT';
        
    case 3
        % with bandpass filters
        ff    = (1:0.5:opt.maxf);
        nfreq = length(ff);

        Envc = zeros(nfreq,length(x));
        megc = Envc;
        
        for k = 1:nfreq
            if k==1
                Envc(k,:) = hilbert(ft_preproc_lowpassfilter(x, fs, ff(k)+1, [], 'firws'));
                megc(k,:) = hilbert(ft_preproc_lowpassfilter(y, fs, ff(k)+1, [], 'firws'));
            else
                Envc(k,:) = hilbert(ft_preproc_bandpassfilter(x, fs, [ff(k)-1 ff(k)+1], [], 'firws'));
                megc(k,:) = hilbert(ft_preproc_bandpassfilter(y, fs, [ff(k)-1 ff(k)+1], [], 'firws'));
            end
        end
        
        % create fieldtrip style freq structure
        freq               = [];
        freq.label         = {'Env'; 'MEG'};
        freq.dimord        = 'rpttap_chan_freq';
        freq.freq          = ff;
        freq.fourierspctrm = zeros(size(Envc,2),2,size(Envc,1));
        freq.fourierspctrm(:,1,:) = Envc.';
        freq.fourierspctrm(:,2,:) = megc.';
        freq.cumsumcnt     = ones(size(Envc,2),1);
        freq.cumtapcnt     = ones(size(Envc,2),1);
        freq.meth          = 'BF';

    case 4
        win = hanning(opt.length*fs);
        nov = fs*opt.overlap*opt.length;
        [Envc, fr, t] = spectrogram(x, win, nov, (1:0.5:opt.maxf), fs);
        [megc, fr, t] = spectrogram(y, win, nov, (1:0.5:opt.maxf), fs);
        
        % create fieldtrip style freq structure
        freq               = [];
        freq.label         = {'Env'  'MEG'};
        freq.dimord        = 'rpttap_chan_freq';
        freq.freq          = fr;
        freq.fourierspctrm = zeros(size(Envc,2),2,size(Envc,1));
        freq.fourierspctrm(:,1,:) = Envc';
        freq.fourierspctrm(:,2,:) = megc';
        freq.cumsumcnt     = 2*fs*ones(size(Envc,2),1);
        freq.cumtapcnt     = ones(size(Envc,2),1);
        freq.meth          = 'SG';
end
end

