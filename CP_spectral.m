function freq = CP_spectral(x,y,fs,meth,opt)
%function to compute cortico-peripheral connectivity
%opt.maxf highest frequency to compute

if nargin<5
  opt = [];
end
if ~isfield(opt, 'maxf')
  opt.maxf = 10;
end

nsamp=length(x);
switch meth
    case 1
        % multitaper analysis; needs opt.tap (number of tapers); opt.measure {'coh' 'plv' 'wppc'};
        data=[];
        data.label{1}='Env';
        data.label{2}='MEG';
        data.fsample=fs;
        data.dimord='chan_time';
        data.trial{1}=[x; y];
        data.time{1}=[1:nsamp]/fs;
        cfg = [];
        cfg.length = 2;
        cfg.overlap=0.5;
        data2 = ft_redefinetrial(cfg,data);
        ft_warning off
        ft_info off
        
        cfg=[];
        cfg.channel={'all'};
        cfg.channelcmb = {'Env' 'MEG';};
        cfg.method='mtmfft';
        cfg.feedback='no';
        cfg.output='fourier'; %fourier
        cfg.pad='nextpow2';
        if opt.tap == 0,
            cfg.taper='hanning';
        else
            cfg.taper='dpss';
            cfg.tapsmofrq=opt.tap;
        end
        cfg.foi=[1:0.5:opt.maxf];
        cfg.keeptrials='yes';
        freq=ft_freqanalysis(cfg,data2);
        freq.meth=['FT' num2str(opt.tap)];
               
    case 2
        %fb = cwtfilterbank('SignalLength',numel(x),'SamplingFrequency',100,'FrequencyLimits',[1 10],'WaveletParameters',[3,20]); %3 60
        fb = cwtfilterbank('SignalLength',numel(x),'SamplingFrequency',fs,'FrequencyLimits',[1 opt.maxf],'Wavelet','amor'); %3 60

        [Envc,fw]=wt(fb,x);
        [megc,fw]=wt(fb,y);
        freq=[];
        freq.label={'Env'  'MEG'};
        freq.dimord='rpttap_chan_freq';
        freq.freq=flipud(fw);
        freq.fourierspctrm=zeros(size(Envc,2),2,size(Envc,1));
        freq.fourierspctrm(:,1,:)=flip(Envc,1)';
        freq.fourierspctrm(:,2,:)=flip(megc,1)';
        freq.cumsumcnt=ones(size(Envc,2),1);
        freq.cumtapcnt=ones(size(Envc,2),1);
        freq.meth='CWT';
        
    case 3
        % with bandpass filters
%         [~, ftdir] = ft_version;
%         s=which('butter');
%         if ~strcmp(s(1:12),'/Applic.ZIV/')
%             rmpath([ftdir '/external/signal/'])
%         end
        ff=[1:0.5:opt.maxf];
        nfreq=length(ff);

        Envc=zeros(nfreq,length(x));
        megc=Envc;
        
        for k=1:nfreq,
%             if k==1,
%                 [b,a]=butter(4,2*k/fs,'low');
%             else
%                 [b,a]=butter(4,2*[ff(k)-1 ff(k)+1]/fs);
%             end
%             Envc(k,:)=(hilbert(filtfilt(b,a,x)));
%             megc(k,:)=(hilbert(filtfilt(b,a,y)));
            
            if k==1,
            Envc(k,:)=hilbert(ft_preproc_lowpassfilter(x,fs,ff(k)+1,[],'firws'));
            megc(k,:)=hilbert(ft_preproc_lowpassfilter(y,fs,ff(k)+1,[],'firws'));
            else
            Envc(k,:)=hilbert(ft_preproc_bandpassfilter(x,fs,[ff(k)-1 ff(k)+1],[],'firws'));
            megc(k,:)=hilbert(ft_preproc_bandpassfilter(y,fs,[ff(k)-1 ff(k)+1],[],'firws'));
            end
        end
        freq=[];
        freq.label={'Env'  'MEG'};
        freq.dimord='rpttap_chan_freq';
        freq.freq=ff;
        freq.fourierspctrm=zeros(size(Envc,2),2,size(Envc,1));
        freq.fourierspctrm(:,1,:)=Envc';
        freq.fourierspctrm(:,2,:)=megc';
        freq.cumsumcnt=ones(size(Envc,2),1);
        freq.cumtapcnt=ones(size(Envc,2),1);
        freq.meth='BF';

    case 4
        win=hanning(2*fs);
        nov=fs;
        [Envc,fr,t] = spectrogram(x,win,nov,[1:0.5:opt.maxf],fs);
        [megc,fr,t] = spectrogram(y,win,nov,[1:0.5:opt.maxf],fs);
        freq=[];
        freq.label={'Env'  'MEG'};
        freq.dimord='rpttap_chan_freq';
        freq.freq=fr;
        freq.fourierspctrm=zeros(size(Envc,2),2,size(Envc,1));
        freq.fourierspctrm(:,1,:)=Envc';
        freq.fourierspctrm(:,2,:)=megc';
        freq.cumsumcnt=2*fs*ones(size(Envc,2),1);
        freq.cumtapcnt=ones(size(Envc,2),1);
        freq.meth='SG';
end
end

