homeDir = '/analyse/cdhome/';
homeDir = '/home/chrisd/';
addpath([homeDir '/fieldtrip-20170315/'])
addpath(genpath([homeDir '/wavelet/']))
addpath(genpath([homeDir '/gcmi-master/']))

proj0012Dir = '/analyse/Project0012/';

fs = 100;
nSamp = 6000;
nTrials = 500;
allNoiseAmps = [0.5 1 1.5 2 2.5];
allSpecEstMeths = {'1','2','3','4','5','6'};
allConnMeths = {'PLV','gcmi','circ_R','wppc','coh','entropy'};
noiseTypes = {'GWN','1/f'};

gtPassBand = [3 6];
gtFilterOrder = 4;

% set up filter for ground truth effect
[b,a] = butter(gtFilterOrder,2*gtPassBand./fs);

allPerf = zeros(numel(allSpecEstMeths),numel(allConnMeths),numel(allNoiseAmps),3,numel(noiseTypes),nTrials);
allRes = cell(numel(noiseTypes),nTrials);

initparclus(16)

c = clock;

for nn = 1:numel(noiseTypes)
    
    parfor tr = 1:nTrials

        % preallocate perf variable for this trial
        perf = zeros(numel(allSpecEstMeths),numel(allConnMeths),numel(allNoiseAmps),3); % freqs, measure, SNR, perfmeasure
        tmpRes = cell(numel(allConnMeths),numel(allSpecEstMeths),numel(allNoiseAmps));
        
        for na = 1:numel(allNoiseAmps)

            % sample noise
            n1 = randn(1,nSamp);
            % generate signal as filtered noise
            sig0 = filtfilt(b,a,n1);

            % simulate source and target by adding different noise to signal
            switch noiseTypes{nn}
                case 'GWN'
                    sourceNoise = randn(1,nSamp);
                    targetNoise = randn(1,nSamp);
                case '1/f'
                    sourceNoise = arbssnoise(nSamp,-1/2);
                    targetNoise = arbssnoise(nSamp,-1/2);   
            end
            
            x = sig0 + allNoiseAmps(na)*sourceNoise;
            y = sig0 + allNoiseAmps(na)*targetNoise;


            opt = [];
            opt.tap = 0;
            opt.maxf = 10;
            freq1 = CP_spectral(x,y,fs,1,opt);

            opt.tap = 1;
            freq2 = CP_spectral(x,y,fs,1,opt);

            opt.tap = 2;
            freq3 = CP_spectral(x,y,fs,1,opt);
            freq4 = CP_spectral(x,y,fs,2,opt);
            freq5 = CP_spectral(x,y,fs,3,opt);
            freq6 = CP_spectral(x,y,fs,4,opt);

            disp(['trial ' num2str(tr) ' SNR ' num2str(na) ' ' datestr(clock,'HH:MM:SS')])

            for cc = 1:numel(allConnMeths)

                se = 1;
                res = CP_measure(freq1,cc);
                [perf(se,cc,na,1), perf(se,cc,na,2)] = plot_shaded(res,0,0);
                perf(se,cc,na,3) = res.etime;
                tmpRes{cc,se,na} = res;

                se = 2;
                res = CP_measure(freq2,cc);
                [perf(se,cc,na,1), perf(se,cc,na,2)] = plot_shaded(res,0,0);
                perf(se,cc,na,3) = res.etime;
                tmpRes{cc,se,na} = res;

                se = 3;
                res = CP_measure(freq3,cc);
                [perf(se,cc,na,1), perf(se,cc,na,2)] = plot_shaded(res,0,0);
                perf(se,cc,na,3) = res.etime;
                tmpRes{cc,se,na} = res;

                se = 4;
                res = CP_measure(freq4,cc);
                [perf(se,cc,na,1), perf(se,cc,na,2)] = plot_shaded(res,0,0);
                perf(se,cc,na,3) = res.etime;
                tmpRes{cc,se,na} = res;

                se = 5;
                res = CP_measure(freq5,cc);
                [perf(se,cc,na,1), perf(se,cc,na,2)] = plot_shaded(res,0,0);
                perf(se,cc,na,3) = res.etime;
                tmpRes{cc,se,na} = res;

                se = 6;
                res = CP_measure(freq6,cc);
                [perf(se,cc,na,1), perf(se,cc,na,2)] = plot_shaded(res,0,0);
                perf(se,cc,na,3) = res.etime;
                tmpRes{cc,se,na} = res;
                
            end
        end

        allPerf(:,:,:,:,nn,tr) = perf;
        allRes{nn,tr} = tmpRes;
        
    end
    
end

save([proj0012Dir 'chrisd/DI/data/2021/CPconnect/simulation_GWN_1overF.mat'], ...
    'allPerf','allRes','fs','noiseTypes','allNoiseAmps','allSpecEstMeths', ...
    'gtFilterOrder','gtPassBand','allConnMeths','nSamp','nTrials','-v7.3')
