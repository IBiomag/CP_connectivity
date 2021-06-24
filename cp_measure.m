function res = cp_measure(freq, meth, opt)

% cp_measure computes 'cortico-peripheral' connectivity from a fieldtrip
% style freq structure, consisting of 2 channels, and formatted as if
% cfg.output = 'fourier', was used in ft_freqanalysis. Although most of the
% connectivity metrics are implemented, and accessible using
% ft_connectivityanalysis, the below code uses its own implementation for
% speed considerations, since the overhead in data bookkeeping is bypassed.
% Next to computing the connectivity metric, 200 surrogates are computed,
% based on a circular shift (across time) for one of the channels, prior to
% recomputation of the connectivity.
%
% Use as
%
%  res = cp_measure(freq, meth)
%
% with input arguments
%  freq = fieldtrip style structure containing fourier coefficients for 2
%         channels
%  meth = scalar number indicating the metric to be computed
%
%  1. PLV 
%  2. GCMI
%  3. R-test
%  4. WPPC
%  5. magnitude-squared coherence
%  6. entropy

if nargin<3
    opt = struct([]);
end

opt(1).uo = ft_getopt(opt, 'uo',   'taper');
if ismember(meth, [1 3 6])
    opt.norm = ft_getopt(opt, 'norm', opt.uo); % original behavior
    assert(isequal(opt.norm, opt.uo));
elseif ismember(meth, [2 4 5])
    opt.norm = ft_getopt(opt, 'norm', 'no');
elseif ismember(meth, [7:11])
    error('to be implemented');
end

nsurr    = 200; %number of randomizations, hard coded
minshift = round(0.1*size(freq.cumsumcnt,1));
sh       = randi([minshift length(freq.cumsumcnt)-minshift],1,nsurr)*freq.cumtapcnt(1);
nf       = size(freq.fourierspctrm,3);
tmp      = zeros(nsurr,nf);

% start creating the output structure
res          = [];
res.f        = freq.freq(:)';
res.specmeth = freq.meth;

switch meth
    case 1
        connfun  = @plv;
         
    case 2
        connfun = @gcmi;
        
        % copula transform is best done here, to avoid redundant
        % computations
        freq.fourierspctrm(:,1,:) = copnorm(squeeze(real(freq.fourierspctrm(:,1,:)))) + 1i.*copnorm(squeeze(imag(freq.fourierspctrm(:,1,:))));
        freq.fourierspctrm(:,2,:) = copnorm(squeeze(real(freq.fourierspctrm(:,2,:)))) + 1i.*copnorm(squeeze(imag(freq.fourierspctrm(:,2,:))));
        
    case 3
        connfun = @rtest;
    case 4
        connfun = @wppc;
    case 5
        connfun = @coh;
    case 6
        connfun = @ent;
    case 7
    case 8
    case 9
    case 10
    case 11
    case 12
    otherwise
end

switch opt.norm
    case 'no'
        % no normalisation
    case 'taper'
        % normalise each taper
        freq.fourierspctrm = freq.fourierspctrm./abs(freq.fourierspctrm);
        
    case 'trial'
        % normalise each trial
        powtap = abs(freq.fourierspctrm).^2;
        
        assert(all(freq.cumtapcnt(:)==freq.cumtapcnt(1)));
        xindx = repmat((1:numel(freq.cumtapcnt)), freq.cumtapcnt(1), 1);
        xindx = xindx(:);
        yindx = (1:size(freq.fourierspctrm,1))';
        zindx = ones(numel(yindx),1)./freq.cumtapcnt(1);
        P     = sparse(xindx, yindx, zindx); % fast averaging matrix
        pow1(1,:,:) = P * squeeze(powtap(:,1,:));
        pow2(1,:,:) = P * squeeze(powtap(:,2,:));
        
        pow(:,1,:) = reshape(repmat(pow1, freq.cumtapcnt(1), 1), [], numel(freq.freq));
        pow(:,2,:) = reshape(repmat(pow2, freq.cumtapcnt(1), 1), [], numel(freq.freq));
        
        freq.fourierspctrm = freq.fourierspctrm./sqrt(pow);
    otherwise
        error('unknown normalisation option');
        
end
args = {freq.fourierspctrm};

switch opt.uo
    case 'taper'
        % original implementation, treat each taper as a unit of
        % observation
    case 'trial'
        % treat each input epoch as a unit of observation, i.e. average the
        % 'csd' across tapers before computing the actual metric
        assert(all(freq.cumtapcnt(:)==freq.cumtapcnt(1)));
        xindx = repmat((1:numel(freq.cumtapcnt)), freq.cumtapcnt(1), 1);
        xindx = xindx(:);
        yindx = (1:size(freq.fourierspctrm,1))';
        zindx = ones(numel(yindx),1)./freq.cumtapcnt(1);
        P     = sparse(xindx, yindx, zindx); % fast averaging matrix

        args = cat(2, args, {P});
    otherwise
        error('unknown unit of observation option');
end


tic;
res.cp = connfun(args{:});
for k = 1:nsurr
    args{1}(:,1,:) = circshift(squeeze(freq.fourierspctrm(:,1,:)), sh(k), 1);
    tmp(k,:)       = connfun(args{:});
end
res.etime = toc;
res.surr  = tmp;
res.meth  = func2str(connfun);


%%%%%%%%%%%%%%%%%%%%%%%
function out = plv(dat, P)
        
% plv
          
% get the phase differences
phdiff = squeeze(dat(:,1,:).*conj(dat(:,2,:)));

if nargin>1
    phdiff = P*phdiff;
end

% actual computation
out = abs(mean(phdiff));

%%%%%%%%%%%%%%%%%%%%%%%%
function out = gcmi(dat)

% gcmi

% copula transform needs to be done only once, do it separately for the 
% real and imaginary values, treat as pair of bivariate time series
x1 = real(squeeze(dat(:,1,:)));
x2 = imag(squeeze(dat(:,1,:)));
y1 = real(squeeze(dat(:,2,:)));
y2 = imag(squeeze(dat(:,2,:)));
         
% actual computation
nf  = size(dat,3);
out = zeros(1,nf);
for k = 1:nf
    out(k) = mi_gg([x1(:,k) x2(:,k)], [y1(:,k) y2(:,k)], false, false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
function out = rtest(dat, P)

% circular R-test

% get the phase differences
phdiff = squeeze(dat(:,1,:).*conj(dat(:,2,:)));
  
if nargin>1
    phdiff = P*phdiff;
end

% actual computation
out = size(phdiff,1)*(abs(mean(phdiff))).^2;
 
%%%%%%%%%%%%%%%%%%%%%%%%
function out = wppc(dat, P)

% wpcc
input    = squeeze(dat(:,1,:).*conj(dat(:,2,:)));

if nargin>1
    input = P*input;
end
outsum   = sum(input); 
outssq   = sum(input.*conj(input));
outsumw  = sum(abs(input));
out      = (outsum.*conj(outsum) - outssq)./(outsumw.*conj(outsumw) - outssq); %

%%%%%%%%%%%%%%%%%%%%%%%
function out = coh(dat, P)

% coherence -> magnitude squared coherence

x = squeeze(dat(:,1,:));
y = squeeze(dat(:,2,:));

% compute the cross-spectral and auto-spectral densities
xy = abs(mean(x.*conj(y)));
px = mean(abs(x).^2);
py = mean(abs(y).^2);

if nargin>1
    xy = P*xy;
    px = P*x;
    py = P*y;
end


out = (xy.^2)./(px.*py);

%%%%%%%%%%%%%%%%%%%%%%%
function out = ent(dat, P)

% entropy

nbins   = 20;
lnbin   = log(nbins);
phsbins = linspace(-pi,pi,nbins);
nsamp   = size(dat, 1);
         
% get the phase differences
phdiff = angle(squeeze(dat(:,1,:).*conj(dat(:,2,:))));
   
if nargin>1
    phdiff = P*phdiff;
end

%phdiff = angle(exp(1i*(x-y)));
%phdiff = mod(x - y + pi, 2*pi)-pi; % the same as above, but ~3 times as fast
phhist = hist(phdiff, phsbins);
phhist = phhist/nsamp;
         
tmp         = phhist;
tmp(tmp==0) = 1; %disappears anyway in sum because of 0 in phhist
rho         = -sum(phhist.*log(tmp));
out         = (lnbin-rho)/lnbin;  %normalization of rho

%     case 7
%         %wpli
%         tic
%         
%         input    = squeeze(freq.fourierspctrm(:,1,:).*conj(freq.fourierspctrm(:,2,:)));
%         input    = imag(input);       % take the imaginary part
%         outsum   = sum(input,1);      % compute the sum; this is 1 x size(2:end)
%         outsumW  = sum(abs(input),1); % normalization of the WPLI
%         ci       = outsum./outsumW;
%         
%         for is = 1:nsurr 
%           xshifted = circshift(squeeze(freq.fourierspctrm(:,1,:)), sh(is), 1);
%           input    = xshifted.*squeeze(conj(freq.fourierspctrm(:,2,:)));
%           
%           input    = imag(input);
%           outsum   = nansum(input,1);
%           outsumW  = nansum(abs(input),1);
%           tmp(is,:) = outsum./outsumW;
%         end
%         
%         res.cp    = ci';
%         res.surr  = tmp;
%         res.meth  = 'wpli';
%         res.etime = toc;
% 
%     case 8
%              %distance correlation
%         freq.fourierspctrm=(freq.fourierspctrm./abs(freq.fourierspctrm));
%         %x=angle(squeeze(freq.fourierspctrm(:,1,:)));
%         %y=angle(squeeze(freq.fourierspctrm(:,2,:)));
%         
%         nf=size(freq.fourierspctrm,3);
%         pl=zeros(1,nf);
%         for k=1:nf
%             %pl(k)=dcor_dc(x(:,k),y(:,k));
%             x=[real(freq.fourierspctrm(:,1,k)) imag(freq.fourierspctrm(:,1,k))];
%             y=[real(freq.fourierspctrm(:,2,k)) imag(freq.fourierspctrm(:,2,k))];
%             pl(k)=dcor_dc(x,y);
% 
%         end
%         tmp=zeros(nsurr,nf);
%         sh=randi([minshift size(freq.fourierspctrm,1)-minshift],1,nsurr);
%         for k=1:nf 
%             x = [real(freq.fourierspctrm(:,1,k)) imag(freq.fourierspctrm(:,1,k))];
%             y = [real(freq.fourierspctrm(:,2,k)) imag(freq.fourierspctrm(:,2,k))];
%             
%             for is=1:nsurr 
%                 %tmp(is,k)=dcor_dc(circshift(x(:,k),sh(is),1),y(:,k));
%                 tmp(is,k)=dcor_dc(circshift(x,sh(is),1),y);
%             end
%         end
%         
%         res.cp   = pl;
%         res.surr = (tmp);
%         res.meth = 'DC';
%         
%     case 9 
%         % MI with binning requires IBTB toolbox
%         nb=20;
%         
%         tic
%         opts.nt=[];
%         opts.method = 'gs'; % dr is good;
%         opts.bias   = 'naive'; %dr: qe ok, pt less noisy;naive ok
%         opts.btsp = 0;
%         freq.fourierspctrm=(freq.fourierspctrm./abs(freq.fourierspctrm));
%         x=angle(squeeze(freq.fourierspctrm(:,1,:)));
%         y=angle(squeeze(freq.fourierspctrm(:,2,:)));
%         nsamp=size(freq.fourierspctrm,1);
%         nf=size(freq.fourierspctrm,3);
%         mi=zeros(1,nf);
%         tmp=zeros(nsurr,nf);
%         for k=1:nf,
%             ph1=binr(x(:,k)',nsamp,nb,'eqpop');
%             R=zeros(1,2,nb);opt.nt=[];
%             for ib=1:nb,
%                 fi=find(ph1 == ib-1);
%                 opts.nt(ib)=length(fi);
%                 R(1,1:opts.nt(ib),ib)=y(fi,k);
%             end
%             %R2=binr(R,opts.nt',nb,'eqpop');
%             I =information(R,opts,'I');
%             mi(k)=I(1);
%             %tmp(:,k)=I(2:end);
%          end
%         
%         tmp=zeros(nsurr,nf);
%         sh=randi([minshift length(freq.cumsumcnt)-minshift],1,nsurr)*freq.cumtapcnt(1);
%         for is=1:nsurr,
%             for k=1:nf,
%                 ph1=binr(circshift(x(:,k),sh(is),1)',nsamp,nb,'eqpop');
%                 R=zeros(1,2,nb);opt.nt=[];
%                 for ib=1:nb,
%                     fi=find(ph1 == ib-1);
%                     opts.nt(ib)=length(fi);
%                     R(1,1:opts.nt(ib),ib)=y(fi,k);
%                 end
%                 %R2=binr(R,opts.nt',nb,'eqpop');
%                 tmp(is,k)=information(R,opts,'I');
%             end;
%         end
% %    
%         res.f = freq.freq;
%         res.cp = (mi);
%         res.surr = (tmp);
%         res.meth = 'mi2';
%         res.specmeth=freq.meth;
%         res.etime=toc;
%     case 10
% 
%          %conditional probability
%         freq.fourierspctrm=(freq.fourierspctrm./abs(freq.fourierspctrm));
%         x=angle(squeeze(freq.fourierspctrm(:,1,:)));
%         y=angle(squeeze(freq.fourierspctrm(:,2,:)));
%         
%         nf=size(freq.fourierspctrm,3);
%         pl=phlk_si(x,y,1,1,size(freq.fourierspctrm,1),'cprob');
%         tmp=zeros(nsurr,nf);
%         sh=randi([minshift size(freq.fourierspctrm,1)-minshift],1,nsurr);
%         for is=1:nsurr,
%             %ri=randperm(nsamp);
%             tmp(is,:)=phlk_si(circshift(x,sh(is),1),y,1,1,size(freq.fourierspctrm,1),'cprob');
%         end
%         res.f = freq.freq;
%         res.cp = pl;
%         res.surr = (tmp);
%         res.meth = 'SI (entropy)';
%         res.specmeth=freq.meth;
%         
%     case 11 
%         %entropy 2
%          nb=8;
%         % Mutual information (copula)
%         tic
%         opts.nt=[];
%         opts.method = 'gs'; % dr is good;
%         opts.bias   = 'naive'; %dr: qe ok, pt less noisy;naive ok
%         %opts.btsp = 20;
%         freq.fourierspctrm=(freq.fourierspctrm./abs(freq.fourierspctrm));
%         x=angle(squeeze(freq.fourierspctrm(:,1,:)));
%         y=angle(squeeze(freq.fourierspctrm(:,2,:)));
%         phdiff=angle(exp(i*(x-y)));
% 
%         nsamp=size(freq.fourierspctrm,1);
%         nf=size(freq.fourierspctrm,3);
%         mi=zeros(1,nf);
%         R2=zeros(1,1,nsamp);
%         opts.nt=1;
%         for k=1:nf,
%            % R(1,1,:)=phdiff(:,k);
%             R=binr(phdiff(:,k)',nsamp,nb,'eqspace',[-pi pi]);
%             %R2(1,1,:)=R;
%             mi(k)=entropy(R,opts,'HR');
%          end
%         
%         tmp=zeros(nsurr,nf);
%         sh=randi([minshift length(freq.cumsumcnt)-minshift],1,nsurr)*freq.cumtapcnt(1);
%         for is=1:nsurr,
%         phdiff=angle(exp(i*(circshift(x,sh(is),1)-y)));
%             for k=1:nf,
%             %R(1,1,:)=phdiff(:,k);
%              R=binr(phdiff(:,k)',nsamp,nb,'eqspace',[-pi pi]);
%                 tmp(is,k)=entropy(R,opts,'HR');
%             end;
%         end
%    
%         res.f = freq.freq;
%         res.cp = (mi);
%         res.surr = (tmp);
%         res.meth = 'entropy';
%         res.specmeth=freq.meth;
%         res.etime=toc;
% 
% end
% 
% 
