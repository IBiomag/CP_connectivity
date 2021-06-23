function res = CP_measure(freq,meth,opt)
%function to compute cortico-peripheral connectivity
%PLV GCMI R-test wppc Coh ent

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
        %plv
        tic
        
        % normalise the fourier coefficients to length 1 complex vectors
        freq.fourierspctrm = freq.fourierspctrm./abs(freq.fourierspctrm);
       
        x = angle(squeeze(freq.fourierspctrm(:,1,:)));
        y = angle(squeeze(freq.fourierspctrm(:,2,:)));
        
        % actual computation
        plv = abs(mean(exp(1i*(x-y))));
        
        % surrogate data
        for is = 1:nsurr
            xshifted  = circshift(x, sh(is), 1);
            tmp(is,:) = abs(mean(exp(1i*(xshifted-y))));
        end

        res.cp    = plv;
        res.surr  = tmp;
        res.meth  = 'plv';
        res.etime = toc;

    case 2   
        %gcmi
        tic
        
        % copula transform needs to be done only once, do it separately for
        % the real and imaginary values, treat as pair of bivariate time
        % series
        x1 = copnorm(real(squeeze(freq.fourierspctrm(:,1,:))));
        x2 = copnorm(imag(squeeze(freq.fourierspctrm(:,1,:))));
        y1 = copnorm(real(squeeze(freq.fourierspctrm(:,2,:))));
        y2 = copnorm(imag(squeeze(freq.fourierspctrm(:,2,:))));
        
        % actual computation
        mi = zeros(1,nf);
        for k = 1:nf
            mi(k) = mi_gg([x1(:,k) x2(:,k)],[y1(:,k) y2(:,k)],false,false);
        end
        
        % surrogate data
        for is = 1:nsurr
            x1shifted = circshift(x1, sh(is), 1);
            x2shifted = circshift(x2, sh(is), 1);
            for k = 1:nf
                tmp(is,k) = mi_gg([x1shifted(:,k) x2shifted(:,k)],[y1(:,k) y2(:,k)],false,false);
            end
        end
        
        res.cp    = mi;
        res.surr  = tmp;
        res.meth  = 'gcmi';
        res.etime = toc;

             
    case 3
        %circ_R test
        tic
        
        % normalise the fourier coefficients to length 1 complex vectors ->
        % this is not needed
        freq.fourierspctrm = freq.fourierspctrm./abs(freq.fourierspctrm);
        
        x = angle(squeeze(freq.fourierspctrm(:,1,:)));
        y = angle(squeeze(freq.fourierspctrm(:,2,:)));
        
        % actual computation
        ci = (length(x)*(abs(mean(exp(1i*(x-y)),1))).^2);
        
        % surrogate data
        for is=1:nsurr
            xshifted  = circshift(x, sh(is), 1);
            tmp(is,:)=length(x)*(abs(mean(exp(1i*(xshifted-y))))).^2;
        end
        
        res.cp    = ci;
        res.surr  = tmp;
        res.meth  = 'R-test';
        res.etime = toc;

    case 4
        %wppc
        tic
        
%         input=zeros(size(freq.fourierspctrm,1),size(freq.fourierspctrm,2),size(freq.fourierspctrm,2),size(freq.fourierspctrm,3));
%         for k1=1%:2,
%             for k2=2%1:2,
%                 input(:,k1,k2,:)=freq.fourierspctrm(:,k1,:).*conj(freq.fourierspctrm(:,k2,:));
%             end
%         end
        input    = squeeze(freq.fourierspctrm(:,1,:).*conj(freq.fourierspctrm(:,2,:)));
        outsum   = sum(input); % normalization of the WPLI
        outssq   = sum(input.*conj(input));
        outsumw  = sum(abs(input));
        ci       = (outsum.*conj(outsum) - outssq)./(outsumw.*conj(outsumw) - outssq); %
        
        for is = 1:nsurr 
            xshifted = circshift(squeeze(freq.fourierspctrm(:,1,:)), sh(is), 1);
            input    = xshifted.*squeeze(conj(freq.fourierspctrm(:,2,:)));
            outsum   = sum(input); 
            outssq   = sum(input.*conj(input));
            outsumw  = sum(abs(input));
                
            tmp(is,:) = (outsum.*conj(outsum) - outssq)./(outsumw.*conj(outsumw) - outssq); %
        end
        
        res.cp   = ci;
        res.surr = tmp;
        res.meth = 'wppc';
        res.etime = toc;
        
  case 5
        %coherence -> magnitude squared coherence
        tic
        
        x = squeeze(freq.fourierspctrm(:,1,:));
        y = squeeze(freq.fourierspctrm(:,2,:));
        
        siz   = size(freq.fourierspctrm);
        input = zeros(siz([1 2 2 3]));

        % compute the cross-spectral and auto-spectral densities
        input(:,1,1,:) = abs(x).^2;
        input(:,1,2,:) = x.*conj(y);
        input(:,2,1,:) = conj(input(:,1,2,:));
        input(:,2,2,:) = abs(y).^2;
        
        ci = squeeze(abs(mean(input(:,1,2,:),1)).^2)./squeeze(mean(abs(input(:,1,1,:)),1).*mean(abs(input(:,2,2,:)),1));

        for is=1:nsurr
            xshifted = circshift(x, sh(is), 1);
            input(:,1,2,:) = xshifted.*conj(y);
            input(:,2,1,:) = conj(input(:,1,2,:));
            % the auto terms for x are only shifted, but should not yield a
            % difference in power (i.e. the average)
            % the auto terms for y is unchanges, so does not need to be
            % recomputed
            
            co2 = squeeze(abs(mean(input(:,1,2,:),1)).^2)./squeeze(mean(abs(input(:,1,1,:)),1).*mean(abs(input(:,2,2,:)),1));
            tmp(is,:) = squeeze(co2);                
        end
        
        res.cp    = ci';
        res.surr  = tmp;
        res.meth  = 'coh';
        res.etime = toc;
      
    case 6
        %entropy
        tic
        nbins=20;lnbin=log(nbins);
        phsbins=linspace(-pi,pi,nbins);
        nsamp=size(freq.fourierspctrm,1);
        freq.fourierspctrm=(freq.fourierspctrm./abs(freq.fourierspctrm));
        x=angle(squeeze(freq.fourierspctrm(:,1,:)));
        y=angle(squeeze(freq.fourierspctrm(:,2,:)));
        phdiff=angle(exp(1i*(x-y)));
        phhist=hist(phdiff,phsbins);
        phhist=phhist/nsamp;
        tmp=phhist;
        tmp(tmp==0)=1;        %disappears anyway in sum because of 0 in phhist
        rho=-sum(phhist.*log(tmp));
        pl=(lnbin-rho)/lnbin;  %normalization of rho

        nf=size(freq.fourierspctrm,3);
        %pl=phlk_si(x,y,1,1,size(freq.fourierspctrm,1),'entropy');
        tmp=zeros(nsurr,nf);
        sh=randi([minshift size(freq.fourierspctrm,1)-minshift],1,nsurr);
        for is=1:nsurr,
            %phdiff=angle(exp(i*(circshift(x,sh(is),1)-y)));
            phdiff=mod(circshift(x,sh(is),1)-y+pi,2*pi)-pi; % the same as the above, but ~3 times as fast
            phhist=hist(phdiff,phsbins);
            phhist=phhist/nsamp;
            tmp1=phhist;
            tmp1(tmp1==0)=1;        %disappears anyway in sum because of 0 in phhist
            rho=-sum(phhist.*log(tmp1));
            tmp(is,:)=(lnbin-rho)/lnbin;  %normalization of rho
            %tmp(is,:)=phlk_si(circshift(x,sh(is),1),y,1,1,size(freq.fourierspctrm,1),'entropy');
        end
        res.f = freq.freq;
        res.cp = pl;
        res.surr = (tmp);
        res.meth = 'ent';
        res.specmeth=freq.meth;
        res.etime=toc;
        
        
        
         case 7
        %wpli
        tic
        nf=size(freq.fourierspctrm,3);
        input=zeros(size(freq.fourierspctrm,1),size(freq.fourierspctrm,2),size(freq.fourierspctrm,2),size(freq.fourierspctrm,3));
        for k1=1:2,
            for k2=1:2,
                input(:,k1,k2,:)=freq.fourierspctrm(:,k1,:).*conj(freq.fourierspctrm(:,k2,:));
            end
        end
        input    = imag(input);          % make everything imaginary
        outsum   = nansum(input,1);      % compute the sum; this is 1 x size(2:end)
        outsumW  = nansum(abs(input),1); % normalization of the WPLI
        %outssq   = nansum(input.^2,1);
        %ci     = (outsum.^2 - outssq)./(outsumW.^2 - outssq); % do the pairwise thing in a handy way
        ci     = outsum./outsumW;
        ci=squeeze(ci(1,2,1,:));
        
        tmp=zeros(nsurr,nf);
        sh=randi([minshift length(freq.cumsumcnt)-minshift],1,nsurr)*freq.cumtapcnt(1);
        for is=1:nsurr,
            x=circshift(freq.fourierspctrm,sh(is),1);
            for k=1:nf,
                for k1=1:2,
                    for k2=1:2,
                        input(:,k1,k2,:)=x(:,k1,:).*conj(freq.fourierspctrm(:,k2,:));
                    end
                end
                input    = imag(input);          % make everything imaginary
                outsum   = nansum(input,1);      % compute the sum; this is 1 x size(2:end)
                outsumW  = nansum(abs(input),1); % normalization of the WPLI                ci2     = outsum./outsumW;
                %ci2     = (outsum.^2 - outssq)./(outsumW.^2 - outssq); % do the pairwise thing in a handy way
                 ci2     = outsum./outsumW;
                tmp(is,:)=squeeze(ci2(1,2,1,:));
            end
        end
        res.f = freq.freq;
        res.cp = ci';
        res.surr = (tmp);
        res.meth = 'wpli';
        res.specmeth=freq.meth;
        res.etime=toc;

        
        

    case 8
             %distance correlation
        freq.fourierspctrm=(freq.fourierspctrm./abs(freq.fourierspctrm));
        %x=angle(squeeze(freq.fourierspctrm(:,1,:)));
        %y=angle(squeeze(freq.fourierspctrm(:,2,:)));
        
        nf=size(freq.fourierspctrm,3);
        pl=zeros(1,nf);
        for k=1:nf,
            %pl(k)=dcor_dc(x(:,k),y(:,k));
            x=[real(freq.fourierspctrm(:,1,k)) imag(freq.fourierspctrm(:,1,k))];
            y=[real(freq.fourierspctrm(:,2,k)) imag(freq.fourierspctrm(:,2,k))];
            pl(k)=dcor_dc(x,y);

        end
        tmp=zeros(nsurr,nf);
        sh=randi([minshift size(freq.fourierspctrm,1)-minshift],1,nsurr);
        for k=1:nf,
            x=[real(freq.fourierspctrm(:,1,k)) imag(freq.fourierspctrm(:,1,k))];
            y=[real(freq.fourierspctrm(:,2,k)) imag(freq.fourierspctrm(:,2,k))];
            
            for is=1:nsurr,
                %tmp(is,k)=dcor_dc(circshift(x(:,k),sh(is),1),y(:,k));
                tmp(is,k)=dcor_dc(circshift(x,sh(is),1),y);
            end
        end
        res.f = freq.freq;
        res.cp = pl;
        res.surr = (tmp);
        res.meth = 'DC';
        res.specmeth=freq.meth;
     
        
         case 9 
        %MI with binning
        %make this optional
        %freq.fourierspctrm=(freq.fourierspctrm./abs(freq.fourierspctrm));
        nb=20;
        % Mutual information (copula)
        tic
        opts.nt=[];
        opts.method = 'gs'; % dr is good;
        opts.bias   = 'naive'; %dr: qe ok, pt less noisy;naive ok
        opts.btsp = 0;
        freq.fourierspctrm=(freq.fourierspctrm./abs(freq.fourierspctrm));
        x=angle(squeeze(freq.fourierspctrm(:,1,:)));
        y=angle(squeeze(freq.fourierspctrm(:,2,:)));
        nsamp=size(freq.fourierspctrm,1);
        nf=size(freq.fourierspctrm,3);
        mi=zeros(1,nf);
        tmp=zeros(nsurr,nf);
        for k=1:nf,
            ph1=binr(x(:,k)',nsamp,nb,'eqpop');
            R=zeros(1,2,nb);opt.nt=[];
            for ib=1:nb,
                fi=find(ph1 == ib-1);
                opts.nt(ib)=length(fi);
                R(1,1:opts.nt(ib),ib)=y(fi,k);
            end
            %R2=binr(R,opts.nt',nb,'eqpop');
            I =information(R,opts,'I');
            mi(k)=I(1);
            %tmp(:,k)=I(2:end);
         end
        
        tmp=zeros(nsurr,nf);
        sh=randi([minshift length(freq.cumsumcnt)-minshift],1,nsurr)*freq.cumtapcnt(1);
        for is=1:nsurr,
            for k=1:nf,
                ph1=binr(circshift(x(:,k),sh(is),1)',nsamp,nb,'eqpop');
                R=zeros(1,2,nb);opt.nt=[];
                for ib=1:nb,
                    fi=find(ph1 == ib-1);
                    opts.nt(ib)=length(fi);
                    R(1,1:opts.nt(ib),ib)=y(fi,k);
                end
                %R2=binr(R,opts.nt',nb,'eqpop');
                tmp(is,k)=information(R,opts,'I');
            end;
        end
%    
        res.f = freq.freq;
        res.cp = (mi);
        res.surr = (tmp);
        res.meth = 'mi2';
        res.specmeth=freq.meth;
        res.etime=toc;
    case 10

         %conditional probability
        freq.fourierspctrm=(freq.fourierspctrm./abs(freq.fourierspctrm));
        x=angle(squeeze(freq.fourierspctrm(:,1,:)));
        y=angle(squeeze(freq.fourierspctrm(:,2,:)));
        
        nf=size(freq.fourierspctrm,3);
        pl=phlk_si(x,y,1,1,size(freq.fourierspctrm,1),'cprob');
        tmp=zeros(nsurr,nf);
        sh=randi([minshift size(freq.fourierspctrm,1)-minshift],1,nsurr);
        for is=1:nsurr,
            %ri=randperm(nsamp);
            tmp(is,:)=phlk_si(circshift(x,sh(is),1),y,1,1,size(freq.fourierspctrm,1),'cprob');
        end
        res.f = freq.freq;
        res.cp = pl;
        res.surr = (tmp);
        res.meth = 'SI (entropy)';
        res.specmeth=freq.meth;
        
    case 11 
        %entropy 2
         nb=8;
        % Mutual information (copula)
        tic
        opts.nt=[];
        opts.method = 'gs'; % dr is good;
        opts.bias   = 'naive'; %dr: qe ok, pt less noisy;naive ok
        %opts.btsp = 20;
        freq.fourierspctrm=(freq.fourierspctrm./abs(freq.fourierspctrm));
        x=angle(squeeze(freq.fourierspctrm(:,1,:)));
        y=angle(squeeze(freq.fourierspctrm(:,2,:)));
        phdiff=angle(exp(i*(x-y)));

        nsamp=size(freq.fourierspctrm,1);
        nf=size(freq.fourierspctrm,3);
        mi=zeros(1,nf);
        R2=zeros(1,1,nsamp);
        opts.nt=1;
        for k=1:nf,
           % R(1,1,:)=phdiff(:,k);
            R=binr(phdiff(:,k)',nsamp,nb,'eqspace',[-pi pi]);
            %R2(1,1,:)=R;
            mi(k)=entropy(R,opts,'HR');
         end
        
        tmp=zeros(nsurr,nf);
        sh=randi([minshift length(freq.cumsumcnt)-minshift],1,nsurr)*freq.cumtapcnt(1);
        for is=1:nsurr,
        phdiff=angle(exp(i*(circshift(x,sh(is),1)-y)));
            for k=1:nf,
            %R(1,1,:)=phdiff(:,k);
             R=binr(phdiff(:,k)',nsamp,nb,'eqspace',[-pi pi]);
                tmp(is,k)=entropy(R,opts,'HR');
            end;
        end
   
        res.f = freq.freq;
        res.cp = (mi);
        res.surr = (tmp);
        res.meth = 'entropy';
        res.specmeth=freq.meth;
        res.etime=toc;
end
end


