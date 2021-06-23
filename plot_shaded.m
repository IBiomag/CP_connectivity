function [auc,dVal]=plot_shaded(res,p,realdat)
%compute boostrap Confidence Interval and plot (if p=1;
%realdat=1 for real data and 0 for simulation

fw = res.f;
mi = res.cp;
surr = res.surr;

nf = length(fw);
sm = mean(surr,1);
sd = std(surr,0,1);
bc = bootci(500,@(x)[mean(x) std(x)],surr);

% z-score relative to surrogate distribution
z = (mi-sm)./sd;

% (x - lo_mean(surrogate)) / lo_std(surrogate)? = upper bound of z?
upperBound = (mi-bc(1,1:nf))./bc(1,nf+1:end);
% (x - hi_mean(surrogate)) / hi_std(surrogate)? = lower bound of z?
lowerBound = (mi-bc(2,1:nf))./bc(2,nf+1:end);
% this eventually works with shadederror, but more or less by accident: this
% makes both bounds negative, and shadedErrorBar thus adds the negative
% lower bound and subtracts the negative upper bound; better code would be
% to have confl = [upperBound-z z-lowerBound]
confl = [(lowerBound-z); (z-upperBound)]';

% I would suggest to change this to prctile(max(surr,[],2),99), because the
% way it is right now there is no control for multiple comparisons
thresh = prctile(surr,99,1);
thresh = (thresh-sm)./sd;
thresh = mean(thresh); %standardize

p1 = zeros(size(fw));
p1(fw>=3 & fw<=6) = 1;

p2 = z-thresh;
p2(p2<0) = 0; % ??

p2b = p2;
p2b(p2<0) = NaN; % after setting everything in p2 to 0 that is  0, this line is superfluous?

if realdat==1
    dVal = nanmean(p2b); % "average distance from 99th percentile" -- I don't see why the negative values aren't included here. Wouldn't it make more sense to e.g. count the area >thresh & < z-curve?
    auc = 0;
else
    dVal = nanmean(p2b(fw>=3 & fw<=6)); %average distance from 99th percentile
    [x,y,t,auc] = perfcurve(p1,p2,1); %check if correct
end

if p==1
    
    shadedErrorBar((fw)',(z)',(confl));
    hold;
    plot(fw,thresh*ones(size(fw)),'k--');
    hold
    
    if realdat==1
        title([res.specmeth ' ' res.meth ' D: ' num2str(dVal,2)],'Interpreter','none');
    else
        title([res.specmeth ' ' res.meth ' AUC: ' num2str(auc,2) ' D: ' num2str(dVal,2)],'Interpreter','none');
    end
end


