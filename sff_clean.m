function sff = sff_clean(sff)
%SFF_CLEAN 
%   

sff.gp1.cln = ffClean(sff.gp1,sff.prm);
sff.gp2.cln = ffClean(sff.gp2,sff.prm);

end



function ffo = ffClean(ffi,prm)
% clean detected flashes by removing clusters


%% pass and keep raw data
%ffo = ffi;
ffo.n = ffi.n;
ffo.xyt = ffi.xyt;
ffo.aeit = ffi.aeit;
ffo.xy = ffi.xy;
ffo.aei = ffi.aei;


%% truncate beginning and end
N = length(ffi.n);

t = ffi.xyt(:,3);
rmt = (t <= prm.cln.initlBufferFrm | t >= (N-prm.cln.finalBufferFrm+1));

ffo.xyt(rmt,:) = [];
ffo.aeit(rmt,:) = [];


%% remove too dark
isDark = (ffo.aeit(:,3) < prm.cln.flashMinBrightUI8); % parametrize etc.

ffo.xyt(isDark,:) = [];
ffo.aeit(isDark,:) = [];

% %%%
% for i=1:N
%     f = find(ffo.xyt(:,3) == i);
%     if length(f) > 3
%         [~,d] = knnsearch(ffo.xyt(f,1:2),ffo.xyt(f,1:2),'K',2);
%         d = d(:,2);
%         ff = find(d<100);
%         ffo.xyt(f(ff),:) = [];
%         ffo.aeit(f(ff),:) = [];
%     end
% end
% 
% %%%
% 
% %% remove too frequent bins
% [n,~,~,binX,binY] = histcounts2(ffo.xyt(:,1),ffo.xyt(:,2),100,'Normalization','probability');
% [row,col] = find(n>0.001);
% isTooMany = ismember([binX binY],[row,col],'rows');
% ffo.xyt(isTooMany,:) = [];
% ffo.aeit(isTooMany,:) = [];


%% returns n, xy{} and aei{} from xyt
ffo.n = histcounts(ffo.xyt(:,3),0.5:N+0.5);
ffo.xy = mat2cell(ffo.xyt(:,1:2),ffo.n(:));
ffo.aei = mat2cell(ffo.aeit(:,1:2),ffo.n(:));


%% determine size of clusters to remove

n = ffo.n;
nn = n>0;

rp = regionprops(nn,n,'Area','MeanIntensity','PixelList');

s = vertcat(rp.Area);
i = vertcat(rp.MeanIntensity);

% chunk area (under the N-curve)
a = s.*i;

isOut = isoutlier(log(a),'mean','ThresholdFactor',prm.cln.areaOutlierThr);
tOut = vertcat(rp(isOut).PixelList);
if ~isempty(tOut)
    tOut = tOut(:,1);
end

ffo.n(tOut) = 0;
ffo.xy(tOut) = cell(length(tOut),1);
ffo.aei(tOut) = cell(length(tOut),1);
t = ffo.xyt(:,3);
ffo.xyt(ismember(t,tOut),:) = [];
ffo.aeit(ismember(t,tOut),:) = [];

% shows removed clusters
ffo.isOut = isOut;
ffo.tOut = tOut;
ffo.nOut = 0*ffo.n;
ffo.nOut(isOut) = n(isOut);

% save processing code
ffo.code = fileread([mfilename('fullpath') '.m']);

end