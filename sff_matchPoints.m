function sff = sff_matchPoints(sff)
%SFF_MATCHPOINTS 
%   
% Raphael Sarfati, 05/2022

xyt1 = sff.gp1.cln.xyt;
xyt2 = sff.gp2.cln.xyt;
dk = sff.clb.dk;
frameDim = sff.prm.mov.frameDim;
stereo360Params = sff.clb.stereo360Params;
distThresh = sff.prm.trg.matchPointsThrChord;

[matchedAlpha1t,matchedAlpha2t] = matchPoints(xyt1,xyt2,dk,frameDim,stereo360Params,distThresh);

sff.trg.matchedAlpha1t = matchedAlpha1t;
sff.trg.matchedAlpha2t = matchedAlpha2t;


end

function [matchedAlpha1t,matchedAlpha2t] = matchPoints(xyt1,xyt2,dk,frameDim,stereo360Params,distThresh)
%MATCHPOINTS Matches pairs of x,y coordinates in equirectangular frames given t and R.
%
%   xyt1 (N1x3) and xyt2 (N2x3) contain list of xy coordinates and
%   corresponding frame number t; they need not have same size
%   v informs about frame size, see xy2alpha
%   stereo360Params contains t and R
%   dThreshold is cost threshold, see matchAlphaPoints 
%   
%
% Raphael Sarfati, 03/2020

xyt2(:,3) = xyt2(:,3)-dk; % verify THAT -- could be a plus

matchedAlpha1t = [];
matchedAlpha2t = [];

for t = unique(xyt1(:,3))'
    alpha1 = xy2alpha(xyt1(xyt1(:,3)==t,1:2),frameDim);
    alpha2 = xy2alpha(xyt2(xyt2(:,3)==t,1:2),frameDim);
    
    if ~isempty(alpha1) && ~isempty(alpha2)
        
        [ma1,ma2] = matchAlphaPoints(alpha1,alpha2,stereo360Params,distThresh);
        matchedAlpha1t = vertcat(matchedAlpha1t,[ma1 repmat(t,size(ma1,1),1)]);
        matchedAlpha2t = vertcat(matchedAlpha2t,[ma2 repmat(t,size(ma1,1),1)]);
        
    end
end



end

function [alpha] = xy2alpha(xy,v)
%XY2ALPHA Converts the xy (Nx2) coordinates in the equirectangular frame to
%   spherical projections alpha (Nx3)
%
%   xy is list of positions (Nx(2+p), subsequent p columns added at the end)
%   v informs about size of equirectangular frames; it is either: 
%       - a movie structure 
%       - an array containing [width height]
%
%   alpha is Nx(3+p) array (alpha projections are first 3 columns)
%   
%
% Raphael Sarfati, 03/2020

[theta,phi] = xy2thetaphi(xy,v);

alpha = horzcat(cos(theta).*sin(phi),sin(theta).*sin(phi),cos(phi),xy(:,3:end));

end

function [theta,phi] = xy2thetaphi(xy,v)
%XY2THETAPHI Converts the xy coordinates in the equirectangular frame to
%   polar angle theta and azimuthal angle phi
%
%   xy is a list of positions (Nx2, subsequent colums are not considered)
%   v informs about size of equirectangular frames; it is either: 
%       - a movie structure 
%       - an array containing [width height]
%   
%
% Raphael Sarfati, 03/2020

if nargin == 1
    error('You need to specify the size of equirectangular frames.')
else
    if isobject(v)
        w = v.Width;
        h = v.Height;
    elseif isnumeric(v)
        w = v(1);
        h = v(2);
        if w~=2*h
            error('Frames dimensions are not consistent; width must be twice height.')
        end
    end
end

theta = xy(:,1)*(2*pi)/w;
phi = xy(:,2)*pi/h;

end

function [matchedAlpha1,matchedAlpha2,m,cij] = matchAlphaPoints(alpha1,alpha2,stereo360Params,dThreshold)
%MATCHALPHAPOINTS Finds optimal overall matching of points given t and R.
%
%   alpha1, alpha2 (Nx3) are the spherical projections in a given frame
%   stereo360Params contains t and R
%   dThreshold is maximum allowed cost for pairing (recommended: 0.1)
%
%   Cost Function
%   Function calculates the optimal distances r1, r2 assuming points i and j are
%   matched, and calculates the distance between the two reconstituted
%   points.
%   Then applies matching algorithm.
%
% Raphael Sarfati 03/2020
% rev 05/2022

% default value for dThreshold
if nargin == 3
    dThreshold = 0.1;
end

% initialization
t = stereo360Params.t(:);
R = stereo360Params.R;
N1 = size(alpha1,1);
N2 = size(alpha2,1);
beta2 = (R'*alpha2')';
cij = NaN(N1,N2);

% calculates cij for each pair ij
for i = 1:N1    
    for j = 1:N2
        
        % discounts points they are too close to the cameras' connecting
        % line
        if abs(dot(alpha1(i,:)',t))>0.99 || abs(dot(beta2(j,:)',t))>0.99
            cij(i,j) = Inf;
            
        % Calculates the optimal distances r1, r2 assuming points i and j are
        % matched, and calculates the distance between the two reconstituted
        % points.
        else
            C = [alpha1(i,:)' -beta2(j,:)'];
            r1r2 = lsqnonneg(C,t);
            cij(i,j) = vecnorm(C*r1r2-t);
            
        end
        
    end
end

% apply hungarian/munkres assignment algorithm;
% can be replaced with Matlab FEX's munkres function if old Matlab version
% https://www.mathworks.com/matlabcentral/fileexchange/20328-munkres-assignment-algorithm
m = matchpairs(cij,dThreshold);

matchedAlpha1 = alpha1(m(:,1),:);
matchedAlpha2 = alpha2(m(:,2),:);

end



