function sff = sff_estimate360CameraParameters(sff)
%SFF_ESTIMATE360CAMERAPARAMETERS 
%   
% Raphael Sarfati, 05/2022

calAlpha1 = xy2alpha(sff.clb.calPoints.p1,sff.prm.mov.frameDim);
calAlpha2 = xy2alpha(sff.clb.calPoints.p2,sff.prm.mov.frameDim);

sff.clb.stereo360ParamsRANSAC = estimate360CameraParameters(calAlpha1(:,1:3),calAlpha2(:,1:3),'RANSAC');
sff.clb.stereo360Params = estimate360CameraParameters(calAlpha1(:,1:3),calAlpha2(:,1:3),sff.prm.clb.estMethod);
%sff.stereo360Params.t = -sff.stereo360Params.t; %%%%%%%%%%%%%%% why
end


%% estimateCameraParameters
function stereo360Params = estimate360CameraParameters(matchedAlpha1,matchedAlpha2,estMethod)
%ESTIMATE360CAMERAPARAMETERS Estimates camera pose (t,R) from set of matched
%spherical projections alpha.
% Estimation methods estMethod are: 
%   - RANSAC, MSAC, LTS, Norm8Point, LMedS for Fundamental Matrix
%   - or minSearch for optimization search
% (see corresponding functions for details)
%
%   
% Raphael Sarfati, 03/2020
% raphael.sarfati@aya.yale.edu

% default method is fundamental matrix estimation using RANSAC
if nargin == 2
    estMethod = 'RANSAC';
end

if strcmp(estMethod,'RANSAC')
    stereo360Params = tRestimate_fundamentalMatrix(matchedAlpha1,matchedAlpha2,'RANSAC');
    
elseif strcmp(estMethod,'MSAC')
    stereo360Params = tRestimate_fundamentalMatrix(matchedAlpha1,matchedAlpha2,'MSAC');
    
elseif strcmp(estMethod,'LTS')
    stereo360Params = tRestimate_fundamentalMatrix(matchedAlpha1,matchedAlpha2,'LTS');
    
elseif strcmp(estMethod,'Norm8Point')
    stereo360Params = tRestimate_fundamentalMatrix(matchedAlpha1,matchedAlpha2,'Norm8Point');
    
elseif strcmp(estMethod,'LMedS')
    stereo360Params = tRestimate_fundamentalMatrix(matchedAlpha1,matchedAlpha2,'LMedS');
    
elseif strcmp(estMethod,'minSearch')
    stereo360Params = tRestimate_minSearch(matchedAlpha1,matchedAlpha2);
    
end

end


%% tRestimate_fundamentalMatrix
function stereo360Params = tRestimate_fundamentalMatrix(alpha1,alpha2,estMethod)
%TRESTIMATE_FUNDAMENTALMATRIX Estimates t and R from calculation of F.
%   Uses built-in function estimateFundamentalMatrix to estimate F.
%   Calculates t and R from F.
%   Different algorithms can be used to estimate F.
%
% Raphael Sarfati, 03/2020
% raphael.sarfati@aya.yale.edu


% set the default estimation method to RANSAC
% other options are 'Norm8Point', 'LMedS', 'MSAC', 'LTS'
% see doc for more details
if nargin == 2
    estMethod = 'RANSAC';
end

% renormalize to get third coordinate equal to 1 (see Matlab documentation)
points1 = alpha1(:,1:2)./alpha1(:,3);
points2 = alpha2(:,1:2)./alpha2(:,3);

% calculate fundamental matrix
% for consistency, points of camera 2 enter as first argument
F = estimateFundamentalMatrix(points2,points1,'Method',estMethod);

% estimate t and R from F
tR = F2tR(F);

stereo360Params.t = tR.t;
stereo360Params.R = tR.R;
stereo360Params.F = F;
stereo360Params.Method = estMethod;

end


%% F2tR
function [out] = F2tR(F,tguess,Rguess)
%F2TR Returns best estimate for t and R, given F (= Tx*R')
%   F is the fundamental matrix, given by estimateFundamentalMatrix
%   tguess, Rguess are estimates useful to resolve ambiguities on t, R
%   Method is based on single-value decomposition (SVD) and W, Z matrices.
%   See link: http://www.maths.lth.se/matematiklth/personal/calle/datorseende13/notes/forelas6.pdf
%
% RS, 10/2019

% if no estimate for t, R, assume t = tx and R = Id.
if nargin == 1
    tguess = [1 0 0];
    Rguess = eye(3);
elseif nargin == 2
    Rguess = eye(3);
end

% make tguess vertical
tguess = tguess(:);

% Decomposition using SVD and W, Z as explained in link above.
[U,~,V] = svd(F);

W = [0 -1  0 ; 
     1  0  0 ;
     0  0  1]; 
 
Z = [0 1 0; 
    -1 0 0; 
     0 0 0];
 
S1 = -U*Z*U';
S2 = U*Z*U';

R1 = U*W'*V';
R2 = U*W*V';

e1 = null(S1);
out.e1 = e1;

% corrects to get real rotations
if det(R1)<0
    out.R1 = -R1';
else
    out.R1 = R1';
end

if det(R2)<0
    out.R2 = -R2';
else
    out.R2 = R2';
end

if norm(out.R1-Rguess) < norm(out.R2-Rguess)
    out.R = out.R1;
else
    out.R = out.R2;
end

if norm(e1-tguess) < norm(-e1-tguess)
    out.t = e1;
else
    out.t = -e1;
end


end


%% tRestimate_minSearch
function [out] = tRestimate_minSearch(alpha1,alpha2,x0)
%TRESTIMATE_MINSEARCH Calculates t and R by optimization. 
%   Minimizes sum of distances given t and R, represented as a 6D vector
%   x = {tx,ty,tz,phix,phiy,phiz}, with phi the 3 rotation angles.
%   alpha1, alpha2 are matching sets of projected points
%   x0 (optional) is the initial guess vector
%
%   Uses fmincon, since we have the constraint |t| = 1.
%   The cost function is defined below.
%
% Raphael Sarfati, 03/2020

% if guess vector not provided, assumes t = tx and R = Id
if nargin == 2
    x0 = [0 1 0 0 0 0]';
end

% arguments for fmincon
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon = @fcon;
options = optimset('MaxFunEvals',1e6,'MaxIter',1e6,'Display','off','PlotFcn',{@optimplotx,@optimplotfval});

% optimization, applying fmincon to the cost function defined separately
% and with the constraint |t| = 1
xmin = fmincon(@(x) tRcostFunction(alpha1,alpha2,x),x0,...
    A,b,Aeq,beq,lb,ub,nonlcon,options);


out.t = xmin(1:3);
out.phi = xmin(4:6);
out.R = rotx(xmin(4))*roty(xmin(5))*rotz(xmin(6));
out.F = tR2F(out.t,out.R);

end

function [c,ceq] = fcon(x)
%FCON Constraint for fmincon.
%   t(1)^2 + t(2)^2 + t(3)^2 = 1

c = [];
ceq = x(1)^2 + x(2)^2 +x(3)^2 - 1;

end


%% Cost Function

function [c] = tRcostFunction(alpha1,alpha2,x)
%TRCOSTFUNCTION Cost (sum of squared distances) given t and R.
%   Calculates sum of squared distances between points in camera 1 and
%   camera 2 given the 6 parameters in x (translation and rotation
%   coordinates).
%   Given the projected points alpha1, alpha2, and a test t, R, there is an
%   optimum set (r1,r2) to mimimize || r1*a1 - (t + R'*r2*a2) ||, 
%   so (r1,r2) don't need to be put into the test vector x.
%
% RS, 10/2019

% number of points
N = size(alpha1,1);

% t and R from test vector x
t = [x(1) x(2) x(3)];
R = rotx(x(4))*roty(x(5))*rotz(x(6));

% calculate optimal r1 and r2
r1 = zeros(N,1);
r2 = zeros(N,1);

for i=1:N
    
    rr = rOptimum(alpha1(i,:),alpha2(i,:),t,R);
    r1(i) = rr(1);
    r2(i) = rr(2);
    
end

% separations between X1 and t+R'*X2 
s = r1.*alpha1 - r2.*(R'*alpha2')' - t;

% sum of distances
%c = sum(vecnorm(s,2,2),1);
c = mean(vecnorm(s,2,2),1);


end


function [rr] = rOptimum(a1,a2,t,R)
%ROPTIMUM Calculates r1 and r2 that minimize N = || r1*a1 - t - R'*r2*a2 ||
% Matrix method based on Ma 2015, but equivalent to setting dN/dr1 = 0 and
% dN/dr2 = 0.

a1 = a1(:);
a2 = a2(:);
t = t(:);

C = [a1 , -R'*a2];

%rr = (C'*C)\C'*t;

A = -eye(2);
b = zeros(2,1);
Aeq = [];
beq = [];
lb = [];
ub = [];
x0 = [];
options = optimoptions(@lsqlin,'Display','off');
rr = lsqlin(C,t,A,b,Aeq,beq,lb,ub,x0,options);

end

%% tR2F
function [F] = tR2F(t,R)
%TR2F Calculates the fundamental matrix F from t and R
%   t is a vector 
%   R should be rotation matrix but can also be a vector of
%   the three rotation angles around the x, y, and z axes respectively.
%   NB: rotx, roty, rotz require Phased Array System Toolbox; if unable to
%   get, uncomment function at the end of this file.
%
% RS, 10/2019

Tx = vcross(t);

if ismatrix(R)
    F = Tx*R';
    
elseif isvector(R)
    phixyz = R;
    R = rotx(phixyz(1))*roty(phixyz(2))*rotz(phixyz(3));
    F = Tx*R';
    
else
    error('R must be either a rotation matrix or a vector of three angles.');
    
end

end

function [Vx] = vcross(v)
%VCROSS Calculate matrix Vx such that Vx*u = cross(v,u).

Vx = [0 -v(3) v(2) ; v(3) 0 -v(1); -v(2) v(1) 0]; 

end


%% Manual definition of rotx, roty, rotz if unable to get Phased Array System Toolbox 
% Uncomment if necessary.
% See https://www.mathworks.com/help/phased/ref/rotx.html for details.
% 
function Rx = rotx(phi)

Rx = [1 0 0 ; 0 cos(phi) -sin(phi) ; 0 sin(phi) cos(phi)];

end

function Ry = roty(phi)

Ry = [cos(phi) 0 sin(phi) ; 0 1 0 ; -sin(phi) 0 cos(phi)];

end

function Rz = rotz(phi)

Rz = [cos(phi) -sin(phi) 0 ; sin(phi) cos(phi) 0 ; 0 0 1];

end

%%
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

