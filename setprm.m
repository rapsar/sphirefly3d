function prm = setprm(cyber2world,v)
%SFFPARAM initializes sff parameters
%   

% real world units
prm.world.horzMtr = cyber2world(1);
try
    prm.world.vertMtr = cyber2world(2);
catch
    prm.world.vertMtr = 0.6;
end

% flags
prm.flag.trk = 0; %track
prm.flag.cln = 0; %clean
prm.flag.clb = 0; %calibrate
prm.flag.trg = 0; %triangulate
prm.flag.stk = 0; %streak
prm.flag.trj = 0; %trajectory

% track
prm.trk.bkgrStackSec = 2;
prm.trk.bwThr = 0.1;
prm.trk.blurRadiusPxl = 1;

% clean
prm.cln.initlBufferSec = 60;
prm.cln.finalBufferSec = 60;
prm.cln.flashMinBrightUI8 = 30;
prm.cln.areaOutlierThr = 4;

% triangulate
prm.clb.estMethod = 'minSearch';
prm.trg.matchPointsThrChord = 0.2;
prm.trg.distThresholdMtr = 20;

% trajectorize
prm.stk.linkRadiusMtr = 0.3; 
prm.trj.linkRadiusMtr = 1; 
prm.trj.linkMinLagSec = 0; 
prm.trj.linkMaxLagSec = 1; 

% intrinsic
prm.mov.frameWidth = v.Width;
prm.mov.frameHeight = v.Height;
prm.mov.frameDim = [v.Width v.Height];
prm.mov.frameRate = v.FrameRate; 

% frame units
prm.trk.bkgrStackFrm = round(prm.trk.bkgrStackSec*prm.mov.frameRate);
prm.cln.initlBufferFrm = round(prm.cln.initlBufferSec*prm.mov.frameRate);
prm.cln.finalBufferFrm = round(prm.cln.finalBufferSec*prm.mov.frameRate);
prm.trj.linkMinLagFrm = round(prm.trj.linkMinLagSec*prm.mov.frameRate);
prm.trj.linkMaxLagFrm = round(prm.trj.linkMaxLagSec*prm.mov.frameRate);


end

