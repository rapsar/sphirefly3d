function sff = sff_extractCalibrationTrajectories(sff)
%SFF_EXTRACTCALIBRATIONTRAJECTORIES 
%   
% RS, 05/2022

xyt1 = sff.gp1.cln.xyt;
xyt2 = sff.gp2.cln.xyt;
dk = sff.clb.dk;

sff.clb.calPoints = extractCalibrationPoints(xyt1,xyt2,dk);
sff.clb.calPoints = spreadCalPoints(sff.clb.calPoints,1000,sff.prm.mov.frameDim);

end


function calPoints = extractCalibrationPoints(xyt1,xyt2,dk)
%EXTRACTCALIBRATIONTRAJECTORIES Extracts calibration trajectories by
% finding synchronized frames with only one flash detected.
% dk returned by dkRobust
%   
% RS, 07/2020

rp = regionprops(xyt1(:,3));
a1 = vertcat(rp.Area);

rp = regionprops(xyt2(:,3)-dk);
a2 = vertcat(rp.Area);

f1 = find(a1==1);
f2 = find(a2==1);

f = intersect(f1,f2);

idx1 = ismember(xyt1(:,3),f);
idx2 = ismember(xyt2(:,3)-dk,f);

calPoints.p1 = xyt1(idx1,:);
calPoints.p2 = xyt2(idx2,:);

end


function calPointsOut = spreadCalPoints(calPointsIn,N,frameDim)
%SPREADCALPOINTS Downsample calibration points
%   

% don't do anything
if N==Inf 
    calPointsOut = calPointsIn;
end

n = length(calPointsIn.p1);

if n < N
    warning('not enough points to downsample')
    calPointsOut = calPointsIn;
else
    q1 = KAZEPoints(calPointsIn.p1(:,1:2));
    u1 = selectUniform(q1,N,frameDim);

    %[~,f] = ismember(u1.Location,calPointsIn.p1(:,1:2),'rows');
    [~,f] = ismember(u1.Location,q1.Location,'rows');

    calPointsOut.p1 = double([u1.Location calPointsIn.p1(f,3)]);
    calPointsOut.p2 = calPointsIn.p2(f,:);

end

end

