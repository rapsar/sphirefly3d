function sff = sff_strk2traj(sff)
%SFF_STRK2TRAJ links streaks into trajectories
%   
% Raphael Sarfati, 05/2022

sff.trj = strk2traj(sff.stk,sff.prm);

sff.xyztkj = sff.trj.j(:,[1 2 3 4 5 6]);
sff.xyztkj(:,4) = sff.xyztkj(:,4)./sff.prm.mov.frameRate; %time in seconds

end



%% engine

function traj = strk2traj(strk,prm)
%STRK2TRAJ distance-based linkage based on (sparse) adjacency matrix
% strk -- strk structure returned by xyzt2strk
% prm -- parameter structure
%    .trajLinkRadius : max distance between streaks (meters)
%    .trajLinkMinDel : min delay between streaks (frames)
%    .trajLinkMaxDel : max delay between streaks (frames)
%
% Raphael Sarfati
% raphael.sarfati@aya.yale.edu
%
% rev: -05/22

% build adjacency matrix
trajKernel = determineTrajKernel(strk);
adj = buildSparseTrajAdj(strk,trajKernel,prm.trj.linkRadiusMtr,prm.trj.linkMinLagFrm,prm.trj.linkMaxLagFrm);

% adjacency matrix to graph
dg = digraph(adj);
trajID = conncomp(dg,'type','weak');

% graph connected streaks, ie trajectories
nTraj = max(trajID);
traj.xyztsr = cell(nTraj,1);

for r = 1:nTraj
    
    traj.xyztsr{r} = vertcat(strk.xyzts{trajID == r});
   
    % sort by increasing time
    traj.xyztsr{r} = sortrows(traj.xyztsr{r},4);

end

traj.nTraj = nTraj;
traj.trajLength = cellfun(@(x) size(x,1), traj.xyztsr);
traj.nStreaks = cellfun(@(x) size(unique(x(:,5)),1), traj.xyztsr);
traj.j = [vertcat(traj.xyztsr{:}) repelem((1:nTraj)',traj.trajLength(:))];

end



function adj = buildSparseTrajAdj(strk,trajKernel,trajLinkRadius,trajLinkMinLagFrm,trajLinkMaxLagFrm)
% builds sparse adjacency matrix (to avoid memory overload)

%% build sparse matrix
nStreaks = strk.nStreaks; 
sp = spdiags(ones(nStreaks,2*trajKernel+1),-trajKernel:trajKernel,nStreaks,nStreaks);
[row,col] = find(sp);

%% flash delays
dt = sparse(row, col, strk.ti(row)-strk.tf(col));
dt = dt';

%% flash distances
dx = sparse(row, col, strk.ri(row,1)-strk.rf(col,1));
dy = sparse(row, col, strk.ri(row,2)-strk.rf(col,2));
dz = sparse(row, col, strk.ri(row,3)-strk.rf(col,3));
dr = sqrt(dx.^2+dy.^2+dz.^2);
dr = dr';


%% distance-based linkage (distance-adjacency matrix)
adjtm = dt > trajLinkMinLagFrm;
adjtM = spfun(@(S) S-trajLinkMaxLagFrm,dt) < 0;
adjrM = spfun(@(S) S-trajLinkRadius,dr) < 0;
adj = adjtm & adjtM & adjrM;

end



function trajKernel = determineTrajKernel(strk)

trajKernel = 100;

end

