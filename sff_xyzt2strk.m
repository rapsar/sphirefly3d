function sff = sff_xyzt2strk(sff)
%SFF_XYZT2STRK Summary of this function goes here
%   Detailed explanation goes here


sff.stk = xyzt2strk(sff.xyzt,sff.prm);
sff.xyztk = sff.stk.xyzti;

end


function strk = xyzt2strk(xyzt,prm)

% build adjacency matrix
strkKernel = determineStrkKernel(xyzt);
adj = buildSparseStrkAdj(xyzt,strkKernel,prm.stk.linkRadiusMtr);

% adjacency matrix to graph
dg = digraph(adj);

% graph connected components, i.e. streaks
[strkID,strkDuration] = conncomp(dg,'type','weak');

%xyzt streakID
xyzti = [xyzt strkID(:)];
xyzti = sortrows(xyzti,[5 4]);

% streaks in cell format
xyzts = mat2cell(xyzti,strkDuration(:));

% streak number of frames
nf = cellfun(@(x) size(x,1), xyzts);

% streak first and last frames
ti = cellfun(@(x) x(1,4), xyzts);
tf = cellfun(@(x) x(end,4), xyzts);

% streak first and last positions
ri = cell2mat(cellfun(@(x) x(1,1:3), xyzts,'UniformOutput',false));
rf = cell2mat(cellfun(@(x) x(end,1:3), xyzts,'UniformOutput',false));

% streak average velocity
v = cellfun(@(x) mean(vecnorm(diff(x(:,1:3),1,1),2,2)), xyzts);

% streak average height
z = cellfun(@(x) mean(x(:,3)), xyzts);

% streak vertical displacement
dz = cellfun(@(x) x(end,3)-x(1,3), xyzts);

% streak horizontal displacement
dxy = cellfun(@(x) vecnorm(x(end,1:2)-x(1,1:2)), xyzts);

% streak distance from reference camera
q = cellfun(@(x) mean(vecnorm(x(:,1:3),2,2)), xyzts);

% out
strk.xyzti = xyzti;
strk.xyzts = xyzts;
strk.nf = nf;
strk.v = v;
strk.ti = ti;
strk.tf = tf;
strk.ri = ri;
strk.rf = rf;
strk.nStreaks = max(strkID);
strk.z = z;
strk.dz = dz;
strk.dxy = dxy;
strk.q = q;

end


function adj = buildSparseStrkAdj(xyzt,strkKernel,strkLinkRadius)
%% number of adjacent streaks to probe for match
x = xyzt(:,1);
y = xyzt(:,2);
z = xyzt(:,3);
t = xyzt(:,4);
p = length(t);

%% build sparse matrix
sp = spdiags(ones(p,2*strkKernel+1),-strkKernel:strkKernel,p,p);
[row,col] = find(sp);

%% flash delays
dt = sparse(row, col, abs(t(row)-t(col)));

%% flash distances
dx = sparse(row, col, x(row)-x(col)+eps); % +eps necessary to avoid numeric zero equated to sparse zero 
dy = sparse(row, col, y(row)-y(col)+eps);
dz = sparse(row, col, z(row)-z(col)+eps);
dr = sqrt(dx.^2+dy.^2+dz.^2);


%% distance-based linkage (distance-adjacency matrix)
adjt = (dt == 1);
adjr = (spfun(@(S) S-strkLinkRadius,dr) < 0);
adj = adjt & adjr;

end

function strkKernel = determineStrkKernel(xyzt)

% simply the maximum number of flashes in a frame
t = xyzt(:,4);
n = histcounts(t,0.5:max(t)+0.5);

strkKernel = max(n);

end

