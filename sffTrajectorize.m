function sff = sffTrajectorize(sff)
%SFFTRAJECTORIZE 
%   
% RS, 05/2022

disp([datestr(now,31) ' -- Trajectorizing started...'])

% streaks
sff = sff_xyzt2strk(sff);
sff.prm.flag.stk = 1;

% trajectories
sff = sff_strk2traj(sff);
sff.prm.flag.trj = 1;

disp([datestr(now,31) ' -- Processing completed.'])


end

