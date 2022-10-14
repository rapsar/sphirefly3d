function sff = sffTrack(sff)
%SFFTRACK track and clean
% input and output: sff structure
%  
% Raphael Sarfati
% raphael.sarfati@aya.yale.edu


%% launch tracking

if ~sff.prm.flag.trk

    disp([datestr(now,31) ' -- Tracking started...'])

    sff = sff_fffabc(sff);

    % flag; return to workspace
    sff.prm.flag.trk = 1;
    assignin('base','sff0',sff)
    disp([datestr(now,31) ' -- Tracking completed.'])
    
end


%% launch cleaning

if ~sff.prm.flag.cln

    disp([datestr(now,31) ' -- Cleaning started...'])

    sff = sff_clean(sff);

    % flag; return to workspace
    sff.prm.flag.cln = 1;
    assignin('base','sff0',sff)
    disp([datestr(now,31) ' -- Cleaning completed.'])

end


end

