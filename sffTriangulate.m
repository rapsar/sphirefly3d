function sff = sffTriangulate(sff)
%SFFTRIANGULATE calibrate and triangulate
% input and output: sff structure
%
% Raphael Sarfati
% raphael.sarfati@aya.yale.edu


%% launch calibration

if ~sff.prm.flag.clb

    disp([datestr(now,31) ' -- Calibration started...'])

    % time calibration
    sff = sff_dk(sff);

    % space calibration
    sff = sff_extractCalibrationTrajectories(sff); 
    sff = sff_estimate360CameraParameters(sff); 

    % flag; return to workspace
    sff.prm.flag.clb = 1;
    assignin('base','sff0',sff)
    disp([datestr(now,31) ' -- Calibration completed.'])

end


%% launch triangulation

if ~sff.prm.flag.trg

    disp([datestr(now,31) ' -- Triangulation started...'])

    sff = sff_matchPoints(sff); 
    sff = sff_triangulate360(sff); 

    % flag; return to workspace
    sff.prm.flag.trg = 1;
    assignin('base','sff0',sff)
    disp([datestr(now,31) ' -- Triangulation completed.'])    

end
 

end

