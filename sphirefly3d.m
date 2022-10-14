function sff = sphirefly3d(varargin)
%SPHIREFLY3D calibration-free 3D reconstruction of firefly spherical videos  
%
% Raphael Sarfati
% raphael.sarfati@aya.yale.edu
%
% spherifly3d(sff)
% spherifly3d(cyber2world)
%   
%
% 05/2022 -v0.1
% 07/2022 -v0.9

sff = initsff(varargin);

sff = sffTrack(sff);

sff = sffTriangulate(sff);

sff = sffTrajectorize(sff);

exitsff(sff);

end


%% initialize 
function sff = initsff(varargin)
% initialize sff structure based on input arguments

argin = varargin{1};
narg = length(argin);

if (narg == 1) && isstruct(argin{1}) && isfield(argin{1},'gp1')
    disp('input is sff')
    sff = argin{1};

elseif (narg == 1) && isstruct(argin{1}) && isfield(argin{1},'flag')
    disp('input is prm')

    [v1,v2] = getv1v2();

    sff.prm = argin{1};

    sff.gp1.mov = v1;
    sff.gp2.mov = v2;

elseif (narg == 1) && isnumeric(argin{1})
    disp('input is cyber2world')

    [v1,v2] = getv1v2();
    
    w = VideoReader(v1{1});
    sff.prm = setprm(argin{1},w);

    sff.gp1.mov = v1;
    sff.gp2.mov = v2;

    

elseif (narg == 2)

elseif (narg == 3) && isstruct(argin{3})
    disp('input is v1, v2, prm')

elseif (narg == 3) && isnumeric(argin{3})
    disp('input is v1, v2, cyber2world')

    sff.prm = setprm(argin{3},argin{1});

    sff.gp1.mov = argin{1};
    sff.gp2.mov = argin{2};

else
    error('Invalid argument(s). Please refer to sphirefly help.')

end

end


%% exit
function exitsff(sff)
% exits by plotting 3D points and saving xyzt as .csv files

figure,
scatter3(sff.xyzt(:,1),sff.xyzt(:,2),sff.xyzt(:,3),10,sff.xyzt(:,4),'filled')
axis equal
xlim([-10 10]), xlabel('x (m)')
ylim([-10 10]), ylabel('y (m)')
zlim([-10 10]), zlabel('z (m)')

% save main output as csv file
if ~isfile('xyztkj.csv')
    csvname = 'xyztkj.csv';
else
    % avoids overwriting previous file
    csvname = strcat('xyztkj','_',datestr(now,30),'.csv');
end
writematrix(sff.xyztkj,csvname)

end


%% get paths of movie files
function [v1,v2] = getv1v2()

root = pwd;

[file,path] = uiget(pwd,'MultiSelect',true);

if length(path) == 1
    [file(2),path(2)] = uiget(pwd,'MultiSelect',true);
elseif length(path) > 2
    error('Select only 2 video files/folders.')
end


if ~isempty(file{1})                    %video files selected
    v1{1} = fullfile(path(1),file(1));
    v2{1} = fullfile(path(2),file(2));

else                                    % video folders selected
    cd(path(1))
    vid = dir('*.MP4');
    for i=1:length(vid)
        v1{i} = fullfile(path(1),vid(i).name);
    end

    cd(path(2))
    vid = dir('*.MP4');
    for i=1:length(vid)
        v2{i} = fullfile(path(2),vid(i).name);
    end

end

% return to root folder
cd(root)

end


%% select multiple files or folders
function [file, path] = uiget(basepath, varargin)
% UIGET generic folder and/or file selection dialog box
% From Matlab File Exhange:
% sco1 (2022). uiget (https://github.com/StackOverflowMATLABchat/uiget), GitHub. Retrieved May 25, 2022.
%
% Syntax:
%     file = uiget()
%     [file, path] = uiget()
%     ___ = uiget(basepath)
%     ___ = uiget(basepath, Name, Value)
%
% Available Name, Value Pairs:
%     MultiSelect      - Specify whether a user can select multiple files or folders
%     ScalarPathOutput - Specify whether a scalar path is output when using MultiSelect
%     Title            - Specify a custom dialog title
%     ExtensionFilter  - Specify a custom file extension filter
%     ForceCharOutput  - Force return of a cell array of char
%
% See README.md for detailed documentation and examples
%
% See also UIGETDIR, UIGETFILE
if nargin == 0
    % Use current working directory if no inputs are passed
    basepath = pwd;
elseif nargin == 1
    % Check for existence of basepath as a directory, default to current
    % directory if it doesn't exist
    if ~exist(basepath, 'dir')
        basepath = pwd;
    end
end
% Parse additional inputs
p = buildParser();
p.parse(varargin{:});
% Initialize JFileChooser window
% https://docs.oracle.com/javase/8/docs/api/javax/swing/JFileChooser.html
jFC = javax.swing.JFileChooser(basepath);
jFC.setFileSelectionMode(jFC.FILES_AND_DIRECTORIES)
jFC.setDialogTitle(p.Results.Title)
% Build file filter
if ~isempty(p.Results.ExtensionFilter)
    extensions = parsefilter(p.Results.ExtensionFilter(:, 1));
    
    nfilters = size(p.Results.ExtensionFilter, 1);
    for ii = 1:nfilters
        if isempty(extensions{ii})
            % Catch invalid extension specs
            continue
        end
        jExtensionFilter = javax.swing.filechooser.FileNameExtensionFilter(p.Results.ExtensionFilter{ii, 2}, extensions{ii});
        jFC.addChoosableFileFilter(jExtensionFilter)
    end
    
    tmp = jFC.getChoosableFileFilters();
    jFC.setFileFilter(tmp(2))
end
if p.Results.MultiSelect
    jFC.setMultiSelectionEnabled(true)
    
    % Change title if default is being used
    if any(strcmp(p.UsingDefaults, 'Title'))
        jFC.setDialogTitle('Select File(s) and/or Folder(s)')
    end
else
    jFC.setMultiSelectionEnabled(false)
end
% Switch over possible responses from JFileChooser
returnVal = jFC.showOpenDialog([]);
switch returnVal
    case jFC.APPROVE_OPTION
        % Selection string will be empty if getSelectedFiles is used when
        % MultiSelect is disabled
        if jFC.isMultiSelectionEnabled
            selectionStr = string(jFC.getSelectedFiles());
        else
            selectionStr = string(jFC.getSelectedFile());
        end
    case jFC.CANCEL_OPTION
        % Short-circuit: Return empty array on cancel
        file = "";
        path = "";
        return
    otherwise
        err = MException("uiget:JFileWindow:unsupportedResult", ...
                         "Unsupported result returned from JFileChooser: %s.\n" + ...
                         "Please consult the documentation for the current MATLAB Java version (%s)", ...
                         returnVal, string(java.lang.System.getProperty("java.version")));
        err.throw()
end
npicked = numel(selectionStr);
file = strings(npicked, 1);
path = strings(npicked, 1);
for ii = 1:npicked
    [path(ii), filename, ext] = fileparts(selectionStr(ii));
    file(ii) = filename + ext;
    
    % Because we can select directories, we want to have them output as a
    % path and not a file
    if verLessThan('matlab','9.4')
        % string inputs to fullfile were silently added in R2018a, use char
        % for R2017a and R2017b
        tmppath = char(path(ii));
        tmpfile = char(file(ii));
        if exist(fullfile(tmppath, tmpfile), 'dir')
            path(ii) = string(fullfile(tmppath, tmpfile));
            file(ii) = "";
        end
    else
        if exist(fullfile(path(ii), file(ii)), 'dir')
            path(ii) = fullfile(path(ii), file(ii));
            file(ii) = "";
        end
    end
end
% Since we've now adjusted file in cases where a folder was selected, warn
% the user if file is going to be empty and they're not requesting path to
% go with it
emptyfiletest = (file == '');
if nargout <= 1 && any(emptyfiletest)
    warning("uiget:uiget:nopathoutputrequested", ...
            "One or more paths have been selected without requesting the path output.\n" + ...
            "Please specify a second output to uiget to receive these paths." ...
            );
end
% Simplify path output if needed
if p.Results.ScalarPathOutput
    % Check for number of unique paths
    % If more than one is present, use the first & throw a warning
    uniquepaths = unique(path);
    nuniquepaths = numel(uniquepaths);
    if nuniquepaths == 1
        path = uniquepaths;
    elseif nuniquepaths > 1
        path = uniquepaths(1);
        
        warning("uiget:ScalarPathOutput:multipleuniquepaths", ...
                "Multiple unique paths selected, ignoring %u extra selections.", ...
                nuniquepaths - 1);
    end
end
% Convert to cell array of char if flag is set
if p.Results.ForceCharOutput
    file = cellstr(file);
    path = cellstr(path);
end
end
function p = buildParser()
% Validate input Name,Value pairs
% Initialize verbosely, since inputParser apparently doesn't have a
% constructor that takes inputs...
p = inputParser();
p.FunctionName = 'uiget';
p.CaseSensitive = false;
p.KeepUnmatched = true;
p.PartialMatching = false;
% Add Name,Value pairs
p.addParameter('MultiSelect', false, @(x)islogical(x))
p.addParameter('ScalarPathOutput', false, @(x)islogical(x))
p.addParameter('Title', 'Select File or Folder', @(x)validateattributes(x, {'char', 'string'}, {'scalartext'}))
p.addParameter('ExtensionFilter', [], @(x)validateattributes(x, {'cell'}, {'ncols', 2}))
p.addParameter('ForceCharOutput', false, @(x)islogical(x))
end
function extensions = parsefilter(incell)
% Parse the extension filter extensions into a format usable by 
% javax.swing.filechooser.FileNameExtensionFilter
% 
% Since we're keeping with the uigetdir-esque extension syntax
% (e.g. *.extension), we need strip off '*.' from each for compatibility
% with the Java component.
extensions = cell(size(incell));
for ii = 1:numel(incell)
    exp = '\*\.(\w+)';
    extensions{ii} = string(regexp(incell{ii}, exp, 'tokens'));
end
end



