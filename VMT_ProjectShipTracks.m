function [A,ResultTable] = VMT_ProjectShipTracks(A,Endpoints,Station,Offsets,UTMzone)
%THIS function attempts to build the VMT A struct from input WRII ASC files
%that contain no GPS data. 
%
% To locate the position of each measured ensemble within a transect, a
% method based on ADCP bottom track data that accounts for positional
% errors has been developed. The positions of the ADCP probe at the
% beginning and end of each pass is determined relative to the surveyed
% endpoints of each cross section, accounting for the offset between the
% endpoint and the probe. The position of each ensemble along the transect
% is then computed by distributing the probe’s reported values of linear
% distance traveled (dead reckoning), or “distance made good” (DMG),
% between the beginning and ending locations of the pass. In high flows, an
% actively mobile bed may exist, producing error in the ADCP bottom track
% positions.  In these cases, length scaling adjustments are made to
% ensemble locations by distributing the difference between the measured
% distance and probe-reported distance equally among all ensembles.
% 
% USAGE
% The user must supply the ASCII files, ENDPOINTS in UTM, Starting STATION
% (i.e., channel left or right bank) of each pass, and the  OFFSETS
% distance in meters between the ENDPOINTS and the actual location of the
% ADCP. If the user gives empty inputs "[]", the function will prompt for
% MAT files containing data.
% 
% Data formats are as follows: 
%   Station = [0 1 0 1];                        0 = Channel left, 
%                                                   1 = Channel right
%   Endpoints = [314381.119 4467822.594;...     [X1 Y1; X2 Y2] in UTM (m)
%       314353.275 4467819.957]; 
%   Offsets = [1.4 4.5];                        [LB RB] in meters
%   UTMzone = '16 T';                           UTM Zone code, case
%                                                   sensitive
% 
% A = VMT_ProjectShipTracks([],[],[],[]) Causes function to prompts for all
% inputs and process (easiest). Select all ASC files associated with a
% cross section, and a MAT file containing the cross section STATION,
% ENDPOINTS, OFFSETS, and UTMZONE variables for that cross section.
% 
% If sucessfull, the function should produce a plot (Figure 1) of the cross
% section line, bank stationing, and projected ensemble locations. Green
% ensemble locations are good, whereas red are repeat locations where the
% function has not included them in the resulting A struct. VMT has
% algorithms to replace these missing ensemble locations with interpolated
% valid bottom track.
% 
% Created by: Frank L. Engel (fengel@usgs.gov)
% Last modified: 2014/02/03 (FLE, added comments, descriptions for Ricardo)

% Set to 1 to save A structures automatically
Save = 0;

% Check to see if A was given, prompt user to load ASC files if needed.
% Read in the ASCII data, and make an A struct.
if isempty(A)
    zPathName = uigetdir(pwd,'Select folder containing WR ASCII data');
    Files = dir(zPathName);
    allFiles = {Files.name};
    filefind=strfind(allFiles,'ASC.TXT')';
    filesidx=nan(size(filefind,1),1);
    for i=1:size(filefind,1)
        filesidx(i,1)=size(filefind{i},1);
    end
    filesidx=find(filesidx>0);
    files=allFiles(filesidx);
    
    % Allow user to select which files are to be processed
    selection = listdlg('ListSize',[500 300],'ListString', files,...
        'Name','Select ASC files associated with ONE cross section');
    zFileName = files(selection);
    
    % Determine number of files to be processed
    if  isa(zFileName,'cell')
        z=size(zFileName,2);
        zFileName=sort(zFileName);
    else
        z=1;
        zFileName={zFileName};
    end
    
    for zi = 1:z
        fullName = fullfile(zPathName,zFileName{zi});
        [A(zi)] = tfile(fullName,1,0);
        disp(['ASCII File: ' fullName])
    end
else % Arguement 1 is a cell of full paths to ASC files
    files2load = A;
    z = numel(files2load);
    clear A
    for zi = 1:z
        fullName = files2load{zi};
        [zPathName{zi},zFileName{zi},zExt{zi}] = fileparts(fullName); 
        [A(zi)] = tfile(fullName,1,0);
        disp(['ASCII File: ' fullName])
    end
    
end

% Check to see if the XS Input data were supplied
if isempty(Endpoints) || isempty(Station) || isempty(Offsets)
    [INfile,INpath,~] = uigetfile(...
        '*.mat','Select cross section input-MAT file',zPathName{1});
    load([INpath filesep INfile])
    disp(['Offset File: ' fullfile(INpath,INfile)])
end

% Loop through each transect and process
hf = figure('name','ProjectShiptracks'); clf
set(gca,'DataAspectRatio',[1 1 1])
for zi = 1:z
    disp(['Transect pass ' num2str(zi)])
    
    % If there is GPS data, go ahead and plot it
    isGPS = any(~isnan(A(zi).Nav.lat_deg));
    if isGPS
        [xUTMraw, yUTMraw, utmzone] = deg2utm(...
            A(zi).Nav.lat_deg,A(zi).Nav.long_deg);
        plot(xUTMraw,yUTMraw,'b'); hold on
    end
    
    % Make a string of station for clarity
    if Station(zi) == 1
        sta{zi} = 'right';
    else
        sta{zi} = 'left';
    end
    
    % Compute dx, dy, dl, & m of the line, always going from LB to RB
    x1 = Endpoints(1,1); x2 = Endpoints(2,1);
    y1 = Endpoints(1,2); y2 = Endpoints(2,2);
    dl(zi) = hypot((x2-x1),(y2-y1));
    m = atan((y2-y1)/(x2-x1));
    dmg{zi} = nanmax(A(zi).Nav.dmg);
    
%     disp(['   Endpoint distance:        ' num2str(dl(zi))]);
        
    % Plot the line of the mean cross section
    figure(hf);hold on;plot([x1 x2],[y1 y2],'-'); axis equal
    
    % Overwrite the original endpoints, moving them to where the ADCP probe
    % "started"
    if x1>x2 % left bank is east of right bank
        if strcmpi(sta{zi},'left')
            % Start from LB, subtract offset, find point on line where ADCP
            % "started"
            [x1, y1] = line_projection(x1,y1,-Offsets(zi,1),0,m);
            [x2, y2] = line_projection(x2,y2,Offsets(zi,2),0,m);
        elseif strcmpi(sta{zi},'right')
            % Start from RB, add offset, find point on line where ADCP
            % "started"
            [x1, y1] = line_projection(x1,y1,-Offsets(zi,1),0,m);
            [x2, y2] = line_projection(x2,y2,Offsets(zi,2),0,m);
        end
        dl2(zi) = hypot((x2-x1),(y2-y1));
        
        % Compute positional error and display results
        abs_err(zi) = abs(dl2(zi) - dmg{zi});
        per_err(zi) = (abs_err(zi)/abs(dl2(zi))) * 100;
        probe_length(zi) = nanmax(A(zi).Nav.length);
%         disp(['   Offset distance is:       ' num2str(dl2(zi))]);
%         disp(['   Probe reported DMG:       ' num2str(dmg{zi})]);
%         disp(['   Probe max length:         ' num2str(probe_length(zi))]);
%         disp(['   Absolute error:           ' num2str(abs_err(zi))]);
%         disp(['   Percent error:            ' num2str(per_err(zi))]);
        
        % Plot where the ADCP started on the MCS, accounting for offsets
        figure(hf); hold on; plot([x1 x2],[y1 y2],'o');
        text(x1,y1,'LB'); text(x2,y2,'RB'); hold off
        
        % Scale the DMG to span from the ADCP starting points. This
        % distributes the locational error equally along the transect
        for j = 1:A(zi).Sup.noe
            origLength(j) = A(zi).Nav.length(j);
            Ratio(j)      = origLength(j)/nanmax(A(zi).Nav.length);
            adjLength(j)  = dl2(zi)*Ratio(j);
        end
        
        % Traverse the ADCP probe along the cross section using the
        % adjusted length
        if strcmpi(sta{zi},'left')
            if m > 0
                % Start from LB, traverse along slope, subtracting station
                % as you go
                [xUTM, yUTM] = line_projection(...
                    x1,y1,...
                    -adjLength,zeros(size(adjLength)),...
                    m);
            else
                % Start from LB, traverse along slope, adding station as
                % you go
                [xUTM, yUTM] = line_projection(...
                    x1,y1,...
                    -adjLength,zeros(size(adjLength)),...
                    m);
            end
            
        elseif strcmpi(sta{zi},'right')
            if m>0
                % Start from RB, traverse along slope, adding station as
                % you go
                [xUTM, yUTM] = line_projection(...
                    x2,y2,...
                    adjLength,zeros(size(adjLength)),...
                    m);
            else
                % Start from RB, traverse along slope, subtracting station
                % as you go
                [xUTM, yUTM] = line_projection(...
                    x2,y2,...
                    adjLength,zeros(size(adjLength)),...
                    m);
            end
        end
     
    elseif x1<x2 % left bank west of right bank
        if strcmpi(sta{zi},'left')
            % Start from LB, subtract offset, find point on line where ADCP
            % "started"
            [x1, y1] = line_projection(x1,y1,Offsets(1),0,m);
            [x2, y2] = line_projection(x2,y2,-Offsets(2),0,m);
        elseif strcmpi(sta{zi},'right')
            % Start from RB, add offset, find point on line where ADCP
            % "started"
            [x1, y1] = line_projection(x1,y1,Offsets(1),0,m);
            [x2, y2] = line_projection(x2,y2,-Offsets(2),0,m);
        end
        dl2(zi) = hypot((x2-x1),(y2-y1));
        
        % Compute positional error and display results
        abs_err(zi) = abs(dl2(zi) - dmg{zi});
        per_err(zi) = (abs_err(zi)/abs(dl2(zi))) * 100;
        probe_length(zi) = nanmax(A(zi).Nav.length);
%         disp(['   Offset distance is:       ' num2str(dl2(zi))]);
%         disp(['   Probe reported DMG:       ' num2str(dmg{zi})]);
%         disp(['   Probe max length:         ' num2str(probe_length(zi))]);
%         disp(['   Absolute error:           ' num2str(abs_err(zi))]);
%         disp(['   Percent error:            ' num2str(per_err(zi))]);
        
        % Plot where the ADCP started on the MCS, accounting for offsets
        figure(hf); hold on; plot([x1 x2],[y1 y2],'o');
        text(x1,y1,'LB'); text(x2,y2,'RB'); hold off
        
        % Scale the DMG to span from the ADCP starting points. This
        % distributes the locational error equally along the transect
        for j = 1:A(zi).Sup.noe
            origLength(j) = A(zi).Nav.length(j);
            Ratio(j)      = origLength(j)/nanmax(A(zi).Nav.length);
            adjLength(j)  = dl2*Ratio(j);
        end
        
        % Traverse the ADCP probe along the cross section using the
        % adjusted length
        if strcmpi(sta{zi},'left')
            if m > 0
                % Start from LB, traverse along slope, subtracting station
                % as you go
                [xUTM, yUTM] = line_projection(...
                    x1,y1,...
                    adjLength,zeros(size(adjLength)),...
                    m);
            else
                % Start from LB, traverse along slope, adding station as
                % you go
                [xUTM, yUTM] = line_projection(...
                    x1,y1,...
                    adjLength,zeros(size(adjLength)),...
                    m);
            end
            
        elseif strcmpi(sta{zi},'right')
            if m>0
                % Start from RB, traverse along slope, adding station as
                % you go
                [xUTM, yUTM] = line_projection(...
                    x2,y2,...
                    -adjLength,zeros(size(adjLength)),...
                    m);
            else
                % Start from RB, traverse along slope, subtracting station
                % as you go
                [xUTM, yUTM] = line_projection(...
                    x2,y2,...
                    -adjLength,zeros(size(adjLength)),...
                    m);
            end
        end
    end
    
    
    % Transpose the coordinates to make them have the correct dimensions
    % for what VMT expects
    xUTM = xUTM';
    yUTM = yUTM';
    
    % VMT screens for bad GPS, and in normal situations it will replace bad
    % GPS with valid bottom track data. To make this functionality work
    % when there are dropped ensembles, we have to "fake" a RG bad data tag
    % (-32768) in the totDist variables. To do this, find where projected
    % x,y UTM coodinates repeat, and take only the first repeat ensemble.
    % The rest are censored.
    [C,ia,ic] = unique(xUTM,'first');
    tmp       = nan(size(xUTM)); 
    tmp(ia)   = C;
    bvbt      = isnan(tmp);
    A(zi).Nav.totDistEast(bvbt)  = -32768;
    A(zi).Nav.totDistNorth(bvbt) = -32768;
    
    % Plot the bottom track ensemble locations
    figure(1); hold on; plot(xUTM,yUTM,'r.'); hold on
    plot(xUTM(ic),yUTM(ic),'g.')
    
    % Write the result back to A struct as lat/lon and UTM
    A(zi).Nav.xUTM = nan(size(xUTM));
    A(zi).Nav.yUTM = nan(size(yUTM));
    A(zi).Nav.xUTM(ic) = xUTM(ic);
    A(zi).Nav.yUTM(ic) = yUTM(ic);
    utmz = repmat(UTMzone,size(xUTM));
    [...
        A(zi).Nav.lat_deg,...
        A(zi).Nav.long_deg] = utm2deg(xUTM',yUTM',utmz);
    clear utmz
    
    
    % Ensure that the screening of bad data worked
    badidx = -32768;
    A(zi).Nav.depth(A(zi).Nav.depth==badidx)     = nan;
    A(zi).Nav.bvEast(A(zi).Nav.bvEast==badidx)   = nan;
    A(zi).Nav.bvError(A(zi).Nav.bvError==badidx) = nan;
    A(zi).Nav.bvNorth(A(zi).Nav.bvNorth==badidx) = nan;
    
    
       
    clear j origLength Ratio adjLength
end

% Create table with results
EndpointDistance = dl';
DistOffsetApplied = dl2';
ADCP_DMG = cell2mat(dmg)';
ADCPLength = probe_length';
AbsoluteError = abs_err';
PercentError = per_err';

ResultTable = table(...
    EndpointDistance,...
    DistOffsetApplied,...
    ADCP_DMG,...
    ADCPLength,...
    AbsoluteError,...
    PercentError,...
    'RowNames',zFileName);
ResultTable.Properties.VariableUnits = {'m' 'm' 'm' 'm' 'm' 'percent'};
ResultTable

if Save
    [OUTfile,OUTpath,~] = uiputfile(...
        '*.mat','Save processed VMT A structure as',zPathName);
    save(fullfile(OUTpath,OUTfile), 'A')
end


%%%%%%%%%%%%%%%%
% SUBFUNCTIONS %
%%%%%%%%%%%%%%%%
function [rx, ry] = line_projection(X1,Y1,X2,Y2,phi)
% Project a line onto a given coordinate system

% Rotation matrix
Rz = [cos(phi) -sin(phi);...
    sin(phi) cos(phi)];

% Rotate every element in the matrix
for i=1:size(X2,2)
    for j = 1:size(X2,1)
        XY = [X2(j,i);Y2(j,i)];
        Rotated = Rz*XY;
        Rotated(1) = Rotated(1) + X1;
        Rotated(2) = Rotated(2) + Y1;
        
        rx(j,i) = Rotated(1);
        ry(j,i) = Rotated(2);
        
    end
end