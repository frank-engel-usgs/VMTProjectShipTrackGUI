function varargout = VMT_ProjectShipTrackGUI(varargin)
% VMT_PROJECTSHIPTRACKGUI MATLAB code for VMT_ProjectShipTrackGUI.fig
%      VMT_PROJECTSHIPTRACKGUI, by itself, creates a new VMT_PROJECTSHIPTRACKGUI or raises the existing
%      singleton*.
%
%      H = VMT_PROJECTSHIPTRACKGUI returns the handle to a new VMT_PROJECTSHIPTRACKGUI or the handle to
%      the existing singleton*.
%
%      VMT_PROJECTSHIPTRACKGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VMT_PROJECTSHIPTRACKGUI.M with the given input arguments.
%
%      VMT_PROJECTSHIPTRACKGUI('Property','Value',...) creates a new VMT_PROJECTSHIPTRACKGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before VMT_ProjectShipTrackGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to VMT_ProjectShipTrackGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help VMT_ProjectShipTrackGUI

% Last Modified by GUIDE v2.5 28-Jan-2016 14:22:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @VMT_ProjectShipTrackGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @VMT_ProjectShipTrackGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before VMT_ProjectShipTrackGUI is made visible.
function VMT_ProjectShipTrackGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to VMT_ProjectShipTrackGUI (see VARARGIN)

% Choose default command line output for VMT_ProjectShipTrackGUI
handles.output = hObject;

% Ensure path to utils & docs is available
% ----------------------------------------
if ~isdeployed
    filesep = '\'; % windows
    utilspath = [pwd filesep 'utils'];
    docspath  = [pwd filesep 'doc'];
    toolspath = [pwd filesep 'tools'];
    addpath(utilspath,docspath,toolspath)
end

% Update handles structure
guidata(hObject, handles);

% Create GUI Parameters
guiparams.horizontal_smoothing_window   = 1;
guiparams.vertical_smoothing_window     = 1;
guiparams.water_surface_elevation       = 0;
guiparams.set_cross_section_endpoints   = false;
guiparams.unit_discharge_correction     = false;
guiparams.mcs_id                        = cell(500,1);
guiparams.full_path_to_ascii_file       = [];
guiparams.horizontal_grid_node_spacing  = 1;
guiparams.vertical_grid_node_spacing    = 0.4;
guiparams.data_folder                   = pwd;
guiparams.data_files                    = [];
guiparams.table_data                    = [];
guiparams.shiptracks                    = [];
guiparams.station                       = [];
guiparams.endpoints                     = [];
guiparams.offsets                       = [];
guiparams.left_latitude_dms             = [];
guiparams.left_longitude_dms            = [];
guiparams.right_latitude_dms            = [];
guiparams.right_latitude_dms            = [];
guiparams.left_latitude_dd              = [];
guiparams.left_longitude_dd             = [];
guiparams.right_latitude_dd             = [];
guiparams.right_latitude_dd             = [];

% Store the appication data
setappdata(handles.figure1,'guiparams',guiparams)
% UIWAIT makes VMT_ProjectShipTrackGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = VMT_ProjectShipTrackGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in ProcessTransect.
function ProcessTransect_Callback(hObject, eventdata, handles)
% hObject    handle to ProcessTransect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get the Application data:
% -------------------------
guiparams = getappdata(handles.figure1,'guiparams');
mcs_id = guiparams.mcs_id;
horizontal_grid_node_spacing = guiparams.horizontal_grid_node_spacing;
vertical_grid_node_spacing = guiparams.vertical_grid_node_spacing;
data_folder = guiparams.data_folder;
data_files = guiparams.data_files;
fullfiles = {};
for i = 1:numel(data_files)
    fullfiles(i) = {fullfile(data_folder,data_files{i})};
end
table_data = guiparams.table_data;
Map = [];
Lat = [guiparams.left_latitude_dd guiparams.right_latitude_dd];
Lon = [guiparams.left_longitude_dd guiparams.right_longitude_dd];
%zt      = unique(mcs_id);

% Massage inputs for processing function call
Station = table_data(:,4)';
Station(strcmpi('L',Station)) = {0};
Station(strcmpi('R',Station)) = {1};
Station = cell2mat(Station);

[x,y,utmzone] = deg2utm(Lat,Lon);
Endpoints = [x y];
UTMzone   = utmzone(1,:);
Offsets   = cell2mat(table_data(:,2:3));

A = VMT_ProjectShipTracks(fullfiles,Endpoints,Station,Offsets,UTMzone);

% 4. Add the following variables to the workspace:
  z               = length(A); % number of ASCII files in the cross section
  setends         = 0;         % turn off manual enpoint control 
  unitQcorrection = 0;         % turn off Hoitink correction
  A(1).hgns       = 1;         % Horizontal Grid node spacing in meters
  A(1).vgns       = 0.4;       % Vertical Grid node spacing in meters
  A(1).wse        = 0;         % water surface elev in meters
  Map             = [];        % this would indicate if there was a saved
                               % Map file loaded. Set to empty matrix. 
 
 % Adjust vertical grid size by bin size
 if A(1).Sup.wm ~= 3 % RG
%      set(handles.VerticalGridNodeSpacing,'String',double(A(1).Sup.binSize_cm(1))/100)
     guiparams.vertical_grid_node_spacing = double(A(1).Sup.binSize_cm(1))/100;
 else % Older file, must be RR or M9
%      set(handles.VerticalGridNodeSpacing,'String',0.4)
     guiparams.vertical_grid_node_spacing = 0.4;
 end
 
%%
% Run the PreProcessing Engine
  A = VMT_PreProcess(z,A);

% Run the Processing Engine
start_bank = 'auto';
  [A,V,log_text] = VMT_ProcessTransects(z,A,setends,unitQcorrection,start_bank);
 
% Save the workspace as a MAT file. This MAT file can be processed
% normally by VMT v4. 
  [FileName,PathName,~] = uiputfile(...
       '*.mat','Save processed VMT v4 compliant file as',pwd);
  save(fullfile(PathName,FileName),'A','V','z','Map')




% --- Executes on button press in SelectFiles.
function SelectFiles_Callback(hObject, eventdata, handles)
% hObject    handle to SelectFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get the Application data:
% -------------------------
guiparams = getappdata(handles.figure1,'guiparams');
% guiprefs = getappdata(handles.figure1,'guiprefs');

% Ask the user to select files:
% -----------------------------
current_file = pwd; %fullfile(guiprefs.ascii_path,guiprefs.ascii_file{1});
[filename,pathname] = uigetfile({'*_ASC.TXT','ASCII (*_ASC.TXT)'}, ...
    'Select the ASCII Output Files', ...
    current_file, ...
    'MultiSelect','on');

if ischar(pathname) % The user did not hit "Cancel"
    guiparams.data_folder = pathname;
    if ischar(filename)
        filename = {filename};
    end
    guiparams.data_files = filename;
%     guiparams.mat_file = '';
    
    setappdata(handles.figure1,'guiparams',guiparams)
    
% Populate the table
% Ensure UItable is empty before filling it with current selection
ClearTable_Callback(hObject, eventdata, handles)

% set(handles.FileList,'data',single.empty(500,2,0));

% Construct table
numXS = length(filename);
emptytable = {[] [] []};
table_data = [filename' repmat(emptytable,numXS,1)];
guiparams.mcs_id = num2cell(ones(numel(filename),1));

% Push it to the UItable
set(handles.FileList,'data',table_data);

% Store parameters
guiparams.table_data = table_data;
setappdata(handles.figure1,'guiparams',guiparams)


%     % Update the preferences:
%     % -----------------------
%     guiprefs = getappdata(handles.figure1,'guiprefs');
%     guiprefs.ascii_path = pathname;
%     guiprefs.ascii_file = filename;
%     setappdata(handles.figure1,'guiprefs',guiprefs)
%     store_prefs(handles.figure1,'ascii')
    
end




% --- Executes on button press in LoadTable.
function LoadTable_Callback(hObject, eventdata, handles)
% hObject    handle to LoadTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get application data
guiparams = getappdata(handles.figure1,'guiparams');

[filename,pathname] = uigetfile('*.xlsx','Load Batch File',guiparams.data_folder);
[ndata, text, alldata] = xlsread(fullfile(pathname,filename));
tabledata = alldata(8:end,1:4);
left_latitude_dms   = alldata{4,2};
left_longitude_dms  = alldata{5,2};
right_latitude_dms  = alldata{4,3};
right_longitude_dms = alldata{5,3};

set(handles.LeftLatitude,   'String', left_latitude_dms)
set(handles.LeftLongitude,  'String', left_longitude_dms)
set(handles.RightLatitude,  'String', right_latitude_dms)
set(handles.RightLongitude, 'String', right_longitude_dms)
set(handles.FileList,'data',tabledata);


guiparams.data_folder         = alldata{1,2};
guiparams.data_files          = tabledata(:,1);
guiparams.table_data          = get(handles.FileList,'data');
guiparams.left_latitude_dms   = left_latitude_dms;
guiparams.left_longitude_dms  = left_longitude_dms;
guiparams.right_latitude_dms  = right_latitude_dms;
guiparams.right_longitude_dms = right_longitude_dms;
guiparams.left_latitude_dd    = dms2dd(cellfun(@str2num,strsplit(left_latitude_dms)));
guiparams.left_longitude_dd   = dms2dd(cellfun(@str2num,strsplit(left_longitude_dms)));
guiparams.right_latitude_dd   = dms2dd(cellfun(@str2num,strsplit(right_latitude_dms)));
guiparams.right_longitude_dd  = dms2dd(cellfun(@str2num,strsplit(right_longitude_dms)));


setappdata(handles.figure1,'guiparams',guiparams)



% --- Executes on button press in SaveTable.
function SaveTable_Callback(hObject, eventdata, handles)
% hObject    handle to SaveTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get application data
guiparams = getappdata(handles.figure1,'guiparams');

% Get whatever is in the UItable
data = get(handles.FileList,'data');
numXS = size(data,1);

% Write to an Excel File
[filename,pathname] = uiputfile('*.xlsx','Save Batch File As',guiparams.data_folder);
data = [...
    {'Data Path:'         guiparams.data_folder       []                             []};...
    cell(1,4);...
    {[]                   'Left Endpoint'             'Right Endpoint'               []};...
    {'Latitude:'          guiparams.left_latitude_dms  guiparams.right_latitude_dms  '(+/-dd mm ss.sss)'};...
    {'Longitude:'         guiparams.left_longitude_dms guiparams.right_longitude_dms '(+/-ddd mm ss.sss)'};...
    cell(1,4);...
    {'ASCII File Name(s)' 'Left Offset (m)'            'Right Offset (m)'            'Start Station (L/R)'};...
    data];
xlswrite(fullfile(pathname,filename),data);




% --- Executes on button press in ClearTable.
function ClearTable_Callback(hObject, eventdata, handles)
% hObject    handle to ClearTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Ensure UItable is empty before filling it with current selection
nrows      = 4;
set(handles.FileList,'data',single.empty(nrows,4,0));


function LeftLatitude_Callback(hObject, eventdata, handles)
% hObject    handle to LeftLatitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LeftLatitude as text
%        str2double(get(hObject,'String')) returns contents of LeftLatitude as a double
%  Get Application Data
guiparams = getappdata(handles.figure1,'guiparams');
left_latitude_dms = get(handles.LeftLatitude,'String');

% Check for valid input format
if ValidateCoordinates(left_latitude_dms)
    C = cellfun(@str2num,strsplit(left_latitude_dms));
    guiparams.left_latitude_dms = left_latitude_dms;
    guiparams.left_latitude_dd  = dms2dd(C);
else
    % Reject and replace with default
    set(handles.LeftLatitude,'String',guiparams.left_latitude_dms)
end

  
% Store result
setappdata(handles.figure1,'guiparams',guiparams)


function LeftLongitude_Callback(hObject, eventdata, handles)
% hObject    handle to LeftLongitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LeftLongitude as text
%        str2double(get(hObject,'String')) returns contents of LeftLongitude as a double
%  Get Application Data
guiparams = getappdata(handles.figure1,'guiparams');
left_longitude_dms = get(handles.LeftLongitude,'String');

% Check for valid input format
if ValidateCoordinates(left_longitude_dms)
    C = cellfun(@str2num,strsplit(left_longitude_dms));
    guiparams.left_longitude_dms = left_longitude_dms;
    guiparams.left_longitude_dd  = dms2dd(C);
else
    % Reject and replace with default
    set(handles.LeftLongitude,'String',guiparams.left_longitude_dms)
end
  
% Store result
setappdata(handles.figure1,'guiparams',guiparams)


function RightLatitude_Callback(hObject, eventdata, handles)
% hObject    handle to RightLatitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RightLatitude as text
%        str2double(get(hObject,'String')) returns contents of RightLatitude as a double
%  Get Application Data
guiparams = getappdata(handles.figure1,'guiparams');
right_latitude_dms = get(handles.RightLatitude,'String');

% Check for valid input format
if ValidateCoordinates(right_latitude_dms)
    C = cellfun(@str2num,strsplit(right_latitude_dms));
    guiparams.right_latitude_dms = right_latitude_dms;
    guiparams.right_latitude_dd  = dms2dd(C);
else
    % Reject and replace with default
    set(handles.RightLatitude,'String',guiparams.right_latitude_dms)
end
  
% Store result
setappdata(handles.figure1,'guiparams',guiparams)


function RightLongitude_Callback(hObject, eventdata, handles)
% hObject    handle to RightLongitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RightLongitude as text
%        str2double(get(hObject,'String')) returns contents of RightLongitude as a double
%  Get Application Data
guiparams = getappdata(handles.figure1,'guiparams');
right_longitude_dms = get(handles.RightLongitude,'String');

% Check for valid input format
if ValidateCoordinates(right_longitude_dms)
    C = cellfun(@str2num,strsplit(right_longitude_dms));
    guiparams.right_longitude_dms = right_longitude_dms;
    guiparams.right_longitude_dd  = dms2dd(C);
else
    % Reject and replace with default
    set(handles.RightLongitude,'String',guiparams.right_longitude_dms)
end
  
% Store result
setappdata(handles.figure1,'guiparams',guiparams)


% --- Executes on button press in ClearCoordinates.
function ClearCoordinates_Callback(hObject, eventdata, handles)
% hObject    handle to ClearCoordinates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.LeftLatitude,   'String', '')
set(handles.LeftLongitude,  'String', '')
set(handles.RightLatitude,  'String', '')
set(handles.RightLongitude, 'String', '')

function valid = ValidateCoordinates(coordinates)

% Check for valid input format
if any(isspace(coordinates))
    C = cellfun(@str2num,strsplit(coordinates));
    if numel(C) ~= 3
        % Reject and replace with default
        valid = false;
    end
    valid = true;
else
    % Reject and replace with default
    valid = false;
end


% --- Executes when entered data in editable cell(s) in FileList.
function FileList_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to FileList (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
%  Get Application Data
guiparams = getappdata(handles.figure1,'guiparams');
table_data = hObject.Data;

changed = eventdata.Indices;
if changed(2) == 2 || changed(2) == 3  % Left/Right Offset
    if isnan(eventdata.NewData)
        table_data{changed(1), changed(2)} = eventdata.PreviousData;
        hObject.Data = table_data; % reset back to original value
        warndlg('Offset distances must be a number, and in meters.')
    end
elseif changed(2) == 4 % Start Station
    if ~strcmpi(eventdata.NewData,'l') && ...
            ~strcmpi(eventdata.NewData,'r')
        table_data{changed(1), changed(2)} = eventdata.PreviousData;
        hObject.Data = table_data; % reset back to original value
        warndlg('Start Station must be ''L'' or ''R'' (case insensitive).')
    else
        table_data{changed(1), changed(2)} = upper(eventdata.NewData);
        hObject.Data = table_data;
    end
end

guiparams.table_data = table_data;

% Store result
setappdata(handles.figure1,'guiparams',guiparams)

        




