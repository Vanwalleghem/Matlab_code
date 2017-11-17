function varargout = NbRepeats(varargin)
% NBREPEAT MATLAB code for NbRepeat.fig
%      NBREPEAT, by itself, creates a new NBREPEAT or raises the existing
%      singleton*.
%
%      H = NBREPEAT returns the handle to a new NBREPEAT or the handle to
%      the existing singleton*.
%
%      NBREPEAT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NBREPEAT.M with the given input arguments.
%
%      NBREPEAT('Property','Value',...) creates a new NBREPEAT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NbRepeats_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NbRepeats_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NbRepeat

% Last Modified by GUIDE v2.5 22-Jun-2015 10:04:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NbRepeats_OpeningFcn, ...
                   'gui_OutputFcn',  @NbRepeats_OutputFcn, ...
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


% --- Executes just before NbRepeat is made visible.
function NbRepeats_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NbRepeat (see VARARGIN)

% Choose default command line output for NbRepeat
handles.output = hObject;


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes NbRepeat wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = NbRepeats_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.output;


% --- Executes on button press in OK.
function OK_Callback(hObject, eventdata, handles)
% hObject    handle to OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
NbRepeat_Callback(hObject, eventdata, handles)
NumbRep=handles.NbRepeat;
close(NbRepeats);



% --- Executes on button press in Cancel.
function Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(NbRepeats);


function NbRepeat_Callback(hObject, eventdata, handles)
% hObject    handle to NbRepeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NbRepeat as text
%        str2double(get(hObject,'String')) returns contents of NbRepeat as a double
num = str2double(get(hObject,'String')); 
if isnan(num)
    errordlg('You must enter a numeric value','Bad Input','modal')
else
handles.NbRepeat = num;
guidata(hObject,handles);
return
end


% --- Executes during object creation, after setting all properties.
function NbRepeat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NbRepeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
