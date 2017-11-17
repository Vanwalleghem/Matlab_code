function varargout = repeattest(varargin)
% REPEATTEST MATLAB code for repeattest.fig
%      REPEATTEST, by itself, creates a new REPEATTEST or raises the existing
%      singleton*.
%
%      H = REPEATTEST returns the handle to a new REPEATTEST or the handle to
%      the existing singleton*.
%
%      REPEATTEST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REPEATTEST.M with the given input arguments.
%
%      REPEATTEST('Property','Value',...) creates a new REPEATTEST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before repeattest_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to repeattest_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help repeattest

% Last Modified by GUIDE v2.5 22-Jun-2015 09:46:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @repeattest_OpeningFcn, ...
                   'gui_OutputFcn',  @repeattest_OutputFcn, ...
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


% --- Executes just before repeattest is made visible.
function repeattest_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to repeattest (see VARARGIN)

% Choose default command line output for repeattest
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes repeattest wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = repeattest_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function RepeatNb_Callback(hObject, eventdata, handles)
% hObject    handle to RepeatNb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RepeatNb as text
%        str2double(get(hObject,'String')) returns contents of RepeatNb as a double
user_entry_X = get(RepeatNb, 'string');
if isnan(user_entry_X)
    errordlg('You must enter a numeric value','Bad Input','modal')
    uicontrol(hObject)
else
    user_entry_X=str2num(user_entry_X);
return
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(repeattest);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(repeattest);
