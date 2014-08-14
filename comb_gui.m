function varargout = comb_gui(varargin)
% COMB_GUI MATLAB code for comb_gui.fig
%      COMB_GUI, by itself, creates a new COMB_GUI or raises the existing
%      singleton*.
%
%      H = COMB_GUI returns the handle to a new COMB_GUI or the handle to
%      the existing singleton*.
%
%      COMB_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COMB_GUI.M with the given input arguments.
%
%      COMB_GUI('Property','Value',...) creates a new COMB_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before comb_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to comb_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help comb_gui

% Last Modified by GUIDE v2.5 14-Aug-2014 14:04:21

% Begin initialization code - DO NOT EDIT

clc;
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @comb_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @comb_gui_OutputFcn, ...
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
end

% --- Executes just before comb_gui is made visible.
function comb_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to comb_gui (see VARARGIN)
global timesteps_cme;
global hbar;
global c;
global pi;

global modes_number;
global pump_freq;
global fsr;
global d2;
global d3;
global d4;
global d5;
global nms_a;
global nms_b;
global linewidth;
global reslist;
global pump;
global detuning_profile;
global initial;
global detuning_string;
global initial_string;
global pump_string;
global plotsteps;
global solver;
global cavity_linewidths;
global cavity_linewidths_string;
global noise_to_pump;
global coupling;
global refr_index;
global nonlin_index;
global mode_area;
global sweep_speed;


% constants
hbar = 1.05457148e-34;
c = 299792458;
pi = 3.14159;
timesteps_cme = 2^20;
plotsteps = 2048;

% default parameters
solver='Runge Kutta adaptive';
modes_number=201;
lambda=1553*10^-9; % in m
pump_freq=c/lambda;
fsr=35.2e9;
d2=10e3;
d3=-130;
d4=0;
d5=0;
nms_a=0;
nms_b=0;

coupling=0.5;
refr_index=1.37;
nonlin_index=0.9e-20;
mode_area=5.6e-10;
sweep_speed=0.1;

cavity_linewidths_string='cavity_linewidths=450e3*ones(1,modes_number)';
% cavity_linewidths_string='for kk=1:modes_number ind=kk-round(modes_number/2); cavity_linewidths(kk)=1e6+9e6/(round(modes_number/2))^2*ind^2; end';
%TODO: should be similar to initial profile, otherwise user has to know the
%meaning of variable
eval(cavity_linewidths_string);
reslist = eigenmodes(modes_number, pump_freq, fsr, d2, d3, d4, d5, nms_a, nms_b, cavity_linewidths(round(modes_number/2)));

pump_string='100e-3*ones(1,timesteps_cme)';
% pump_string='[10^-3*linspace(0,50,timesteps_cme/2) 50e-3*ones(1,timesteps_cme/2)]';
pump=eval(pump_string);

detuning_string='linspace(-3,15,timesteps_cme)';
% detuning_string='[linspace(-5,5,timesteps_cme/2) 5*ones(1,timesteps_cme/2)]';
detuning_profile=eval(detuning_string);

initial_string='randn(1,modes_number)+1i*randn(1,modes_number)';
% initial_string='sech(linspace(-pi,pi,modes_number))+1i*zeros(1,modes_number)';
initial=eval(initial_string);
% initial(1:4)=0;
% initial(modes_number-5:modes_number)=0;

noise_to_pump=false;
% Choose default command line output for comb_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% set(handles.inject,'Enable','off');

% This sets up the initial plot - only do when we are invisible
% so window can get raised using comb_gui.
if strcmp(get(hObject,'Visible'),'off')
% constants & parameters
end
switch solver
    case 'ode45'
        solver_value=1;
    case 'Runge Kutta adaptive'
        solver_value=2;
    case 'ode23s'
        solver_value=3;
end
set(handles.solver,'Value',solver_value);
set(handles.coupling,'String',num2str(coupling));
set(handles.refr_index,'String',num2str(refr_index));
set(handles.nonlin_index,'String',num2str(nonlin_index));
set(handles.mode_area,'String',num2str(mode_area));
set(handles.sweep_speed,'String',num2str(sweep_speed));

imshow('epfl_small_logo.png','Parent',handles.axes12);
imshow('rqc_logo.png','Parent',handles.axes14);
axes(handles.axes);
end

% UIWAIT makes comb_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = comb_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

% --- Executes on button press in solve_cme.
function solve_cme_Callback(hObject, eventdata, handles)
% hObject    handle to solve_cme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global modes_number;
global fsr;
global d1;
global linewidth;
global coupling;
global refr_index;
global nonlin_index;
global mode_area;
global filename;
global progress;
global snapshot;
global reslist;
global sweep_speed;
global pump;
global pump_profile;
global initial_conditions;
global initial;
global initial_string;
global detuning_profile;
global timesteps_cme;
global kappa;
global omega;
global omega0;
global hbar; global c; global pi;
global pump_freq; global d2, global d3, global d4, global d5, global nms_a, global nms_b
global injection_detuning;
global injection_force;
global cavity_linewidths;
global kappas;
global norm_factor;

% disable injection in simulation
injection_detuning=0;
injection_force=0;

progress=hObject;
set(hObject,'Enable','off');
set(handles.slider3,'Enable','off');
% set(handles.inject,'Enable','off');
pause(.01);

if ~isempty(get(handles.linewidth,'String'))
    cavity_linewidths=str2num(get(handles.linewidth,'String'))*ones(1,modes_number);
end
if length(cavity_linewidths)~=modes_number
    display('Wrong cavity linewidths array length. Using default cavity linewidth = 1e6 Hz')
    cavity_linewidths=1e6*ones(1,modes_number);
end

filename_apndx=datestr(fix(clock),'yyyymmddHHMMSS');
% save parameters in txt file
names = {'modes_number';'fsr';'pump_frequency';'d2';'d3';'d4';'d5';'distortion_a';'distortion_b';'linewidth';'coupling';'refr_index';'nonlin_index';'mode_area';'sweep_speed';'pump';'detuning_start';'detuning_end'};
values = [modes_number;fsr;pump_freq;d2;d3;d4;d5;nms_a;nms_b;cavity_linewidths(round(modes_number/2));coupling;refr_index;nonlin_index;mode_area;sweep_speed;pump(1);detuning_profile(1);detuning_profile(end)];
tbl = table(values,'RowNames',names);
writetable(tbl,strcat('params',filename_apndx,'.txt'),'WriteRowNames',true);

filename=strcat('coupledeq',filename_apndx);
% snapshot=detuning_profile(1)+(detuning_profile(end)-detuning_profile(1))/2;
snapshot=0.5;
% adjust slider to snapshot
min_slider=get(handles.slider3,'Min');
max_slider=get(handles.slider3,'Max');
slider_value=snapshot*(max_slider-min_slider)+min_slider;

if ~isempty(get(handles.detuning_start,'String')) && ~isempty(get(handles.detuning_end,'String'))
    detuning_profile=linspace(str2num(get(handles.detuning_start,'String')),str2num(get(handles.detuning_end,'String')),timesteps_cme);
end

omega = 2*pi*reslist; % resonance frequencies
d1 = 2*pi*fsr;
% TODO: should be in agreement with import eigenmodes
kappa = cavity_linewidths(round(modes_number/2))*2*pi;
kappas=cavity_linewidths*2*pi;
omega0 = omega(round(modes_number/2)); % central pumped frequency
eta = coupling; % coupling coefficient
n0 = refr_index; % refractive index
n2 = nonlin_index; % nonlinear refractive index
Veff=mode_area*c/n0/fsr;
g = hbar*omega0^2*c*n2/n0^2/Veff; % nonlinear coupling coefficient
pump_profile=sqrt(8*eta*g/kappa^2*pump/hbar/omega0); % normalized amplitutde of input field
norm_factor=sqrt(2*g/kappa)*0.5;

if length(initial)==modes_number
    initial_conditions=sqrt(2*g/kappa)*0.5*initial; % normalized initial conditions
else
    display('Wrong initial conditions length. Using default')
    initial=randn(1,modes_number)+1i*randn(1,modes_number);
    initial_conditions=sqrt(2*g/kappa)*0.5*initial;
    initial_string=mat2str(initial_conditions);    
end

simulate()
set(hObject,'String','Solve CMEs');
set(hObject,'Enable','on');
plotcomb(filename,snapshot,'all');
set(handles.slider3,'Enable','on');
set(handles.slider3,'Value',slider_value);
% set(handles.inject,'Enable','on');
end

function res = coupled_modes_equations(t,a)

global kappa;
global kappas;
global omega;
global omega0;
global d1;
global start_time;
global end_time;
global done;	
global modes_number;
global detuning_profile;
global pump_profile;
global injection_detuning;
global injection_force;
global injection_progress;
global noise_to_pump;
global norm_factor;

% TODO: injection_detuning might be 0
if injection_detuning == 0
    detuning_indx=round(length(detuning_profile)*(t-start_time)/(end_time-start_time));
    if detuning_indx == 0 
        detuning = detuning_profile(1);
    else 
        detuning = detuning_profile(detuning_indx);
    end
else
    detuning=injection_detuning;
end
if injection_force == 0
    force_indx=round(length(pump_profile)*(t-start_time)/(end_time-start_time));
    if force_indx == 0
        force=pump_profile(1);
    else
        force=pump_profile(force_indx);
    end
else
    force=injection_force;
end
% nonlinear term
fa = fft(a);
fNL = fa.*conj(fa).*fa;
NL = ifft(fNL);
res = zeros(modes_number,1);

for k = 1:modes_number
%   res(k)=-(1+1i*2/kappa*double((omega(k)-omega0+(detuning*kappa)-(k-round(modes_number/2))*d1)))*a(k)+1i*NL(k);
%   res(k)=-(1+1i*2/kappas(k)*double((omega(k)-omega0+(detuning*kappa)-(k-round(modes_number/2))*d1)))*a(k)+1i*NL(k)+norm_factor*(randn(1)+1i*randn(1));
    res(k)=-(1+1i*2/kappas(k)*double((omega(k)-omega0+(detuning*kappa)-(k-round(modes_number/2))*d1)))*a(k)+1i*NL(k);
end
% adding pump to the central mode (pumped mode)
if noise_to_pump == true
    res(round(modes_number/2)) = res(round(modes_number/2)) + force + 0.001*randn*force;
else
    res(round(modes_number/2)) = res(round(modes_number/2)) + force;
end

% estimate elapsed progress
done_tmp = 100*(t-start_time)/(end_time-start_time);
global progress;
if done_tmp - done > 1
    done = 100*(t-start_time)/(end_time-start_time);
    if injection_force==0
        set(progress,'String',[num2str(round(done)+1) ' %']);
    else
        set(injection_progress,'String',[num2str(round(done)+1) ' %']);
    end
    pause(.001);
end
end

function simulate()

global done;
global detuning_profile;
global sweep_speed;
global timesteps_cme;
global initial_conditions;
global modes_number;
global filename;
global start_time;
global end_time;
global plotsteps;
global injection_force;
global injection_filename;
global solver;

tic
start_time = 0;
% TODO: not correct for nonlinear detuning
if injection_force ==0
    end_time = (detuning_profile(end)-detuning_profile(1))/sweep_speed+start_time;
else
    % we use hardcoded value for time in soliton injection simulation
    end_time=50;
end
done = 0;
timepoints = linspace(start_time,end_time,plotsteps);
options = odeset('RelTol', 1e-6);
switch solver
    case 'Runge Kutta adaptive'
        tstep=(end_time-start_time)/plotsteps;
        h=0.0001;
        a=initial_conditions.';
        t=start_time;
        Y=zeros(plotsteps,modes_number);
        for kk=1:plotsteps
            tmax=t+tstep;
            if tmax>end_time
                tmax=end_time;
            end
            while t<tmax
                if t+h>tmax 
                    h=tmax-t;
                end
                 [a,t,h]=rk5adapt(a,t,h);         
            end
            Y(kk,:)=a;
        end
    case 'ode45'        
        [~,Y] = ode45(@coupled_modes_equations,timepoints,initial_conditions,options);
    case 'ode23s'
        [~,Y] = ode23s(@coupled_modes_equations,timepoints,initial_conditions,options);
end
toc 

Y_wo_pump = Y;
Y_wo_pump(:,round(modes_number/2)) = zeros(size(Y_wo_pump,1),1);
dB=zeros(plotsteps,modes_number);
for kk=1:plotsteps
    dB(kk,:)=20*log10(abs(Y(kk,:))/max(abs(Y(kk,:))));
end
% dB=20*log10(abs(Y));
% dB=dB-max(max(dB))+100;
spectrum_plot=dB;
% spectrum_plot=abs(Y_wo_pump);
spectrum_plot_transposed=spectrum_plot.';

globals = who('global');
for kk = 1:numel(globals)
  eval(sprintf('global %s', globals{kk}));
end

if injection_force==0
    save([filename '.mat']);
else
    save([injection_filename '.mat']);
end
end

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global modes_number;
global fsr;
global linewidth;
global coupling;
global refr_index;
global nonlin_index;
global mode_area; 
global filename;
global snapshot;
global reslist;
global sweep_speed;
global pump_profile;
global pump;
global initial_conditions;
global initial;
global detuning_profile;
global timesteps_cme;
global plotsteps;
global progress;
global cavity_linewidths;
global pump_freq;

[filename,~,~] = uigetfile('*.mat');
file = load(filename);

modes_number =file.modes_number;
fsr=file.fsr;
linewidth=file.linewidth;
coupling=file.coupling;
refr_index=file.refr_index;
nonlin_index=file.nonlin_index;
mode_area=file.mode_area; 
reslist=file.reslist;
sweep_speed=file.sweep_speed;
pump_profile=file.pump_profile;
pump=file.pump;
initial_conditions=file.initial_conditions;
initial=file.initial;
detuning_profile=file.detuning_profile;
timesteps_cme=file.timesteps_cme;
plotsteps=file.plotsteps;
progress=file.progress;
cavity_linewidths=file.cavity_linewidths;
pump_freq=file.pump_freq;

set(handles.figure1,'Name',filename);
set(handles.linewidth,'String',num2str(cavity_linewidths(round(modes_number/2))));
set(handles.refr_index,'String',num2str(file.refr_index));
set(handles.nonlin_index,'String',num2str(file.nonlin_index));
set(handles.mode_area,'String',num2str(file.mode_area));
set(handles.coupling,'String',num2str(file.coupling));
set(handles.sweep_speed,'String',num2str(file.sweep_speed));

snapshot=0.5;
% adjust slider to snapshot
min_slider=get(handles.slider3,'Min');
max_slider=get(handles.slider3,'Max');
slider_value=snapshot*(max_slider-min_slider)+min_slider;
plotcomb(filename,snapshot,'all');
set(handles.slider3,'Enable','on');
set(handles.slider3,'Value',slider_value);
% set(handles.inject,'Enable','on');
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)
end

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)
end

function modes_number_Callback(hObject, eventdata, handles)
% hObject    handle to modes_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global modes_number;
modes_number=str2double(get(hObject,'String'));
% Hints: get(hObject,'String') returns contents of modes_number as text
%        str2double(get(hObject,'String')) returns contents of modes_number as a double
end

function linewidth_Callback(hObject, eventdata, handles)
% hObject    handle to linewidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global cavity_linewidths;
global modes_number;

cavity_linewidths=str2double(get(hObject,'String'))*ones(1,modes_number);
% Hints: get(hObject,'String') returns contents of linewidth as text
%        str2double(get(hObject,'String')) returns contents of linewidth as a double
end

% --- Executes during object creation, after setting all properties.
function linewidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to linewidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global linewidth;
global linewidth_handler;
linewidth_handler=hObject;
linewidth=str2double(get(hObject,'String'));
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function refr_index_Callback(hObject, eventdata, handles)
% hObject    handle to refr_index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global refr_index;
refr_index=str2double(get(hObject,'String'));
% Hints: get(hObject,'String') returns contents of refr_index as text
%        str2double(get(hObject,'String')) returns contents of refr_index as a double
end

% --- Executes during object creation, after setting all properties.
function refr_index_CreateFcn(hObject, eventdata, handles)
% hObject    handle to refr_index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global refr_index;
global refr_index_handler;
refr_index_handler=hObject;
refr_index=str2double(get(hObject,'String'));
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function nonlin_index_Callback(hObject, eventdata, handles)
% hObject    handle to nonlin_index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global nonlin_index;
nonlin_index=str2double(get(hObject,'String'));
% Hints: get(hObject,'String') returns contents of nonlin_index as text
%        str2double(get(hObject,'String')) returns contents of nonlin_index as a double
end

% --- Executes during object creation, after setting all properties.
function nonlin_index_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nonlin_index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global nonlin_index;
global nonlin_index_handler;
nonlin_index_handler=hObject;
nonlin_index=str2double(get(hObject,'String'));
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function mode_area_Callback(hObject, eventdata, handles)
% hObject    handle to mode_area (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global mode_area;
mode_area=str2double(get(hObject,'String'));
% Hints: get(hObject,'String') returns contents of mode_area as text
%        str2double(get(hObject,'String')) returns contents of mode_area as a double
end

% --- Executes during object creation, after setting all properties.
function mode_area_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mode_area (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global mode_area;
global mode_area_handler;
mode_area_handler=hObject;
mode_area=str2double(get(hObject,'String'));
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function coupling_Callback(hObject, eventdata, handles)
% hObject    handle to coupling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global coupling;
coupling=str2double(get(hObject,'String'));
% Hints: get(hObject,'String') returns contents of coupling as text
%        str2double(get(hObject,'String')) returns contents of coupling as a double
end

% --- Executes during object creation, after setting all properties.
function coupling_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coupling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global coupling_handler;
coupling_handler=hObject;
global coupling;    
coupling=str2double(get(hObject,'String'));
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
end

function sweep_speed_Callback(hObject, eventdata, handles)
% hObject    handle to sweep_speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sweep_speed;
sweep_speed=str2double(get(hObject,'String'));
% Hints: get(hObject,'String') returns contents of sweep_speed as text
%        str2double(get(hObject,'String')) returns contents of sweep_speed as a double
end

% --- Executes during object creation, after setting all properties.
function sweep_speed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sweep_speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global sweep_speed_handler;
sweep_speed_handler=hObject;
global sweep_speed;
sweep_speed=str2double(get(hObject,'String'));
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global filename;
global snapshot;

min_slider=get(hObject,'Min');
max_slider=get(hObject,'Max');
value_slider=get(hObject,'Value');
% snapshot=detune_s+(detune_e-detune_s)*(value_slider-min_slider)/(max_slider-min_slider);
snapshot=(value_slider-min_slider)/(max_slider-min_slider);
plotcomb(filename,snapshot,'all')
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
end

% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Enable','off');
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

% --- Executes on button press in edit_modes.
function edit_modes_Callback(hObject, eventdata, handles)
% hObject    handle to edit_modes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global modes_number;
global pump_freq;
global fsr;
global d2;
global d3;
global d4;
global d5;
global nms_a;
global nms_b;
global linewidth;
global reslist;
global c;
global cavity_linewidths;
global cavity_linewidths_string;

prompt={'Modes number',...
    'Wavelength (m)',...
    'FSR (Hz)', ...
    'D2 (Hz)', ...
    'D3 (Hz)', ...
    'D4 (Hz)', ...
    'D5 (Hz)', ...
    'Distortion, a', ...
    'Distortion, b'
    };
defaultanswer={num2str(modes_number),num2str(c/pump_freq),num2str(fsr),num2str(d2),num2str(d3),num2str(d4),num2str(d5),num2str(nms_a),num2str(nms_b)};
options.Resize='on';
options.WindowStyle='normal';
answer=inputdlg(prompt,'Dispersion',1,defaultanswer,options);
modes_number=str2double(answer{1});
% lambda=str2double(answer{2})*10^-9; % in nm
lambda=str2double(answer{2});
pump_freq=c/lambda;
fsr=str2double(answer{3});
d2=str2double(answer{4});
d3=str2double(answer{5});
d4=str2double(answer{6});
d5=str2double(answer{7});
nms_a=str2double(answer{8});
nms_b=str2double(answer{9});
if length(cavity_linewidths)~=modes_number
    display('Wrong cavity linewidths array length. Using default 1e6 Hz')
    cavity_linewidths_string='cavity_linewidths=1e6*ones(1,modes_number)';
    eval(cavity_linewidths_string);
%     cavity_linewidths=1e6*ones(1,modes_number);
end
reslist = eigenmodes(modes_number, pump_freq, fsr, d2, d3, d4, d5, nms_a,nms_b,cavity_linewidths(round(modes_number/2)));
end

% --- Executes on button press in import_eigenmodes.
function import_eigenmodes_Callback(hObject, eventdata, handles)
% hObject    handle to import_eigenmodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global modes_number;
global reslist;
global pump_freq;
global fsr;
global c;
[FileName,~,~] = uigetfile('*.*');
% format long g
% reslist=dlmread(FileName);
reslist = csvread(FileName);
modes_number = length(reslist);
% pump_freq=reslist(round(modes_number/2));
prompt={'FSR (Hz)', 'Wavelength (m)'};
defaultanswer={'200e9',num2str(c/pump_freq)};
options.Resize='on';
options.WindowStyle='normal';
answer=inputdlg(prompt,'FSR and Pump Wavelength',1,defaultanswer,options);
fsr=str2double(answer{1});
lambda=str2double(answer{2});
pump_freq=c/lambda;
end

% --- Executes on button press in pump_edit.
function pump_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pump_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pump;
global timesteps_cme;
global pump_string;
prompt={'Pump power (W)'};
defaultanswer={pump_string};
options.Resize='on';
options.WindowStyle='normal';
answer=inputdlg(prompt,'Pump profile',[1 50],defaultanswer,options);
pump=eval(answer{1});
pump_string=answer{1};
set(handles.pump_constant,'String','');
end

% --- Executes on button press in pump_import.
function pump_import_Callback(hObject, eventdata, handles)
% hObject    handle to pump_import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pump;
global pump_string;
[FileName,~,~] = uigetfile('*.*');
pump = csvread(FileName);
pump_string=mat2str(pump);
end

% --- Executes on button press in edit_detuning.
function edit_detuning_Callback(hObject, eventdata, handles)
% hObject    handle to edit_detuning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global detuning_profile;
global  timesteps_cme;
global detuning_string;
prompt={'Detuning'};
defaultanswer={detuning_string};
options.Resize='on';
options.WindowStyle='normal';
answer=inputdlg(prompt,'Detuning profile',[1 50],defaultanswer,options);
detuning_profile=eval(answer{1});
detuning_string=answer{1};
set(handles.detuning_start,'String','');
set(handles.detuning_end,'String','');
end

% --- Executes on button press in import_detuning.
function import_detuning_Callback(hObject, eventdata, handles)
% hObject    handle to import_detuning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global detuning_profile;
global detuning_string;
[FileName,~,~] = uigetfile('*.*');
detuning_profile = csvread(FileName);
detuning_string=mat2str(detuning_profile);
set(handles.detuning_start,'String','');
set(handles.detuning_end,'String','');
end

% --- Executes on button press in seeding_edit.
function seeding_edit_Callback(hObject, eventdata, handles)
% hObject    handle to seeding_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global initial;
global modes_number;
global initial_string;
prompt={'Initial conditions'};
defaultanswer={initial_string};
options.Resize='on';
options.WindowStyle='normal';
answer=inputdlg(prompt,'Seeding',[1 50],defaultanswer,options);
initial=eval(answer{1});
initial_string=answer{1};
end

% --- Executes on button press in seeding_import.
function seeding_import_Callback(hObject, eventdata, handles)
% hObject    handle to seeding_import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global initial;
global initial_string;
[FileName,~,~] = uigetfile('*.*');
initial = csvread(FileName);
initial_string=mat2str(initial);
end

% --- Executes on selection change in popup_plot.
function popup_plot_Callback(hObject, eventdata, handles)
% hObject    handle to popup_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_plot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_plot
global snapshot;
global filename;
global reslist;
% global pump_profile;
global pump;
global modes_number;
global detuning_profile;
% global initial_conditions;
global initial;
global cavity_linewidths;
global fsr;
global pump_freq;
% global d2;
% global d3;

modes_centered=(1:modes_number)-round(modes_number/2);
contents = cellstr(get(hObject,'String'));
switch contents{get(hObject,'Value')}
    case 'Eigenmodes'
        list=zeros(modes_number,1,'double');
        for k=1:modes_number
            list(k) = reslist(k)-double(int64(k-round(size(list,1)/2))*int64(fsr))-pump_freq;
        end
        figure
        plot(modes_centered,list)
        title ('deviation from linear grid with mean FSR');
        xlabel('mode number');
        ylabel('$\omega_\mu-\omega_0-\mu D_1 (Hz)$','interpreter','latex');
    case 'Pump Profile'
        figure
        plot(linspace(1,length(pump),length(pump)),pump)
%         plot(linspace(1,length(pump_profile),length(pump_profile)),pump_profile)
        title ('Pump Profile');
        xlabel('timestep');
        ylabel('W');
    case 'Initial Profile'
        figure
        plot(modes_centered,abs(initial))
%         plot(linspace(1,length(initial_conditions),length(initial_conditions)),abs(initial_conditions))
        title ('Initial Profile');
        xlabel('mode number');
        ylabel('a.u.');
    case 'Detuning Profile'
        figure
        plot(linspace(1,length(detuning_profile),length(detuning_profile)),detuning_profile)
        title ('Detuning Profile');
        xlabel('timestep');
        ylabel('linewidths');
    case 'Cavity Linewidths'
        figure
        plot(modes_centered,cavity_linewidths)
        title ('Cavity Linewidths');
        xlabel('mode number');
        ylabel('Hz');
    case 'Total Field'
        plotcomb(filename,snapshot,'total_field');        
    case 'Amplitudes'
        plotcomb(filename,snapshot,'amps');
    case 'Spectrum' 
        plotcomb(filename,snapshot,'spectrum');
    case 'Waveform'
        plotcomb(filename,snapshot,'waveform');
end
end

% --- Executes during object creation, after setting all properties.
function popup_plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function detuning_start_Callback(hObject, eventdata, handles)
% hObject    handle to detuning_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of detuning_start as text
%        str2double(get(hObject,'String')) returns contents of detuning_start as a double
global detuning_start;
detuning_start=str2double(get(hObject,'String'));

end

% --- Executes during object creation, after setting all properties.
function detuning_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to detuning_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function detuning_end_Callback(hObject, eventdata, handles)
% hObject    handle to detuning_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of detuning_end as text
%        str2double(get(hObject,'String')) returns contents of detuning_end as a double
global detuning_end;
detuning_end=str2double(get(hObject,'String'));
end

% --- Executes during object creation, after setting all properties.
function detuning_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to detuning_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function res = eigenmodes(modes_number,pump_freq,fsr,d2over2pi,d3over2pi,d4over2pi,d5over2pi,nms_a,nms_b,linewidth)
res = zeros(modes_number,1);
for k=1:size(res,1)
    res(k) = pump_freq ...
        + (k-round(size(res,1)/2))*fsr ...
        + (k-round(size(res,1)/2))^2 * d2over2pi/2 ...
        + (k-round(size(res,1)/2))^3 * d3over2pi/6 ...
        + (k-round(size(res,1)/2))^4 * d4over2pi/24 ...
        + (k-round(size(res,1)/2))^5 * d5over2pi/120 ...
        + nms_a*linewidth/4/(k-round(size(res,1)/2)-nms_b-0.5);
end
% res(round(size(res,1)/2))=res(round(size(res,1)/2))-5*linewidth;
end

% --- Executes on button press in inject.
function inject_Callback(hObject, eventdata, handles)
% hObject    handle to inject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global snapshot;
global detuning_profile;
global pump_profile;
global modes_number;
global injection_detuning;
global injection_force;
global initial_conditions;
global progress;
global injection_progress;
global injection_filename;
global d2;
global kappa;

injection_filename=strcat('injection',datestr(fix(clock),'yyyymmddHHMMSS'));
injection_progress=hObject;
set(progress,'Enable','off');
set(hObject,'Enable','off');
injection_detuning=detuning_profile(round(length(detuning_profile)*snapshot));
injection_force=pump_profile(round(length(pump_profile)*snapshot));
% see Supplementary information (nphoton.2013.343-s1)
theta=linspace(-pi,pi,modes_number);
zeta0=2*injection_detuning;
f=injection_force;
if ((2/27*zeta0*(zeta0^2+9)-f^2)<0) || (zeta0<sqrt(3))
    display('wrong solition approximation');
end
b=sqrt(zeta0/d2/2/pi*kappa);
% b=sqrt(zeta0*2);
% phi0=acos(sqrt(8*zeta0)/pi/f);
a=4*zeta0/pi/f+1i*sqrt(2*zeta0-(4*zeta0/pi/f)^2);
if abs(2*zeta0^3+18*zeta0-27*f^2)/(2*(zeta0^2-3)^(3/2))<1
    abspsi0sqr=2/3*(zeta0-sqrt(zeta0^2-3)*cos(1/3*acos((2*zeta0^3+18*zeta0-27*f^2)/(2*(zeta0^2-3)^(3/2)))));
else
    abspsi0sqr=2/3*(zeta0-sqrt(zeta0^2-3)*cosh(1/3*acosh((2*zeta0^3+18*zeta0-27*f^2)/(2*(zeta0^2-3)^(3/2)))));
end
psi0=1i*f/(abspsi0sqr-zeta0+1i);
% psi1=psi0+b*exp(1i*phi0)*sech(b*theta);
psi=psi0+a*sech(b*theta);
initial_conditions=ifftshift(ifft(psi));
simulate();
set(progress,'Enable','on');
set(injection_progress,'Enable','on');
set(hObject,'String','Inject Soliton');
plotcomb(injection_filename,0.01,'injection');
end

function [a_out,t_out,h_out]=rk5adapt(a,t,h)

b21=0.2; b31=0.075; b41=0.3; b51=-11/54; b61=1631/55296;
b32=0.225; b42=-0.9; b52=2.5; b62=175/512;
b43=1.2; b53=-70/27; b63=575/13824;
b54=35/27; b64=44275/110592;
b65=253/4096;
c01=37/378; c03=250/621; c04=125/594; c06=512/1771;
c11=2825/27648; c13=18575/48384; c14=13525/55296; c15=277/14336; c16=0.25;
myerr=1e-7;
hmin=0.001;

k1=h*coupled_modes_equations(t,a);
atmp=a+b21*k1;
k2=h*coupled_modes_equations(t+0.2*h,atmp);
atmp=a+b31*k1+b32*k2;	
k3=h*coupled_modes_equations(t+0.3*h,atmp);
atmp=a+b41*k1+b42*k2+b43*k3;
k4=h*coupled_modes_equations(t+0.6*h,atmp);
atmp=a+b51*k1+b52*k2+b53*k3+b54*k4;
k5=h*coupled_modes_equations(t+h,atmp);
atmp=a+b61*k1+b62*k2+b63*k3+b64*k4+b65*k5;
k6=h*coupled_modes_equations(t+0.875*h,atmp);
z1=a+c11*k1+c13*k3+c14*k4+c15*k5+c16*k6;
atmp=a+c01*k1+c03*k3+c04*k4+c06*k6;
err=max(abs(atmp-z1));
if (myerr>=err || h <= hmin)
    h_out=h*0.99*(myerr/err)^0.2;
    t_out=t+h_out;
    a_out=atmp;
else
    hfactor=0.99*(myerr/err)^0.25;
    if hfactor>0.1*h
        h=h*hfactor;
    else
        h=h*0.1;
    end
    [a_out,t_out,h_out]=rk5adapt(a,t,h);
end
end


% --- Executes on selection change in solver.
function solver_Callback(hObject, eventdata, handles)
% hObject    handle to solver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns solver contents as cell array
%        contents{get(hObject,'Value')} returns selected item from solver
global solver;
contents = cellstr(get(hObject,'String'));
solver=contents{get(hObject,'Value')};
end

% --- Executes during object creation, after setting all properties.
function solver_CreateFcn(hObject, eventdata, handles)
% hObject    handle to solver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in cavity_linewdth.
function cavity_linewdth_Callback(hObject, eventdata, handles)
% hObject    handle to cavity_linewdth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global cavity_linewidths;
global modes_number;
global cavity_linewidths_string;
prompt={'Cavity linewiths (Hz)'};
defaultanswer={cavity_linewidths_string};
options.Resize='on';
options.WindowStyle='normal';
answer=inputdlg(prompt,'Cavity linewiths',[1 50],defaultanswer,options);
eval(answer{1});
cavity_linewidths_string=answer{1};
set(handles.linewidth,'String','');
end

function pump_constant_Callback(hObject, eventdata, handles)
% hObject    handle to pump_constant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pump_constant as text
%        str2double(get(hObject,'String')) returns contents of pump_constant as a double
global pump;
global timesteps_cme;
pump=str2double(get(hObject,'String'))*ones(1,timesteps_cme);
end

% --- Executes during object creation, after setting all properties.
function pump_constant_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pump_constant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
global pump_constant;
global pump_constant_handler;
pump_constant_handler=hObject;
pump_constant=str2double(get(hObject,'String'));

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --------------------------------------------------------------------
function open_parameters_Callback(hObject, eventdata, handles)
% hObject    handle to open_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global modes_number;
global fsr;
global linewidth;
global coupling;
global refr_index;
global nonlin_index;
global mode_area; 
% global snapshot;
% global reslist;
global sweep_speed;
% global pump_profile;
% global pump;
% global initial_conditions;
% global initial;
% global detuning_profile;
% global timesteps_cme;
% global plotsteps;
% global progress;
% global cavity_linewidths;
global pump_freq;
global d2;
global d3;
global d4;
global d5;
global nms_a;
global nms_b;

[filename,~,~] = uigetfile('*.txt');
tbl = readtable(filename);
% summary(tbl)
modes_number =tbl{1,2};
fsr=tbl{2,2};
pump_freq=tbl{3,2};
d2=tbl{4,2};
d3=tbl{5,2};
d4=tbl{6,2};
d5=tbl{7,2};
nms_a=tbl{8,2};
nms_b=tbl{9,2};
linewidth=tbl{10,2};
coupling=tbl{11,2};
refr_index=tbl{12,2};
nonlin_index=tbl{13,2};
mode_area=tbl{14,2}; 
sweep_speed=tbl{15,2};

set(handles.linewidth,'String',num2str(linewidth));
set(handles.refr_index,'String',num2str(refr_index));
set(handles.nonlin_index,'String',num2str(nonlin_index));
set(handles.mode_area,'String',num2str(mode_area));
set(handles.coupling,'String',num2str(coupling));
set(handles.sweep_speed,'String',num2str(sweep_speed));
end