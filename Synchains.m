function varargout = Synchains(varargin)
% SYNCHAINS MATLAB code for synchains.fig
%      SYNCHAINS, by itself, creates a new SYNCHAINS or raises the existing
%      singleton*.
%
%      H = SYNCHAINS returns the handle to a new SYNCHAINS or the handle to
%      the existing singleton*.
%
%      SYNCHAINS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SYNCHAINS.M with the given input arguments.
%
%      SYNCHAINS('Property','Value',...) creates a new SYNCHAINS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before synchains_openingfcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to synchains_openingfcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help synchains

% Last Modified by GUIDE v2.5 22-Jul-2013 12:55:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Synchains_OpeningFcn, ...
                   'gui_OutputFcn',  @Synchains_OutputFcn, ...
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


% --- Executes just before synchains is made visible.
function Synchains_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to synchains (see VARARGIN)
set(handles.figure1,'CloseRequestFcn',@closeGUI);

% Choose default command line output for synchains
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes synchains wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Synchains_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function MC_persistenceL_Callback(hObject, eventdata, handles)
% hObject    handle to MC_persistenceL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MC_persistenceL as text
%        str2double(get(hObject,'String')) returns contents of MC_persistenceL as a double


% --- Executes during object creation, after setting all properties.
function MC_persistenceL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MC_persistenceL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_MC_randi.
function checkbox_MC_randi_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_MC_randi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_MC_randi

if (get(hObject,'Value') == get(hObject,'Max'))
   % Checkbox is checked-take appropriate action
   handles.checkbox_MC_randi = 1;
else
   % Checkbox is not checked-take appropriate action
   handles.checkbox_MC_randi = 0;
end
guidata(hObject,handles);


function MC_contourL_mean_Callback(hObject, eventdata, handles)
% hObject    handle to MC_contourL_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MC_contourL_mean as text
%        str2double(get(hObject,'String')) returns contents of MC_contourL_mean as a double


% --- Executes during object creation, after setting all properties.
function MC_contourL_mean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MC_contourL_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MC_contourL_SD_Callback(hObject, eventdata, handles)
% hObject    handle to MC_contourL_SD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MC_contourL_SD as text
%        str2double(get(hObject,'String')) returns contents of MC_contourL_SD as a double


% --- Executes during object creation, after setting all properties.
function MC_contourL_SD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MC_contourL_SD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MC_segmentL_mean_Callback(hObject, eventdata, handles)
% hObject    handle to MC_segmentL_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MC_segmentL_mean as text
%        str2double(get(hObject,'String')) returns contents of MC_segmentL_mean as a double


% --- Executes during object creation, after setting all properties.
function MC_segmentL_mean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MC_segmentL_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MC_segmentL_SD_Callback(hObject, eventdata, handles)
% hObject    handle to MC_segmentL_SD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MC_segmentL_SD as text
%        str2double(get(hObject,'String')) returns contents of MC_segmentL_SD as a double


% --- Executes during object creation, after setting all properties.
function MC_segmentL_SD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MC_segmentL_SD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MC_linewidth_Callback(hObject, eventdata, handles)
% hObject    handle to MC_linewidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MC_linewidth as text
%        str2double(get(hObject,'String')) returns contents of MC_linewidth as a double


% --- Executes during object creation, after setting all properties.
function MC_linewidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MC_linewidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function txt_MCfilename_Callback(hObject, eventdata, handles)
% hObject    handle to txt_MCfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_MCfilename as text
%        str2double(get(hObject,'String')) returns contents of txt_MCfilename as a double


% --- Executes during object creation, after setting all properties.
function txt_MCfilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_MCfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_DoHisto.
function checkbox_DoHisto_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_DoHisto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_DoHisto
if (get(hObject,'Value') == get(hObject,'Max'))
   % Checkbox is checked-take appropriate action
   handles.checkbox_DoHisto = 1;
else
   % Checkbox is not checked-take appropriate action
   handles.checkbox_DoHisto = 0;
end
guidata(hObject,handles);



function MC_Nfib_Callback(hObject, eventdata, handles)
% hObject    handle to MC_Nfib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MC_Nfib as text
%        str2double(get(hObject,'String')) returns contents of MC_Nfib as a double


% --- Executes during object creation, after setting all properties.
function MC_Nfib_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MC_Nfib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MC_hist_Nbins_Callback(hObject, eventdata, handles)
% hObject    handle to MC_hist_Nbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MC_hist_Nbins as text
%        str2double(get(hObject,'String')) returns contents of MC_hist_Nbins as a double


% --- Executes during object creation, after setting all properties.
function MC_hist_Nbins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MC_hist_Nbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_margin.
function checkbox_margin_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_margin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_margin
if (get(hObject,'Value') == get(hObject,'Max'))
   % Checkbox is checked-take appropriate action
   handles.checkbox_margin = 1;
else
   % Checkbox is not checked-take appropriate action
   handles.checkbox_margin = 0;
end
guidata(hObject,handles);



function set_margins_Callback(hObject, eventdata, handles)
% hObject    handle to set_margins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_margins as text
%        str2double(get(hObject,'String')) returns contents of set_margins as a double


% --- Executes during object creation, after setting all properties.
function set_margins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_margins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbut_aboutSyn.
function pushbut_aboutSyn_Callback(hObject, eventdata, handles)
msgbox(sprintf('Guillaume Lamour \n*******************\nUniversity of British Columbia, Canada\n*******************\nlamour99@hotmail.com \n*******************\nEnjoy\n*******************\n2013 '),...
       'About Synchains','Help')







function pushbut_MC_printALL_Callback(hObject, eventdata, handles)
% hObject    handle to pushbut_MC_printALL (see GCBO)

Nfibrils = str2double(get(handles.MC_Nfib,'String'));
P = str2double(get(handles.MC_persistenceL,'String'));
contour = str2double(get(handles.MC_contourL_mean,'String'));
contour_SD = str2double(get(handles.MC_contourL_SD,'String'));
A = normrnd(contour,contour_SD,[1 200]);
segment = str2double(get(handles.MC_segmentL_mean,'String'));
segment_SD = str2double(get(handles.MC_segmentL_SD,'String'));
B = normrnd(segment,segment_SD,[1 200]);
linewidth = str2double(get(handles.MC_linewidth,'String'));

figure
checkboxMC = handles.checkbox_MC_randi;


j = 1;
while j <= Nfibrils
    % Here we select contour length at random in one symmetric interval
    % (1st condition), box checked) or around the mean given the normal
    % distribution  (second)
    if checkboxMC == 1;
        contour = A(randi(200));
    else
        % IT IS essential to update that stuff below
        contour = str2double(get(handles.MC_contourL_mean,'String'));
        contour_SD = str2double(get(handles.MC_contourL_SD,'String'));
        imin = contour - contour_SD;
        if imin <= 0;
            imin = 1;
        else
        end
        imax = contour + contour_SD;
        contour = randi([imin, imax]);
    end
    seg = segment;
    if contour >= seg
        n = round(contour./seg)+1; %define the number of random angles
        sigma = (sqrt(seg./P));
        % generate random numbers with gaussian with standard deviation
        % given by persistence length and segment length
        %R = [];
        R = normrnd(0,sigma,[n 1]);  % R = randi(180,[n 1]);
        seg1 = B(randi(200));
        seg2 = B(randi(200));
        spline_x = [];
        spline_y = [];
        spline_x(1)=0;
        spline_y(1)=0;
        spline_x(2) =seg1;
        spline_y(2) =0;
        spline_x(3) = seg1 + seg2 * cos(R(1,1));
        spline_y(3) = seg2 * sin (R(1,1));
        %i = 4;
        for i = 4:(n+2)
            seg = B(randi(200));
            spline_x(i)= seg * cos(sum(R(1:(i-2)))) + spline_x(i-1);
            spline_y(i)= seg * sin(sum(R(1:(i-2)))) + spline_y(i-1);
        end
    else

    end
    plot(spline_x,spline_y,'LineWidth',linewidth)
    xlabel('nm'); %Write label for x-axis
    ylabel('nm'); %Write label for y-axis
    axis equal
    hold on
    %Store contour length of generated fibrils to be displayed in histogram
    %later on:
    struc_contour(j) = contour; 
    % Store coordinates to adjust graph axes
    minx(j) = min(spline_x);
    miny(j) = min(spline_y);
    maxx(j) = max(spline_x);
    maxy(j) = max(spline_y);
    j = j + 1;
end

checkboxMARGIN = handles.checkbox_margin;
margin = str2double(get(handles.set_margins,'String'));

if checkboxMARGIN ==1;
    % Here we normalize distances to make the graph look good
    minx = min(minx)-margin; maxx = max(maxx)+ margin;
    miny = min(miny)-margin; maxy = max(maxy)+ margin;
    minx = round2(minx,margin); maxx = round2(maxx,margin);
    miny = round2(miny,margin); maxy = round2(maxy,margin);
    midpointx = (maxx - minx)./2 + minx; lengthx = (maxx - minx)./2;
    midpointy = (maxy - miny)./2 + miny; lengthy = (maxy - miny)./2;
    if maxx - minx < maxy - miny;
        minx = midpointx - abs(lengthy);
        maxx = midpointx + abs(lengthy);
    else
        miny = midpointy - abs(lengthx);
        maxy = midpointy + abs(lengthx);
    end
    axis([ minx maxx miny maxy])
else
end

checkboxMC_histo = handles.checkbox_DoHisto;
if checkboxMC_histo == 1;
    mat = struc_contour ;
    nbins = str2double(get(handles.MC_hist_Nbins,'String'));
    figure
    hist(mat,nbins)
    xlabel('Contour Length of Monte-Carlo Polymers (nm)'); %Write label for x-axis
    ylabel('# of events'); %Write label for y-axis
else
end


guidata(hObject,handles);




% --- Executes on button press in pushbut_MC_generateData.
function pushbut_MC_generateData_Callback(hObject, eventdata, handles)

  
Nfib = str2double(get(handles.numchain_generate,'String')); %set
P = str2double(get(handles.MC_persistenceL,'String')); %set (in nm)
contour = str2double(get(handles.MC_contourL_mean,'String'));  %set (in nm)
contour_SD = str2double(get(handles.MC_contourL_SD,'String'));  %set (in nm)
checkboxMC = handles.checkbox_MC_randi; %set to 1 only if you want normal distribution of contour length
%set to any other number will make the distribution comprised between intervals +/ contour_SD
segment = str2double(get(handles.MC_segmentL_mean,'String'));  %set (in nm)
segment_SD = str2double(get(handles.MC_segmentL_SD,'String'));  %set (in nm)

AAA = normrnd(contour,contour_SD,[1 200]); %don't change this
BBB = normrnd(segment,segment_SD,[1 200]);  %don't change this

    bsplines_struct = struct([]);
    bsplines_norm2pl = struct([]);
    mat_length_struct = struct([]);
    mat_intervals_struct = struct([]);
    mat_sd_struct = struct([]);
    deviations_struct = struct([]);
    correlations_struct = struct([]);
    wormlike_struct = struct([]);
%   angle_struct = struct([]);
    
kappa= 1; %kappa= fibril index

while kappa<= Nfib
    
          
    %figure
    if checkboxMC == 1;
        contour = AAA(randi(200));  %don't change this
    else
        contour = randi([contour - contour_SD, contour + contour_SD]);
    end
    
    seg = segment;
    if contour >= seg
        n = round(contour./seg); %define the number of random angles
        sigma = (sqrt(seg./P));
        % generate random numbers with gaussian with standard deviation
        % given by persistence length and segment length
        R = [];
        R = normrnd(0,sigma,[n 1]);  
           
        seg1 = BBB(randi(200));
        seg2 = BBB(randi(200));
        spline_x = [];
        spline_y = [];
        spline_x(1)=0;
        spline_y(1)=0;
        spline_x(2) =seg1;
        spline_y(2) =0;
        spline_x(3) = seg1 + seg2 * cos(R(1,1));
        spline_y(3) = seg2 * sin (R(1,1));
%                       angledata = [];
%                       angledata(1) = R(1,1);

%THat loop is writing the coordinates of one synthetic chain
        ibis = 4;
        for ibis = 4:(n+2)
            seg = BBB(randi(200));
             spline_x(ibis)= seg * cos(sum(R(1:(ibis-2)))) + spline_x(ibis-1);
             spline_y(ibis)= seg * sin(sum(R(1:(ibis-2)))) + spline_y(ibis-1);
            % angledata = [angledata; R(ibis-2)];       
        end
    else
    end
    
    %-----------------------------------------------------------------------
    %calculate intervals
    x_norm = [];    x_norm = spline_x;
    y_norm = [];    y_norm = spline_y;
    Nseg = length(x_norm)-1; % = number of segments 
    intervals = [];
    
    for ippi=1:Nseg
        dist = ( x_norm(ippi) - x_norm(ippi+1) ).^2 + ( y_norm(ippi) - y_norm(ippi+1) ).^2;
        intervals(ippi) = sqrt(dist);
    end
    
    % average_interval between each knot over all the spline
    mean_itv = []; 
    mean_itv = mean (intervals(1,:));
    mean_itv = round (mean_itv);
    % standard deviation
    sd_itv = [];
    sd_itv = std(intervals(1,:));
    sd_itv = round (sd_itv);
    
    %-----------------------------------------------------------------------
    %calculate combo
        
    %MIDPOINT-FLUCT
    %recover data from previous functions and pass them through this one
    midpointX = 0;
    midpointY = 0;
    sec_length = 0;
    lastpoint = length(x_norm);
    
        delta = []; secant = [];
        deltas = []; secants = [];
        A = []; B = []; C = [];
    
        for n = 2:lastpoint - 1
            for j = 1:lastpoint - n
                midpointX = (x_norm(j) + x_norm(j+n))/2;
                midpointY = (y_norm(j) + y_norm(j+n))/2;
                sec_length = sqrt (( y_norm(j+n) - y_norm(j) ).^2 +...
                    (( x_norm(j+n) - x_norm(j) ).^2 ));
                shortest_dist = 1e+10;
                for i = 1:lastpoint
                    dist = (( midpointY - y_norm(i) ).^2 +...
                        (( midpointX - x_norm(i) ).^2));
                    dist = sqrt (dist);
                    if (dist <= shortest_dist)
                        shortest_dist = dist;
                        % Keep track of closest point
                        closest_i = i;
                    end
                end
                %%% Refine the shortest distance
                if closest_i ~= 1 & closest_i <= lastpoint -1
                    ci = closest_i;
                    % Avec, Bvec, Cvec are the length of vectors of a triangle formed between
                    % closest and midpoint, midpoint and next to closest, closest and next to closest,
                    % respectively
                    Avec = sqrt ( (x_norm(ci) - midpointX).^2 + (y_norm(ci) - midpointY).^2 );
                    Bvec = sqrt ( (x_norm(ci+1) - midpointX).^2 + (y_norm(ci+1) - midpointY).^2 );
                    Cvec = sqrt ( (x_norm(ci+1) - x_norm(ci)).^2 + (y_norm(ci+1) - y_norm(ci)).^2 );
                    % theta is the angle between A and C
                    cos_theta1 = (Avec.^2 + Cvec.^2 - Bvec.^2) ./ (2 * Avec * Cvec);
                    theta1 = acos(cos_theta1);
                    shortest_dist1 = Avec * sin(theta1);
                    Bvec2 = sqrt ( (x_norm(ci-1) - midpointX).^2 + (y_norm(ci-1) - midpointY).^2 );
                    Cvec2 = sqrt ( (x_norm(ci-1) - x_norm(ci)).^2 + (y_norm(ci-1) - y_norm(ci)).^2 );
                    cos_theta2 = (Avec.^2 + Cvec2.^2 - Bvec2.^2) ./ (2 * Avec * Cvec2);
                    theta2 = acos(cos_theta2);
                    shortest_dist2 = Avec * sin(theta2);
                    if shortest_dist1 <= shortest_dist2
                        shortest_dist = shortest_dist1;
                    else
                        shortest_dist = shortest_dist2;
                    end
                end
                delta(j)= shortest_dist;
                secant(j) = sec_length;
            end
            deltas(:,n-1) = delta;
            secants(:,n-1) = secant;
        end
    % pick the 2 matrices and make one in 2 dimensions
    A = deltas;
    B = secants;
    p = length(A);
    j = 0;
    for n = 0:(p - 1)
        for k = 1:(p - n)
            u = k + (n*p - j);
            C(u, 1) = A(k, n + 1 );
            C(u, 2) = B(k, n + 1 );
        end
        j = j + n;
    end
    deviat_single = C;
    
    %TANTAN-COREL
    
    tantan_corel = []; contour_length1 = [];
    cosines = []; contours = [];
    D = []; E = []; F = [];
    
    for j = 1:lastpoint - 2
        for i = 1:lastpoint - j - 1
            spX = x_norm(i+1) - x_norm(i);
            spY = y_norm(i+1) - y_norm(i);
            spnextX = x_norm(i+j+1) - x_norm(i+j);
            spnextY = y_norm(i+j+1) - y_norm(i+j);
            tan_i = [ spX, spY ];
            tan_k = [ spnextX, spnextY ];
            scal_prod = dot(tan_i, tan_k) / ( (norm(tan_i)) * (norm(tan_k)) );
            tantan_corel(i) = scal_prod;
            contour_length1(i) = sum(intervals(1,i:i+j-1));
            
        end
        cosines(:,j) = tantan_corel;
        contours(:,j) = contour_length1;
    end
    % pick the 2 matrices and make one in 2 dimensions
    D = cosines;
    E = contours;
    p = length(D);
    j = 0;
    for n = 0:(p - 1)
        for k = 1:(p - n)
            u = k + (n*p - j);
            F(u, 1) = D(k, n + 1 );
            F(u, 2) = E(k, n + 1 );
        end
        j = j + n;
    end
    corel_single = F;
    
    %CONTOUR-ENDEND
    
    end2end_length = []; end2end = [];
    G = []; H = []; I = [];
    
    for j = 1:lastpoint - 2
        for i = 1:lastpoint - j - 1
            end2end_length(i) = sqrt (( y_norm(i+j) - y_norm(i) ).^2 +...
                (( x_norm(i+j) - x_norm(i) ).^2 ));
            contour_length2(i) = sum(intervals(1,i:i+j-1));
        end
        end2end(:,j) = end2end_length;
        %contour(:,j) = contour_length2;
    end
    % pick the 2 matrices and make one in 2 dimensions
    G = end2end;
    H = contours;
    p = length(G);
    j = 0;
    for n = 0:(p - 1)
        for k = 1:(p - n)
            u = k + (n*p - j);
            I(u, 1) = G(k, n + 1 );
            I(u, 2) = H(k, n + 1 );
        end
        j = j + n;
    end
    end2cont_single = I;
    
    % extract the maximum contour length calculated, ie. the fibril length
    goodcontour = max(I(:,2));
    
    %--------------------------------------------------------------------
    %BELOW is the code that fill in the structure with all the elements, pertaining
    % to one fibril, calculated above in this very same function
    bspline_coord = [];
    bspline_coord(1,:) = x_norm;
    bspline_coord(2,:) = y_norm;
    
    bsplines_struct(kappa).bspline = bspline_coord;
%    angle_struct(kappa).ang = angledata; 
    
    mat_length_struct(kappa).contourL = goodcontour;
    mat_intervals_struct(kappa).meanitv_fib = mean_itv;
    mat_sd_struct(kappa).meansd_fib = sd_itv;
    deviations_struct(kappa).deviafib = deviat_single;
    correlations_struct(kappa).corelfib = corel_single;
    wormlike_struct(kappa).wormfib = end2cont_single;
   
   % hello = kappa
    set(handles.NsyntheticFib,'String',kappa) ;
    drawnow expose
    
    %Go to next fibril
    kappa= kappa+1;
    
end

contour_lengths = mat_length_struct;
intervals_means = mat_intervals_struct;
intervals_sd = mat_sd_struct;

bsplines_norm = bsplines_struct;
bsplines_norm2pl = 0;
%angles = angle_struct;

mat_deviations = deviations_struct;
mat_correlations = correlations_struct;
mat_wormlike = wormlike_struct;

filename = get(handles.txt_MCfilename,'String');
output_base = filename;
output_sup = '_syn';
updated_filename=[output_base, output_sup];
save(updated_filename,'contour_lengths','intervals_means',...
    'intervals_sd','bsplines_norm','bsplines_norm2pl','mat_deviations',...
    'mat_correlations','mat_wormlike');

% message confirmation
msgbox(sprintf('A .mat file has been successfully generated.  Please open it using Easyworm2.'),...
       'Confirmation Note','Help')
   
set(handles.MC_Nfib,'string',Nfib);
   
guidata(hObject,handles);


% --- Executes on button press in pushbut_MC_PixandExp.
function pushbut_MC_PixandExp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbut_MC_PixandExp (see GCBO)
filename = get(handles.set_imagename,'String');

A = imread(filename);
A = A(:,:,2); %extract blue channel
A = 1./A; %invert all elements in the matrix


output_sup = '.txt';
updated_filename=[filename, output_sup];
dlmwrite(updated_filename, A, 'delimiter', '\t');

% message confirmation
msgbox(sprintf('The image was converted. Note that this converter works only for single chains printed via Synchains.'),...
       'Confirmation Note','Help')



function set_imagename_Callback(hObject, eventdata, handles)
% hObject    handle to set_imagename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_imagename as text
%        str2double(get(hObject,'String')) returns contents of set_imagename as a double


% --- Executes during object creation, after setting all properties.
function set_imagename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_imagename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numchain_generate_Callback(hObject, eventdata, handles)
% hObject    handle to numchain_generate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numchain_generate as text
%        str2double(get(hObject,'String')) returns contents of numchain_generate as a double


% --- Executes during object creation, after setting all properties.
function numchain_generate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numchain_generate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function closeGUI(src,evnt)
%this function is called when the user attempts to close the GUI window

%this command brings up the close window dialog
selection = questdlg('Close Synchains?',...
                     'Closing Request',...
                     'Yes','No','No');

%if user clicks yes, the GUI window is closed, otherwise
%nothing happens
switch selection,
   case 'Yes',
    delete(gcf)
   case 'No'
     return
end
