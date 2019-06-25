% Pop up Easyworm2 panel 
function varargout = Easyworm2(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @analysis_new_OpeningFcn, ...
                   'gui_OutputFcn',  @analysis_new_OutputFcn, ...
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

% --- Executes just before analysis_new is made visible.
function analysis_new_OpeningFcn(hObject, eventdata, handles, varargin)
set(handles.figure1,'CloseRequestFcn',@closeGUI);
handles.output = hObject;

    set(handles.uipanel_contourE2E,'Visible','Off') 
    set(handles.uipanel_corel,'Visible','Off')
    set(handles.uipanel_deviations,'Visible','Off')
    set(handles.uipanel_inputdata,'Visible','On')
    set(handles.uipanel_plot,'Visible','Off')
    set(handles.uipanel_quickcalc,'Visible','Off') 
    set(handles.uipanel_surfparam,'Visible','Off') 
    set(handles.uipanel_elasticmod,'Visible','Off')  
    set(handles.uipanel_lengths,'Visible','Off')
    set(handles.uipanel_saveoutputs,'Visible','Off')
    
guidata(hObject, handles);

% if strcmp(get(hObject,'Visible'),'off')
%end

% --- Outputs from this function are returned to the command line.
function varargout = analysis_new_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
printdlg(handles.figure1)



% --- Executes on button press in btn_loadsample.
function btn_loadsample_Callback(hObject, eventdata, handles)
% hObject    handle to btn_loadsample (see GCBO)

struc = uiimport

handles.struc.mat_wormlike = struc.mat_wormlike;
handles.struc.mat_correlations = struc.mat_correlations;
handles.struc.mat_deviations = struc.mat_deviations;
handles.struc.bsplines_norm = struc.bsplines_norm;
handles.struc.bsplines_norm2pl = struc.bsplines_norm2pl;
handles.struc.contour_lengths = struc.contour_lengths;
handles.struc.intervals_means = struc.intervals_means;
handles.struc.intervals_sd = struc.intervals_sd;

%handles.struc.angles = struc.angles;

ntot_fibrils = length(struc.mat_wormlike);
handles.ntot_fibrils = ntot_fibrils;
set(handles.ans_nb_fibrils,'String',ntot_fibrils);
set(handles.Nfib_ploti,'String',ntot_fibrils);

set(handles.ans_nfibrilsworm,'String',ntot_fibrils);
set(handles.ans_nfibrilscorel,'String',ntot_fibrils);
set(handles.ans_nfibrilsdevia,'String',ntot_fibrils);

convcont(hObject,handles);

rgm_worm = str2double(get(handles.set_NB_of_RGM_worm,'String'));
handles.rgm_worm = rgm_worm;
bins_worm = str2double(get(handles.set_nb_bins_worm,'String'));
handles.bins_worm = bins_worm;
upmax_worm = str2double(get(handles.set_nb_upmax_worm,'String'));
handles.upmax_worm = upmax_worm;
lowfit_worm = str2double(get(handles.set_nb_lowfit_worm,'String'));
handles.lowfit_worm = lowfit_worm;
upfit_worm = str2double(get(handles.set_nb_upfit_worm,'String'));
handles.upfit_worm = upfit_worm;
figpop_worm = str2double(get(handles.set_nb_figpop_worm,'String'));
handles.figpop_worm = figpop_worm;

rgm_corel = str2double(get(handles.set_NB_of_RGM_corel,'String'));
handles.rgm_corel = rgm_corel;
bins_corel = str2double(get(handles.set_nb_bins_corel,'String'));
handles.bins_corel = bins_corel;
upmax_corel = str2double(get(handles.set_nb_upmax_corel,'String'));
handles.upmax_corel = upmax_corel;
lowfit_corel = str2double(get(handles.set_nb_lowfit_corel,'String'));
handles.lowfit_corel = lowfit_corel;
upfit_corel = str2double(get(handles.set_nb_upfit_corel,'String'));
handles.upfit_corel = upfit_corel;
figpop_corel = str2double(get(handles.set_nb_figpop_corel,'String'));
handles.figpop_corel = figpop_corel;

rgm_devia = str2double(get(handles.set_NB_of_RGM_devia,'String'));
handles.rgm_devia = rgm_devia;
bins_devia = str2double(get(handles.set_nb_bins_devia,'String'));
handles.bins_devia = bins_devia;
upmax_devia = str2double(get(handles.set_nb_upmax_devia,'String'));
handles.upmax_devia = upmax_devia;
lowfit_devia = str2double(get(handles.set_nb_lowfit_devia,'String'));
handles.lowfit_devia = lowfit_devia;
upfit_devia = str2double(get(handles.set_nb_upfit_devia,'String'));
handles.upfit_devia = upfit_devia;
figpop_devia = str2double(get(handles.set_nb_figpop_devia,'String'));
handles.figpop_devia = figpop_devia;


    set(handles.uipanel_contourE2E,'Visible','On') 
    set(handles.uipanel_corel,'Visible','On')
    set(handles.uipanel_deviations,'Visible','On')

guidata(hObject,handles);


function convcont(hObject,handles)

matwo = handles.struc.mat_wormlike;
ntot_fibrils = handles.ntot_fibrils;
new = [];
for i = 1:ntot_fibrils
    new = [new; matwo(1,i).wormfib]; %#ok<AGROW>
end
vecworm = new;
totaldatapoints = length(vecworm);
set(handles.ans_totdatapoints,'String',totaldatapoints);

x = handles.ntot_fibrils;
contour_lengths = handles.struc.contour_lengths;
intervals_means = handles.struc.intervals_means;
intervals_sd = handles.struc.intervals_sd;

for i = 1:x
    contourL_vec(1,i) = contour_lengths(1,i).contourL;
    intervalmeans_vec(1,i) = intervals_means(1,i).meanitv_fib;
    intervalstd_vec(1,i) = intervals_sd(1,i).meansd_fib;
end

contourL_mean = round(mean(contourL_vec));
contourL_std = round(std(contourL_vec));
sum_contL = sum(contourL_vec);
sum_contL = round2(sum_contL/1000,0.1);

maxcontourL = round(max(contourL_vec)); 
set(handles.ans_longestfibL,'String',maxcontourL);
%set(handles.txt_hgfit_uplim,'String',maxcontourL);

intervalmeans_mean = round(mean(intervalmeans_vec));
intervalmeans_std = round(std(intervalmeans_vec));
intervalstd_mean = round(mean(intervalstd_vec));
intervalstd_std = round(std(intervalstd_vec));

set(handles.ans_sumcontourl,'String',sum_contL);
set(handles.ans_fiblength_mean,'String',contourL_mean);
set(handles.ans_fiblength_std,'String',contourL_std);
set(handles.ans_knotintmean_mean,'String',intervalmeans_mean);
set(handles.ans_knotintmean_std,'String',intervalmeans_std);
%set(handles.ans_knotintstd_mean,'String',intervalstd_mean);
%set(handles.ans_knotintstd_std,'String',intervalstd_std);

handles.contourL_vec = contourL_vec;


guidata(hObject,handles);



function Nfib_ploti_Callback(hObject, eventdata, handles)
% hObject    handle to Nfib_ploti (see GCBO)

% --- Executes during object creation, after setting all properties.
function Nfib_ploti_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nfib_ploti (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in checkbox_randomfib.
function checkbox_randomfib_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_randomfib (see GCBO)
if (get(hObject,'Value') == get(hObject,'Max'))
   % Checkbox is checked-take appropriate action
   handles.checkbox_randomfib = 1;
else
   % Checkbox is not checked-take appropriate action
   handles.checkbox_randomfib = 0;
end
guidata(hObject,handles);


% --- Executes on button press in btn_pltfibrilsround.
function btn_pltfibrilsround_Callback(hObject, eventdata, handles)
% hObject    handle to btn_pltfibrilsround (see GCBO)
plotiround(hObject,handles);


function plotiround(hObject,handles)

spline_struc = handles.struc.bsplines_norm2pl;
checksize = size(spline_struc,2);

if checksize == 1
    msgbox(sprintf('The smoothed traces are not available for synthetic chains. Please hit the Raw Trace button.'),...
        'Wrong Button','Warn')
else

%Nth_knot = str2double(get(handles.txt_nthknot,'String'))+1;
Nth_knot = 1+1;
threshold_min = str2double(get(handles.txt_threshold_min,'String'));
threshold_max = str2double(get(handles.txt_threshold_max,'String'));
if threshold_max == 0
    threshold_max = 1e10;
else
end

linewidth = str2double(get(handles.txt_linewidth,'String'));
%Nfibtot = handles.ntot_fibrils;
Nfibtot = str2double(get(handles.ans_nb_fibrils,'String'));
n = str2double(get(handles.Nfib_ploti,'String'));

contour_lengths = handles.struc.contour_lengths;
for i = 1:Nfibtot
    contourL_vec(1,i) = contour_lengths(1,i).contourL;
end

CL = contourL_vec;
figure('Color',[1 1 1], 'Name', 'Chains with initial tangents aligned' );
checkboxrand = handles.checkbox_randomfib ;

% random index generator for selecting different fibrils without
% selecting the same twice
R1 = [1:Nfibtot];
R2 = randperm(length(R1));
for i= 1:length(R1)
    R3(i)=R1(R2(i));
end

i = 1;
while i <= n
    
    if checkboxrand ==1;
        j = R3(i);
    else
        j = i;
    end
    %j = randi(n) ;
    if (CL(j) < threshold_max)
        if (CL(j) > threshold_min)
            new=[spline_struc(1,j).bspline(1,:); spline_struc(1,j).bspline(2,:)];
            %translation of all the fibrils with first point brought to origin
            new1= [new(1,:)-new(1,1); new(2,:)-new(2,1)];
            %defines angle between first vector of the spline and x axis
            a = [new1(1,Nth_knot),new1(2,Nth_knot)];
            b = [1,0];
            costheta = dot(a,b)/(norm(a)*norm(b));
            theta = acos(costheta);
            %rotate the spline according to previous angle
            if new1(2,2) < 0
                M = makehgtform('zrotate',theta);
            else
                M = makehgtform('zrotate',-theta);
            end
            N = [M(1,1) M(1,2) ; M(2,1) M(2,2)];
            new2 = N * new1;
            %new2 is the updated coordinates after translation + rotation of new1
            plot(new2(1,:),new2(2,:),'LineWidth',linewidth)
            xlabel('nm'); %Write label for x-axis
            ylabel('nm'); %Write label for y-axis
            axis equal
            minx(i) = min(new2(1,:));
            miny(i) = min(new2(2,:));
            maxx(i) = max(new2(1,:));
            maxy(i) = max(new2(2,:));
        end
    end
    hold on
    i = i + 1 ;
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

hold off

end

guidata(hObject,handles);



% --- Executes on button press in btn_plotfibrilsraw.
function btn_plotfibrilsraw_Callback(hObject, eventdata, handles)
% hObject    handle to btn_plotfibrilsraw (see GCBO)
plotirawdatafit(hObject,handles)

function plotirawdatafit(hObject,handles)

spline_struc = handles.struc.bsplines_norm;
%Nth_knot = str2double(get(handles.txt_nthknot,'String'))+1;
Nth_knot = 1+1;
threshold_min = str2double(get(handles.txt_threshold_min,'String'));
threshold_max = str2double(get(handles.txt_threshold_max,'String'));
if threshold_max == 0
    threshold_max = 1e10;
else
end

linewidth = str2double(get(handles.txt_linewidth,'String'));
%Nfibtot = handles.ntot_fibrils;
Nfibtot = str2double(get(handles.ans_nb_fibrils,'String'));
n = str2double(get(handles.Nfib_ploti,'String'));

contour_lengths = handles.struc.contour_lengths;
for i = 1:Nfibtot
    contourL_vec(1,i) = contour_lengths(1,i).contourL;
end

CL = contourL_vec;
figure('Color',[1 1 1], 'Name', 'Chains with initial tangents aligned' );
checkboxrand = handles.checkbox_randomfib ;

% random index generator for selecting different fibrils without
% selecting the same twice
R1 = [1:Nfibtot];
R2 = randperm(length(R1));
for i= 1:length(R1)
    R3(i)=R1(R2(i));
end

i = 1;
while i <= n
    
    if checkboxrand ==1;
        j = R3(i);
    else
        j = i;
    end
    %j = randi(n) ;
    if (CL(j) < threshold_max)
        if (CL(j) > threshold_min)
            new=[spline_struc(1,j).bspline(1,:); spline_struc(1,j).bspline(2,:)];
            %translation of all the fibrils with first point brought to origin
            new1= [new(1,:)-new(1,1); new(2,:)-new(2,1)];
            %defines angle between first vector of the spline and x axis
            a = [new1(1,Nth_knot),new1(2,Nth_knot)];
            b = [1,0];
            costheta = dot(a,b)/(norm(a)*norm(b));
            theta = acos(costheta);
            %rotate the spline according to previous angle
            if new1(2,2) < 0
                M = makehgtform('zrotate',theta);
            else
                M = makehgtform('zrotate',-theta);
            end
            N = [M(1,1) M(1,2) ; M(2,1) M(2,2)];
            new2 = N * new1;
            %new2 is the updated coordinates after translation + rotation of new1
            plot(new2(1,:),new2(2,:),'LineWidth',linewidth)
            xlabel('nm'); %Write label for x-axis
            ylabel('nm'); %Write label for y-axis
            axis equal
            minx(i) = min(new2(1,:));
            miny(i) = min(new2(2,:));
            maxx(i) = max(new2(1,:));
            maxy(i) = max(new2(2,:));
        end
    end
    hold on
    i = i + 1 ;
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

hold off
    
guidata(hObject,handles);



% --- Executes on button press in checkboxworm.
function checkboxworm_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of checkboxworm
if (get(hObject,'Value') == get(hObject,'Max'))
   % Checkbox is checked-take appropriate action
   handles.checkboxworm = 1;
else
   % Checkbox is not checked-take appropriate action
   handles.checkboxworm = 0;
end
guidata(hObject,handles);

% --- Executes on button press in checkboxcorel.
function checkboxcorel_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of checkboxcorel
if (get(hObject,'Value') == get(hObject,'Max'))
   % Checkbox is checked-take appropriate action
   handles.checkboxcorel = 1;
else
   % Checkbox is not checked-take appropriate action
   handles.checkboxcorel = 0;
end
guidata(hObject,handles);

% --- Executes on button press in checkboxdevia.
function checkboxdevia_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of checkboxdevia
if (get(hObject,'Value') == get(hObject,'Max'))
   % Checkbox is checked-take appropriate action
   handles.checkboxdevia = 1;
else
   % Checkbox is not checked-take appropriate action
   handles.checkboxdevia = 0;
end
guidata(hObject,handles);




function ans_nfibrilsworm_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of ans_nfibrilsworm as text
%        str2double(get(hObject,'String')) returns contents of ans_nfibrilsworm as a double
input_worm = str2num(get(hObject,'String'));

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function ans_nfibrilsworm_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in update_worm.
function update_worm_Callback(hObject, eventdata, handles)
% hObject    handle to update_worm (see GCBO)

try
    ntot_fibrils = str2double(get(handles.ans_nfibrilsworm,'String'));
    handles.ntot_fibrils = ntot_fibrils;
    set(handles.Nfib_ploti,'String',ntot_fibrils);
    set(handles.ans_nfibrilscorel,'String',ntot_fibrils);
    set(handles.ans_nfibrilsdevia,'String',ntot_fibrils);
    convcont(hObject,handles);
    
    % update_nb_of_fib(hObject,handles);
    rgm_worm = str2double(get(handles.set_NB_of_RGM_worm,'String'));
    handles.rgm_worm = rgm_worm;
    bins_worm = str2double(get(handles.set_nb_bins_worm,'String'));
    handles.bins_worm = bins_worm;
    upmax_worm = str2double(get(handles.set_nb_upmax_worm,'String'));
    handles.upmax_worm = upmax_worm;
    lowfit_worm = str2double(get(handles.set_nb_lowfit_worm,'String'));
    handles.lowfit_worm = lowfit_worm;
    upfit_worm = str2double(get(handles.set_nb_upfit_worm,'String'));
    handles.upfit_worm = upfit_worm;
    figpop_worm = str2double(get(handles.set_nb_figpop_worm,'String'));
    handles.figpop_worm = figpop_worm;
    
    launchworm(hObject,handles);
    
catch err
    ;
end


guidata(hObject,handles);


% --- Executes on button press in update_corel.
function update_corel_Callback(hObject, eventdata, handles)
% handles    structure with handles and user data (see GUIDATA)

try
    ntot_fibrils = str2double(get(handles.ans_nfibrilscorel,'String'));
    handles.ntot_fibrils = ntot_fibrils;
    set(handles.Nfib_ploti,'String',ntot_fibrils);
    set(handles.ans_nfibrilsworm,'String',ntot_fibrils);
    set(handles.ans_nfibrilsdevia,'String',ntot_fibrils);
    convcont(hObject,handles);
    
    % update_nb_of_fib(hObject,handles);
    rgm_corel = str2double(get(handles.set_NB_of_RGM_corel,'String'));
    handles.rgm_corel = rgm_corel;
    bins_corel = str2double(get(handles.set_nb_bins_corel,'String'));
    handles.bins_corel = bins_corel;
    upmax_corel = str2double(get(handles.set_nb_upmax_corel,'String'));
    handles.upmax_corel = upmax_corel;
    lowfit_corel = str2double(get(handles.set_nb_lowfit_corel,'String'));
    handles.lowfit_corel = lowfit_corel;
    upfit_corel = str2double(get(handles.set_nb_upfit_corel,'String'));
    handles.upfit_corel = upfit_corel;
    figpop_corel = str2double(get(handles.set_nb_figpop_corel,'String'));
    handles.figpop_corel = figpop_corel;
    
    launchcorel(hObject,handles);
     
catch err
    ;
end

set(handles.uipanel_surfparam,'Visible','On') 

guidata(hObject,handles);

% --- Executes on button press in update_devia.
function update_devia_Callback(hObject, eventdata, handles)
% handles    structure with handles and user data (see GUIDATA)

try
    ntot_fibrils = str2double(get(handles.ans_nfibrilsdevia,'String'));
    handles.ntot_fibrils = ntot_fibrils;
    set(handles.Nfib_ploti,'String',ntot_fibrils);
    set(handles.ans_nfibrilsworm,'String',ntot_fibrils);
    set(handles.ans_nfibrilscorel,'String',ntot_fibrils);
    convcont(hObject,handles);
    
    % update_nb_of_fib(hObject,handles);
    rgm_devia = str2double(get(handles.set_NB_of_RGM_devia,'String'));
    handles.rgm_devia = rgm_devia;
    bins_devia = str2double(get(handles.set_nb_bins_devia,'String'));
    handles.bins_devia = bins_devia;
    upmax_devia = str2double(get(handles.set_nb_upmax_devia,'String'));
    handles.upmax_devia = upmax_devia;
    lowfit_devia = str2double(get(handles.set_nb_lowfit_devia,'String'));
    handles.lowfit_devia = lowfit_devia;
    upfit_devia = str2double(get(handles.set_nb_upfit_devia,'String'));
    handles.upfit_devia = upfit_devia;
    figpop_devia = str2double(get(handles.set_nb_figpop_devia,'String'));
    handles.figpop_devia = figpop_devia;
    
    launchdevia(hObject,handles);
    
catch err
    ;
end


guidata(hObject,handles);



% --- Generates plot and fit from random matrices,  and return P-length.
function launchworm(hObject,handles)


n = str2double(get(handles.ans_nfibrilsworm,'String'));
z = str2double(get(handles.ans_nb_fibrils,'String'));
structwo = handles.struc.mat_wormlike;
rgm = handles.rgm_worm;
numBins = handles.bins_worm;
upmax = handles.upmax_worm;
upfit = handles.upfit_worm;
lowfit = handles.lowfit_worm;
figpop = handles.figpop_worm;
checkboxworm = handles.checkboxworm;

for j = 1:rgm
    new1 = [];
    i = 1;
% THAT code is ok although a bit different: it selects only half of the 
% fibrils and make sure it does not pick the same fibril twice for each
% iteration.  
%     while i <= n/2
%         %R1 = []; R2 = []; R3 = [];
%         R1 = [1:z];
%         R2 = randperm(length(R1));
%         for pos= 1:length(R1)
%             R3(pos)=R1(R2(pos));
%         end
%         new1 = [new1; structwo(1,R3(pos)).wormfib];
%         i = i + 1;
%     end
        while i <= n
            new1 = [new1; structwo(1,randi(z)).wormfib];
            i = i + 1;
        end
            new1(any(isnan(new1),2),:)=[];
            new1 = real(new1);
        
        randworm(j).randwomat = new1;
        
    %BIN WORM DATA
    topEdge = 0;
    x = randworm(j).randwomat(:,2); %split into x and y
    y = randworm(j).randwomat(:,1);
    botEdge = round(min(x)); % define limits
    if upmax >= max(x)
        topEdge = round(max(x)); % define limits
    else
        topEdge = upmax;
    end
    binEdges = linspace(botEdge, topEdge, numBins+1);
    [~,whichBin] = histc(x, binEdges);
    for k = 1:numBins
        flagBinMembers = (whichBin == k);
        binMembers     = y(flagBinMembers);
        binMean_wo(k)  = mean((binMembers).^2);
        %---setup the length as the center of two separate bin edges
        length_wo(k) = binEdges(k) + ( ( (topEdge - botEdge) / numBins) /2 );
    end
    value = num2str(numBins);
    worm_length(j).lengthwo = length_wo;
    worm_binMean(j).binmeanwo = binMean_wo;
    
    % Fit: 'wormfit'.
    [xData, yData] = prepareCurveData( worm_length(1,j).lengthwo(1,:),...
        worm_binMean(1,j).binmeanwo(1,:) );
    
    % Set up fittype and options.
    ft = fittype( '4*p*x*(1-2*p/x*(1-exp(-x/2/p)))', 'independent', 'x',...
        'dependent', 'y' );
    opts = fitoptions( ft );
    opts.Display = 'Off';
    opts.Lower = -Inf;
    opts.StartPoint = 500;
    opts.Upper = Inf;
    
    if upfit <= topEdge
        ex = excludedata( xData, yData, 'domain', [lowfit  upfit] );
        opts.Exclude = ex;
        [fitresult, gof] = fit( xData, yData, ft, opts );
        if j <= figpop
            if checkboxworm == 1
                fi = figure('Color',[0.9 0.9 0.9], 'Name', 'Contour length versus end-to-end distance (Bootstrap)' );
                set(fi,'Position', [10 400  560  420])  %worm
                plot( fitresult, xData, yData, ex );
                xlabel( 'Contour length (nm)' );
                ylabel( 'Mean-square end-to-end distance (nm^2)' );
                set(legend,'units','normalized');
                set(legend,'Position',[0.2 0.77 0.2 0.12]);
                grid on
            end
        end
    else
        [fitresult, gof] = fit( xData, yData, ft, opts );
        if j <= figpop
            if checkboxworm == 1
                fi = figure('Color',[0.9 0.9 0.9], 'Name', 'Contour length versus end-to-end distance (Bootstrap)' );
                set(fi,'Position', [10 400  560  420])  %worm
                plot( fitresult, xData, yData );
                xlabel( 'Contour length (nm)' );
                ylabel( 'Mean-square end-to-end distance (nm^2)' );
                set(legend,'units','normalized');
                set(legend,'Position',[0.2 0.77 0.2 0.12]);
                grid on
            end
        end
    end
    c = fit( xData, yData, ft, opts );
    persisworm(j,1).fit = coeffvalues(c);
    persisworm(j,1).gof = gof;
    
    topEdgemat(1,j) = topEdge;
end

topEdge = max(topEdgemat);

for i = 1:rgm
    persis_vec(i,1) = persisworm(i,1).fit;
    rsquare_vec(i,1) = persisworm(i,1).gof.rsquare;
end

toggle23 = get(handles.toggle_2D3Dswitcher,'String');
TF = strcmp(toggle23,'Switch to Noneq.');
if TF == 1  % equilibrated in 2D
    SIPworm = 1 ;
    pworm_mean = round(mean(persis_vec)) * SIPworm;
    pworm_std = round(std(persis_vec)) * SIPworm;
else           % not equilibrated ---- propagate uncertainties
    SIPworm = 1.5 ;
    lowerB = round(mean(persis_vec));
    pworm_mean = lowerB * SIPworm;
    Err_res = round(std(persis_vec));  %Err_res is error on resampling
    Err_dim = 0.5;  %Err_res is error on fractional dimension
    pworm_std = pworm_mean .*(sqrt( (Err_res ./ lowerB)^2 + (Err_dim ./ SIPworm)^2 ));
    pworm_std = round(pworm_std);
end


rsqworm_mean = round2(mean(rsquare_vec),0.001);
rsqworm_std = round2(std(rsquare_vec),0.001);
zlowfit = lowfit;
zupfit = upfit;

kB = 1.3806503e-23 ;  %Boltzmann constant

Temp = str2double(get(handles.set_Temperature,'String'));
T = 273.15 + Temp; % Temperature in deg Kelvin
CBworm_mean = round2(pworm_mean*kB*T/1E9,10e-30);
CBworm_std = round2(pworm_std*kB*T/1E9,10e-30);

set(handles.ans_pworm_mean,'String',round(pworm_mean));
set(handles.ans_pworm_std,'String',round(pworm_std));
set(handles.ans_CBworm_mean,'String',CBworm_mean);
set(handles.ans_CBworm_std,'String',CBworm_std);
set(handles.ans_rsqworm_std,'String',rsqworm_std);

if rsqworm_mean <= 0
    set(handles.ans_rsqworm_mean,'String',rsqworm_mean, 'foregroundcolor',[1 0 0]);
elseif rsqworm_mean > 0 && rsqworm_mean < 0.9
    set(handles.ans_rsqworm_mean,'String',rsqworm_mean, 'foregroundcolor',[0 0 0]);
elseif rsqworm_mean >= 0.9
    set(handles.ans_rsqworm_mean,'String',rsqworm_mean, 'foregroundcolor',[0 0.5 0]);
    set(handles.uipanel_saveoutputs,'Visible','On')
end

set(handles.set_nb_upmax_worm,'String',topEdge);
if upfit < topEdge
else
    set(handles.set_nb_upfit_worm,'String',topEdge);
end



guidata(hObject,handles);



% --- Generates plot and fit from random matrices,  and return P-length.
function launchcorel(hObject,handles)


n = str2double(get(handles.ans_nfibrilscorel,'String'));
z = str2double(get(handles.ans_nb_fibrils,'String'));
structco = handles.struc.mat_correlations;
rgm = handles.rgm_corel;
numBins = handles.bins_corel;
upmax = handles.upmax_corel;
upfit = handles.upfit_corel;
lowfit = handles.lowfit_corel;
figpop = handles.figpop_corel;
checkboxcorel = handles.checkboxcorel;

for j = 1:rgm
    new2 = [];
    i = 1;
    
%     while i <= n/2
%         %R1 = []; R2 = []; R3 = [];
%         R1 = [1:z];
%         R2 = randperm(length(R1));
%         for pos= 1:length(R1)
%             R3(pos)=R1(R2(pos));
%         end
%         new2 = [new2; structco(1,R3(pos)).corelfib];
%         i = i + 1;
%     end
        while i <= n
            new2 = [new2; structco(1,randi(z)).corelfib];
            i = i + 1;
        end    
        new2(any(isnan(new2),2),:)=[];
        new2 = real(new2);
        randcorel(j).randcomat = new2;
    
    %BIN COREL DATA
    x = randcorel(j).randcomat(:,2); %split into x and y
    y = randcorel(j).randcomat(:,1);
    botEdge = round(min(x)); % define limits
    if upmax >= max(x)
        topEdge = round(max(x)); % define limits
    else
        topEdge = upmax;
    end
    binEdges = linspace(botEdge, topEdge, numBins+1);
    [~,whichBin] = histc(x, binEdges);
    for k = 1:numBins
        flagBinMembers = (whichBin == k);
        binMembers     = y(flagBinMembers);
        binMean_co(k)  = mean((binMembers));
        %---setup the length as the center of two separate bin edges
        length_co(k) = binEdges(k) + ( ( (topEdge - botEdge) / numBins) /2 );
    end
    value = num2str(numBins);
    corel_length(j).lengthco = length_co;
    corel_binMean(j).binmeanco = binMean_co;
    
    % Fit: 'corelfit'.
    [xData, yData] = prepareCurveData( corel_length(1,j).lengthco(1,:),...
        corel_binMean(1,j).binmeanco(1,:) );
    
    % Set up fittype and options.
    ft = fittype( 'exp(-x/2/p)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( ft );
    opts.Display = 'Off';
    opts.Lower = -Inf;
    opts.StartPoint = 500;
    opts.Upper = Inf;
    
    if upfit <= topEdge
        ex = excludedata( xData, yData, 'domain', [lowfit  upfit] );
        opts.Exclude = ex;
        [fitresult, gof] = fit( xData, yData, ft, opts );
        if j <= figpop
            if checkboxcorel == 1
                fi = figure('Color',[0.9 0.9 0.9], 'Name', 'Tangent-tangent correlations (Bootstrap) ' );
                set(fi, 'Position', [580 400  560  420]) %corel
                plot( fitresult, xData, yData, ex );
                xlabel( 'Contour length (nm)' );
                ylabel( 'Mean tangent correlations' );
                grid on
            end
        end
    else
        [fitresult, gof] = fit( xData, yData, ft, opts );
        if j <= figpop
            if checkboxcorel == 1
                fi = figure('Color',[0.9 0.9 0.9], 'Name', 'Tangent-tangent correlations (Bootstrap) ' );
                set(fi, 'Position', [580 400  560  420]) %corel
                plot( fitresult, xData, yData );
                xlabel( 'Contour length (nm)' );
                ylabel( 'Mean tangent correlations' );
                grid on
            end
        end
    end
    c = fit( xData, yData, ft, opts );
    persiscorel(j,1).fit = coeffvalues(c);
    persiscorel(j,1).gof = gof;
    
    topEdgemat(1,j) = topEdge;
end

topEdge = max(topEdgemat);

for i = 1:rgm
    persis_vec(i,1) = persiscorel(i,1).fit;
    rsquare_vec(i,1) = persiscorel(i,1).gof.rsquare;
end

toggle23 = get(handles.toggle_2D3Dswitcher,'String');
TF = strcmp(toggle23,'Switch to Noneq.');
if TF == 1  % equilibrated in 2D
    SIPcorel = 1 ;
    pcorel_mean = round(mean(persis_vec)) * SIPcorel;
    pcorel_std = round(std(persis_vec)) * SIPcorel;
else           % not equilibrated ---- propagate uncertainties
    SIPcorel = 1.5 ;
    lowerB = round(mean(persis_vec));
    pcorel_mean = lowerB * SIPcorel;
    Err_res = round(std(persis_vec));  %Err_res is error on resampling
    Err_dim = 0.5;  %Err_res is error on fractional dimension
    pcorel_std = pcorel_mean .*(sqrt( (Err_res ./ lowerB)^2 + (Err_dim ./ SIPcorel)^2 ));
    pcorel_std = round(pcorel_std);
end

rsqcorel_mean = round2(mean(rsquare_vec),0.001);
rsqcorel_std = round2(std(rsquare_vec),0.001);
zlowfit = lowfit;
zupfit = upfit;

kB = 1.3806503e-23 ;  %Boltzmann constant
Temp = str2double(get(handles.set_Temperature,'String'));
T = 273.15 + Temp; % Temperature in deg Kelvin
CBcorel_mean = round2(pcorel_mean*kB*T/1E9,10e-30);
CBcorel_std = round2(pcorel_std*kB*T/1E9,10e-30);

set(handles.ans_pcorel_mean,'String',round(pcorel_mean));
set(handles.ans_pcorel_std,'String',round(pcorel_std));
set(handles.ans_CBcorel_mean,'String',CBcorel_mean);
set(handles.ans_CBcorel_std,'String',CBcorel_std);
set(handles.ans_rsqcorel_std,'String',rsqcorel_std);

if rsqcorel_mean <= 0
    set(handles.ans_rsqcorel_mean,'String',rsqcorel_mean, 'foregroundcolor',[1 0 0]);
elseif rsqcorel_mean > 0 && rsqcorel_mean < 0.9
    set(handles.ans_rsqcorel_mean,'String',rsqcorel_mean, 'foregroundcolor',[0 0 0]);
elseif rsqcorel_mean >= 0.9
    set(handles.ans_rsqcorel_mean,'String',rsqcorel_mean, 'foregroundcolor',[0 0.5 0]);
    set(handles.uipanel_saveoutputs,'Visible','On')
end

set(handles.set_nb_upmax_corel,'String',topEdge);
if upfit < topEdge
else
    set(handles.set_nb_upfit_corel,'String',topEdge);
end


guidata(hObject,handles);



% --- Generates plot and fit from random matrices,  and return P-length.
function launchdevia(hObject,handles)


n = str2double(get(handles.ans_nfibrilsdevia,'String'));
z = str2double(get(handles.ans_nb_fibrils,'String'));
structdev = handles.struc.mat_deviations;
rgm = handles.rgm_devia;
numBins = handles.bins_devia;
upmax = handles.upmax_devia;
upfit = handles.upfit_devia;
lowfit = handles.lowfit_devia;
figpop = handles.figpop_devia;
checkboxdevia = handles.checkboxdevia;

for j = 1:rgm
    new3 = [];
    i = 1;
    
%     while i <= n/2
%         %R1 = []; R2 = []; R3 = [];
%         R1 = [1:z];
%         R2 = randperm(length(R1));
%         for pos= 1:length(R1)
%             R3(pos)=R1(R2(pos));
%         end
%         new3 = [new3; structdev(1,R3(pos)).deviafib];
%         i = i + 1;
%     end
        while i <= n
            new3 = [new3; structdev(1,randi(z)).deviafib];
            i = i + 1;
        end
        new3(any(isnan(new3),2),:)=[]; %% Remove all rows that contain
        % at least one NaN in matrix
        new3 = real(new3); % Keep only real part of complex numbers
        randdevia(j).randdevmat = new3;
    
    %BIN DEVIA DATA
    x = randdevia(j).randdevmat(:,2); %split into x and y
    y = randdevia(j).randdevmat(:,1);
    botEdge = round(min(x)); % define limits
    if upmax >= max(x)
        topEdge = round(max(x)); % define limits
    else
        topEdge = upmax;
    end
    binEdges = linspace(botEdge, topEdge, numBins+1);
    [~,whichBin] = histc(x, binEdges);
    for k = 1:numBins
        flagBinMembers = (whichBin == k);
        binMembers     = y(flagBinMembers);
        binMean_dev(k)  = mean((binMembers).^2);
        %---setup the length as the center of two separate bin edges
        length_dev(k) = binEdges(k) + ( ( (topEdge - botEdge) / numBins) /2 );
    end
    value = num2str(numBins);
    devia_length(j).lengthwo = length_dev;
    devia_binMean(j).binmeanwo = binMean_dev;
    
    % Fit: 'deviafit'.
    [xData, yData] = prepareCurveData( devia_length(1,j).lengthwo(1,:),...
        devia_binMean(1,j).binmeanwo(1,:) );
    
    % Set up fittype and options.
    ft = fittype( 'x^3/48/p', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( ft );
    opts.Display = 'Off';
    opts.Lower = -Inf;
    opts.StartPoint = 500;
    opts.Upper = Inf;
    
    if upfit <= topEdge
        ex = excludedata( xData, yData, 'domain', [lowfit  upfit] );
        opts.Exclude = ex;
        [fitresult, gof] = fit( xData, yData, ft, opts );
        if j <= figpop
            if checkboxdevia == 1
                fi = figure('Color',[0.9 0.9 0.9], 'Name', 'Deviations from secant midpoint (Bootstrap) ' );
                set(fi, 'Position', [850 400  560  420]) %devia
                plot( fitresult, xData, yData, ex );
                xlabel( 'Secant length (nm)' );
                ylabel( 'Mean-square deviations (nm^2)' );
                set(legend,'units','normalized');
                set(legend,'Position',[0.2 0.77 0.2 0.12]);
                grid on
            end
        end
    else
        [fitresult, gof] = fit( xData, yData, ft, opts );
        if j <= figpop
            if checkboxdevia == 1
                fi = figure('Color',[0.9 0.9 0.9], 'Name', 'Deviations from secant midpoint (Bootstrap)' );
                set(fi, 'Position', [850 400  560  420]) %devia
                plot( fitresult, xData, yData );
                xlabel( 'Secant length (nm)' );
                ylabel( 'Mean-square deviations (nm^2)' );
                set(legend,'units','normalized');
                set(legend,'Position',[0.2 0.77 0.2 0.12]);
                grid on
            end
        end
    end
    c = fit( xData, yData, ft, opts );
    persisdevia(j,1).fit = coeffvalues(c);
    persisdevia(j,1).gof = gof;
    
    topEdgemat(1,j) = topEdge;
end

topEdge = max(topEdgemat);

for i = 1:rgm
    persis_vec(i,1) = persisdevia(i,1).fit;
    rsquare_vec(i,1) = persisdevia(i,1).gof.rsquare;
end

toggle23 = get(handles.toggle_2D3Dswitcher,'String');
TF = strcmp(toggle23,'Switch to Noneq.');
if TF == 1  % equilibrated in 2D
    SIPdevia = 1 ;
    pdevia_mean = round(mean(persis_vec)) * SIPdevia;
    pdevia_std = round(std(persis_vec)) * SIPdevia;
else           % not equilibrated ---- propagate uncertainties
    SIPdevia = 1.5 ;
    lowerB = round(mean(persis_vec));
    pdevia_mean = lowerB * SIPdevia;
    Err_res = round(std(persis_vec));  %Err_res is error on resampling
    Err_dim = 0.5;  %Err_res is error on fractional dimension
    pdevia_std = pdevia_mean .*(sqrt( (Err_res ./ lowerB)^2 + (Err_dim ./ SIPdevia)^2 ));
    pdevia_std = round(pdevia_std);
end

rsqdevia_mean = round2(mean(rsquare_vec),0.001);
rsqdevia_std = round2(std(rsquare_vec),0.001);
zlowfit = lowfit;
zupfit = upfit;

kB = 1.3806503e-23 ;  %Boltzmann constant
Temp = str2double(get(handles.set_Temperature,'String'));
T = 273.15 + Temp; % Temperature in deg Kelvin
CBdevia_mean = round2(pdevia_mean*kB*T/1E9,10e-30);
CBdevia_std = round2(pdevia_std*kB*T/1E9,10e-30);

set(handles.ans_pdevia_mean,'String',round(pdevia_mean));
set(handles.ans_pdevia_std,'String',round(pdevia_std));
set(handles.ans_CBdevia_mean,'String',CBdevia_mean);
set(handles.ans_CBdevia_std,'String',CBdevia_std);
set(handles.ans_rsqdevia_std,'String',rsqdevia_std);

if rsqdevia_mean <= 0
    set(handles.ans_rsqdevia_mean,'String',rsqdevia_mean, 'foregroundcolor',[1 0 0]);
elseif rsqdevia_mean > 0 && rsqdevia_mean < 0.9
    set(handles.ans_rsqdevia_mean,'String',rsqdevia_mean, 'foregroundcolor',[0 0 0]);
elseif rsqdevia_mean >= 0.9
    set(handles.ans_rsqdevia_mean,'String',rsqdevia_mean, 'foregroundcolor',[0 0.5 0]);
    set(handles.uipanel_saveoutputs,'Visible','On')
end


set(handles.set_nb_upmax_devia,'String',topEdge);
if upfit < topEdge
else
    set(handles.set_nb_upfit_devia,'String',topEdge);
end



guidata(hObject,handles);


function set_nb_figpop_devia_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of set_nb_figpop_devia as text
%        str2double(get(hObject,'String')) returns contents of set_nb_figpop_devia as a double


% --- Executes during object creation, after setting all properties.
function set_nb_figpop_devia_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_NB_of_RGM_devia_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of set_NB_of_RGM_devia as text
%        str2double(get(hObject,'String')) returns contents of set_NB_of_RGM_devia as a double


% --- Executes during object creation, after setting all properties.
function set_NB_of_RGM_devia_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_NB_of_RGM_corel_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of set_NB_of_RGM_corel as text
%        str2double(get(hObject,'String')) returns contents of set_NB_of_RGM_corel as a double


% --- Executes during object creation, after setting all properties.
function set_NB_of_RGM_corel_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_nb_figpop_corel_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of set_nb_figpop_corel as text
%        str2double(get(hObject,'String')) returns contents of set_nb_figpop_corel as a double


% --- Executes during object creation, after setting all properties.
function set_nb_figpop_corel_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end







function set_nb_figpop_worm_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of set_nb_figpop_worm as text
%        str2double(get(hObject,'String')) returns contents of set_nb_figpop_worm as a double


% --- Executes during object creation, after setting all properties.
function set_nb_figpop_worm_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_NB_of_RGM_worm_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of set_NB_of_RGM_worm as text
%        str2double(get(hObject,'String')) returns contents of set_NB_of_RGM_worm as a double


% --- Executes during object creation, after setting all properties.
function set_NB_of_RGM_worm_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_nb_bins_worm_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of set_nb_bins_worm as text
%        str2double(get(hObject,'String')) returns contents of set_nb_bins_worm as a double


% --- Executes during object creation, after setting all properties.
function set_nb_bins_worm_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_nb_upmax_worm_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of set_nb_upmax_worm as text
%        str2double(get(hObject,'String')) returns contents of set_nb_upmax_worm as a double


% --- Executes during object creation, after setting all properties.
function set_nb_upmax_worm_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_nb_bins_corel_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of set_nb_bins_corel as text
%        str2double(get(hObject,'String')) returns contents of set_nb_bins_corel as a double


% --- Executes during object creation, after setting all properties.
function set_nb_bins_corel_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_nb_upmax_corel_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of set_nb_upmax_corel as text
%        str2double(get(hObject,'String')) returns contents of set_nb_upmax_corel as a double


% --- Executes during object creation, after setting all properties.
function set_nb_upmax_corel_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_nb_bins_devia_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of set_nb_bins_devia as text
%        str2double(get(hObject,'String')) returns contents of set_nb_bins_devia as a double


% --- Executes during object creation, after setting all properties.
function set_nb_bins_devia_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_nb_upmax_devia_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of set_nb_upmax_devia as text
%        str2double(get(hObject,'String')) returns contents of set_nb_upmax_devia as a double


% --- Executes during object creation, after setting all properties.
function set_nb_upmax_devia_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_nb_lowfit_devia_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of set_nb_lowfit_devia as text
%        str2double(get(hObject,'String')) returns contents of set_nb_lowfit_devia as a double


% --- Executes during object creation, after setting all properties.
function set_nb_lowfit_devia_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_nb_upfit_devia_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of set_nb_upfit_devia as text
%        str2double(get(hObject,'String')) returns contents of set_nb_upfit_devia as a double


% --- Executes during object creation, after setting all properties.
function set_nb_upfit_devia_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_nb_lowfit_corel_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of set_nb_lowfit_corel as text
%        str2double(get(hObject,'String')) returns contents of set_nb_lowfit_corel as a double


% --- Executes during object creation, after setting all properties.
function set_nb_lowfit_corel_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_nb_upfit_corel_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of set_nb_upfit_corel as text
%        str2double(get(hObject,'String')) returns contents of set_nb_upfit_corel as a double


% --- Executes during object creation, after setting all properties.
function set_nb_upfit_corel_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_nb_lowfit_worm_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of set_nb_lowfit_worm as text
%        str2double(get(hObject,'String')) returns contents of set_nb_lowfit_worm as a double


% --- Executes during object creation, after setting all properties.
function set_nb_lowfit_worm_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_nb_upfit_worm_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of set_nb_upfit_worm as text
%        str2double(get(hObject,'String')) returns contents of set_nb_upfit_worm as a double


% --- Executes during object creation, after setting all properties.
function set_nb_upfit_worm_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function bins_test2D_Callback(hObject, eventdata, handles)
% hObject    handle to bins_test2D (see GCBO)
% 

% --- Executes during object creation, after setting all properties.
function bins_test2D_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pushbut_test2D.
function pushbut_test2D_Callback(hObject, eventdata, handles)
% hObject    handle to pushbut_test2D (see GCBO)
test2Dequilibration(hObject, handles);



% ---calculates and plots the ratio of 4th/2nd moments of probability
% ---distribution function in 2D according to Rivetti, Bustamente et al., JMB 1996 
function test2Dequilibration(hObject, handles)

PLcorel_mean = str2double(get(handles.ans_pcorel_mean,'String'));

if PLcorel_mean == 0
    msgbox(sprintf('Before pressing this button, you must collect some data regarding the distribution of angles in your sample. This can be accomplished by pressing the "Launch Fit" button in the "Tangent Correlations" panel.'),...
        'Attention','Warn')
    return
else 
end

z = handles.ntot_fibrils;
mod1 = str2double(get(handles.ans_nfibrilsworm,'String'));
mod2 = str2double(get(handles.ans_nfibrilscorel,'String'));
mod3 = str2double(get(handles.ans_nfibrilsdevia,'String'));

if mod1 < z
    z = mod1;
elseif mod2 < z
    z = mod2;
elseif mod3 < z
    z = mod3;
end

numBins_co = str2double(get(handles.bins_test2D,'String'));
maxval_co1 = str2double(get(handles.set_nb_upfit_corel,'String'));
maxval_co2 = str2double(get(handles.set_nb_upmax_corel,'String')); 
if maxval_co1 < maxval_co2
    maxval_co  = maxval_co1;
else maxval_co  = maxval_co2;
    set(handles.set_nb_upfit_corel, 'String',maxval_co);
end
minval_co = str2double(get(handles.set_nb_lowfit_corel,'String'));

    matcor = handles.struc.mat_correlations;
    new1 = [];
     for i = 1:z
        new1 = [new1; matcor(1,i).corelfib];
     end
     vectan = new1;  % vectan(,2) is the contour length of each fib
 vectan(:,3) = acos(vectan(:,1)); %switch back to radians (from cos values)


%BIN COREL DATA
    a = vectan(:,2); %split into x and y
    b = vectan(:,3);

%This line is dedicated to storing all segments generated.  THe value of 
%each segment is its length and this can be binned later, e.g. to be displayed 
%in a histogram, together with kurtosis value.
handles.kurtosisSegL = a;
%%%---------------------------------------------------------------------    
     
    botEdge = minval_co; % define limits
    topEdge = maxval_co; % define limits
    
    binEdges = linspace(botEdge, topEdge, numBins_co+1);
    [~,whichBin] = histc(a, binEdges);
    for i = 1:numBins_co
        flagBinMembers = (whichBin == i);
        binMembers     = b(flagBinMembers);
        %---MODIFIED EXPRESSION HERE
         binMean_tan(i) = (mean(binMembers.^4))/((mean(binMembers.^2)).^2);
        %---setup the length as the center of two separate bin edges
        length_tan(i) = binEdges(i) + (((topEdge - botEdge)/numBins_co)/2 );
    end
    value = num2str(numBins_co);
    test2D_length = length_tan;
    test2D_binMean = binMean_tan;

    f = figure('Color',[1 1 1], 'Name', 'Test 2D - Kurtosis' );
    set(f, 'Position', [100 100 800 250])
    scatter(test2D_length,test2D_binMean) 
    hold on
    x = botEdge:.1:topEdge;
    y = 3;
    p = plot(x,y);
    set(p,'Color','red','LineWidth',2);
    
    xmax = topEdge;
    axis([0 xmax 0 18])
    grid on
    xlabel('Contour length (nm)')
    ylabel('Kurtosis of \theta distribution')
    
%Note1 = ['   Note:'];
Note2 = ['   Note: For contour length < persistence length,   '];
Note3 = ['   if full equilibration in 2D occurs, then we have: Kurtosis = 3   '];
Note4 = ['   (Set boundaries in Outliers of Tangent Correlations panel)   '];
            
        mTextBox = text;
        set(mTextBox,'units','normalized');
        set(mTextBox,'Position',[0.47 0.85 0]);
        set(mTextBox,'HorizontalAlignment', 'left');
        set(mTextBox,'fontsize',9);
        set(mTextBox,'color',[0 0 0 ]);
        set(mTextBox,'edgecolor',[0.5 0.5 0.5]);
        set(mTextBox,'String',{ Note2, Note3, Note4});
    set(mTextBox,'BackgroundColor',[1 1 1]);
        
    handles.test2D_length = test2D_length;
    handles.test2D_binmean = test2D_binMean;
    
guidata(hObject, handles);



% --- Executes on button press in btn_sumall.
function btn_sumall_Callback(hObject, eventdata, handles)
% hObject    handle to btn_sumall (see GCBO)
sumandbin(hObject,handles);

%--- Concatenate matrices generated by guicalc and bin the data
function sumandbin(hObject,handles)


try
    
    PLworm_mean = str2double(get(handles.ans_pworm_mean,'String'));
    PLcorel_mean = str2double(get(handles.ans_pcorel_mean,'String'));
    PLdevia_mean = str2double(get(handles.ans_pdevia_mean,'String'));
    
    if PLworm_mean == 0 || PLcorel_mean == 0 || PLdevia_mean == 0
        msgbox(sprintf('Please initialize all fitting procedures by pressing all Launch Fit buttons in all of the following panels: Contour / End-to-end, Tangent Correlations, and Deviations / Secant Midpoint.'),...
            'Attention','Warn')
        return
    else
    end
    
    z = handles.ntot_fibrils;
    
    mod1 = str2double(get(handles.ans_nfibrilsworm,'String'));
    mod2 = str2double(get(handles.ans_nfibrilscorel,'String'));
    mod3 = str2double(get(handles.ans_nfibrilsdevia,'String'));
    
    if mod1 < z
        z = mod1;
    elseif mod2 < z
        z = mod2;
    elseif mod3 < z
        z = mod3;
    end
    
    set(handles.ans_nfibrilsworm,'String',z) ;
    set(handles.ans_nfibrilscorel,'String',z) ;
    set(handles.ans_nfibrilsdevia,'String',z) ;
    
    matwo = handles.struc.mat_wormlike;
    matcor = handles.struc.mat_correlations;
    matdev = handles.struc.mat_deviations;
    numBins_wo = str2double(get(handles.set_nb_bins_worm,'String'));
    numBins_co = str2double(get(handles.set_nb_bins_corel,'String'));
    numBins_dev = str2double(get(handles.set_nb_bins_devia,'String'));
    maxval_wo = str2double(get(handles.set_nb_upfit_worm,'String'));
    maxval_dev = str2double(get(handles.set_nb_upfit_devia,'String'));
    maxval_co = str2double(get(handles.set_nb_upfit_corel,'String'));
    
    new1 = [];
    for i = 1:z
        new1 = [new1; matwo(1,i).wormfib];
    end
    new1(any(isnan(new1),2),:)=[];
    new1 = real(new1);
    vecworm = new1;
    
    %BIN WORM DATA
    x = vecworm(:,2); %split into x and y
    y = vecworm(:,1);
    botEdge = round(min(x)); % define limits
    topEdge = maxval_wo; % define limits
    %numBins = 100; % define number of bins
    binEdges = linspace(botEdge, topEdge, numBins_wo+1);
    [~,whichBin] = histc(x, binEdges);
    for i = 1:numBins_wo
        flagBinMembers = (whichBin == i);
        binMembers     = y(flagBinMembers);
        binMean_wo(i)  = mean((binMembers).^2);
        %---setup the length as the center of two separate bin edges
        length_wo(i) = binEdges(i) + ( ( (topEdge - botEdge) / numBins_wo) /2 );
    end
    
    worm_length = length_wo;
    worm_binMean = binMean_wo;
    %     value = num2str(numBins_wo);
    %     filename=[value,'bins_worm'];
    
    
    
    %---------------------------------------------------
    %MAKE FIGURE for the contour/end2end fit
    
    checkboxwo = handles.checkboxworm;
    
    %if checkboxwo == 1
    % Fit: 'deviafit'.
    [xData_wo, yData_wo] = prepareCurveData( worm_length, worm_binMean);
    
    % Set up fittype and options.
    ft_wo = fittype( '4*p*x*(1-2*p/x*(1-exp(-x/2/p)))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( ft_wo );
    opts.Display = 'Off';
    opts.Lower = -Inf;
    opts.StartPoint = 500;
    opts.Upper = Inf;
    [fitresult_wo, gof_wo] = fit( xData_wo, yData_wo, ft_wo, opts );
    
    h1  = figure('Color',[1 1 1], 'Name', 'Contour length versus end-to-end distance (all chains)' );
    set(h1, 'Position', [10 400  560  420])  %worm
    scatter( xData_wo, yData_wo );
    hold on
    plot(fitresult_wo);
    legend('off')
    xlabel( 'Contour length (nm)');
    ylabel( 'Mean square end-to-end distance (nm^2)');
    grid off
    
    PLworm_mean = get(handles.ans_pworm_mean,'String');
    PLworm_std = get(handles.ans_pworm_std,'String');
    
    toggle23 = get(handles.toggle_2D3Dswitcher,'String');
    TF = strcmp(toggle23,'Switch to Noneq.');
    if TF == 1
        SIParam = '2D fit';
    else
        SIParam = 'Noneq. fit';
    end
    PLworm = ['Persistence Length   =   ', PLworm_mean, '  ±  ', PLworm_std, '  nm', '    (', SIParam, ')'];
    BRworm_mean = get(handles.ans_CBworm_mean,'String');
    BRworm_std = get(handles.ans_CBworm_std,'String');
    BRworm = ['Bending Rigidity        =   ', BRworm_mean, '  ±  ', BRworm_std, '  N.m^2'];
    coeffD_worm = get(handles.ans_rsqworm_mean,'String');
    coeffD_worm = ['Determination Coeff   =   ',coeffD_worm];
    
    % mTextBox = uicontrol('style','text');
    %set(mTextBox,'Position',[120 390 390 60]);
    mTextBox_wo = text;
    set(mTextBox_wo,'units','normalized');
    set(mTextBox_wo,'Position',[0.05 0.9 0]);
    set(mTextBox_wo,'HorizontalAlignment', 'left');
    set(mTextBox_wo,'fontsize',9);
    set(mTextBox_wo,'color',[0 0 0 ]);
    set(mTextBox_wo,'edgecolor',[1 1 1]);
    set(mTextBox_wo,'String',{PLworm, BRworm, coeffD_worm });
    %set(mTextBox,'BackgroundColor',[1 1 1]);
    %else
    %end
    
    
    
    % COREL stuff starts here
    new3 = [];
    for i = 1:z
        new3 = [new3; matcor(1,i).corelfib];
    end
    new3(any(isnan(new3),2),:)=[];
    new3 = real(new3);
    vectan = new3;
    
    
    %BIN COREL DATA
    a = vectan(:,2); %split into x and y
    b = vectan(:,1);
    botEdge = 0; % define limits
    topEdge = maxval_co; % define limits
    %numBins = 100; % define number of bins
    binEdges = linspace(botEdge, topEdge, numBins_co+1);
    [~,whichBin] = histc(a, binEdges);
    for i = 1:numBins_co
        flagBinMembers = (whichBin == i);
        binMembers     = b(flagBinMembers);
        binMean_tan(i) = mean(binMembers);  % IT IS NORMAL we don't put the square here
        %---setup the length as the center of two separate bin edges
        length_tan(i) = binEdges(i) + ( ( (topEdge - botEdge) / numBins_co) /2 );
    end
    
    corel_length = length_tan;
    corel_binMean = binMean_tan;
    %     value = num2str(numBins_co);
    %     filename=[value,'bins_corel'];
    
    
    
    %---------------------------------------------------
    %MAKE FIGURE for the corel fit
    
    checkboxcor = handles.checkboxcorel;
    
    %if checkboxcor == 1
    % Fit: 'deviafit'.
    [xData_co, yData_co] = prepareCurveData( corel_length, corel_binMean);
    
    % Set up fittype and options.
    ft_co = fittype( 'exp(-x/2/p)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( ft_co );
    opts.Display = 'Off';
    opts.Lower = -Inf;
    opts.StartPoint = 500;
    opts.Upper = Inf;
    [fitresult_co, gof_co] = fit( xData_co, yData_co, ft_co, opts );
    
    h2 = figure('Color',[1 1 1], 'Name', 'Tangent-tangent correlations (all chains)' );
    set(h2, 'Position', [580 400  560  420])
    scatter( xData_co, yData_co );
    hold on
    plot(fitresult_co);
    legend('off')
    xlabel( 'Contour length (nm)');
    ylabel( '< cos \theta > ');
    grid off
    
    PLcorel_mean = get(handles.ans_pcorel_mean,'String');
    PLcorel_std = get(handles.ans_pcorel_std,'String');
    
    PLcorel = ['Persistence Length   =   ', PLcorel_mean, '  ±  ', PLcorel_std, '  nm', '    (', SIParam, ')'];
    BRcorel_mean = get(handles.ans_CBcorel_mean,'String');
    BRcorel_std = get(handles.ans_CBcorel_std,'String');
    BRcorel = ['Bending Rigidity        =   ', BRcorel_mean, '  ±  ', BRcorel_std, '  N.m^2'];
    coeffD_corel = get(handles.ans_rsqcorel_mean,'String');
    coeffD_corel = ['Determination Coeff   =   ',coeffD_corel];
    
    % mTextBox = uicontrol('style','text');
    %set(mTextBox,'Position',[120 390 390 60]);
    mTextBox_co = text;
    set(mTextBox_co,'units','normalized');
    set(mTextBox_co,'Position',[0.3 0.9 0]);
    set(mTextBox_co,'HorizontalAlignment', 'left');
    set(mTextBox_co,'fontsize',9);
    set(mTextBox_co,'color',[0 0 0 ]);
    set(mTextBox_co,'edgecolor',[1 1 1]);
    set(mTextBox_co,'String',{PLcorel, BRcorel, coeffD_corel });
    %set(mTextBox,'BackgroundColor',[1 1 1]);
    %else
    %end
    
    
    
    % BIN Deviations starts here
    
    new2 = [];
    for i = 1:z
        new2 = [new2; matdev(1,i).deviafib];
    end
    vecdev = new2;
    vecdev = real(vecdev);
    vecdev(any(isnan(vecdev),2),:)=[];
    
    %BIN DEVIATIONS DATA
    v = vecdev(:,2); %split into x and y
    w = vecdev(:,1);
    botEdge = round(min(v)); % define limits
    topEdge = maxval_dev; % define limits
    %numBins = 100; % define number of bins
    binEdges = linspace(botEdge, topEdge, numBins_dev+1);
    [~,whichBin] = histc(v, binEdges);
    for i = 1:numBins_dev
        flagBinMembers = (whichBin == i);
        binMembers     = w(flagBinMembers);
        binMean_dev(i) = mean((binMembers).^2);
        %---setup the length as the center of two separate bin edges
        length_dev(i) = binEdges(i) + ( ( (topEdge - botEdge) / numBins_dev) /2 );
    end
    
    devia_length = length_dev;
    devia_binMean = binMean_dev;
    %value = num2str(numBins_dev);
    %filename=[value,'bins_devia'];
    
    %---------------------------------------------------
    %MAKE FIGURE for the deviations fit
    
    checkboxdev = handles.checkboxdevia;
    
    %if checkboxdev == 1
    % Fit: 'deviafit'.
    [xData, yData] = prepareCurveData( devia_length, devia_binMean);
    
    % Set up fittype and options.
    ft = fittype( 'x^3/48/p', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( ft );
    opts.Display = 'Off';
    opts.Lower = -Inf;
    opts.StartPoint = 500;
    opts.Upper = Inf;
    [fitresult, gof] = fit( xData, yData, ft, opts );
    
    h3 = figure('Color',[1 1 1], 'Name', 'Deviations from secant midpoint (all chains)' );
    set(h3, 'Position', [1150 400  560  420])
    scatter( xData, yData );
    hold on
    plot(fitresult);
    legend('off')
    xlabel( 'Secant length (nm)' );
    ylabel( 'Mean-square deviations (nm^2)' );
    grid off
    
    PLdevia_mean = get(handles.ans_pdevia_mean,'String');
    PLdevia_std = get(handles.ans_pdevia_std,'String');
    
    PLdevia = ['Persistence Length   =   ', PLdevia_mean, '  ±  ', PLdevia_std, '  nm', '    (', SIParam, ')'];
    BRdevia_mean = get(handles.ans_CBdevia_mean,'String');
    BRdevia_std = get(handles.ans_CBdevia_std,'String');
    BRdevia = ['Bending Rigidity        =   ', BRdevia_mean, '  ±  ', BRdevia_std, '  N.m^2'];
    coeffD_devia = get(handles.ans_rsqdevia_mean,'String');
    coeffD_devia = ['Determination Coeff   =   ',coeffD_devia];
    
    %set(mTextBox,'Position',[120 390 390 60]);
    mTextBox_dev = text;
    set(mTextBox_dev,'units','normalized');
    set(mTextBox_dev,'Position',[0.05 0.9 0]);
    set(mTextBox_dev,'HorizontalAlignment', 'left');
    set(mTextBox_dev,'fontsize',9);
    set(mTextBox_dev,'color',[0 0 0 ]);
    set(mTextBox_dev,'edgecolor',[1 1 1]);
    set(mTextBox_dev,'String',{PLdevia, BRdevia, coeffD_devia });
    
    %else
    %end
    
    
    
    % get data and pass them to handles for future saving
    handles.corel_length = corel_length;
    handles.corel_binMean = corel_binMean;
    handles.devia_length = devia_length;
    handles.devia_binMean = devia_binMean;
    handles.worm_length = worm_length;
    handles.worm_binMean = worm_binMean;
    handles.vecworm = vecworm;
    handles.vecdev = vecdev;
    handles.vectan = vectan;
    
    handles.maxworm = round(max(x));
    handles.maxdevia = round(max(v));
    handles.maxcorel = round(max(a));
    
catch err
    ;
end
    
    guidata(hObject,handles);


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)

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


function closeGUI(src,evnt)
%this function is called when the user attempts to close the GUI window

%this command brings up the close window dialog
selection = questdlg('Are you sure you want to exit?',...
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


function z = round2(x,y)
%defensive programming
error(nargchk(2,2,nargin))
error(nargoutchk(0,1,nargout))
if numel(y)>1
  error('Y must be scalar')
end
z = round(x/y)*y;


function ans_nfibrilsdevia_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of ans_nfibrilsdevia as text
%        str2double(get(hObject,'String')) returns contents of ans_nfibrilsdevia as a double

% --- Executes during object creation, after setting all properties.
function ans_nfibrilsdevia_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ans_nfibrilscorel_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of ans_nfibrilscorel as text
%        str2double(get(hObject,'String')) returns contents of ans_nfibrilscorel as a double

% --- Executes during object creation, after setting all properties.
function ans_nfibrilscorel_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in btn_savePLdata_only.
function btn_savePLdata_only_Callback(hObject, eventdata, handles)
% hObject    handle to btn_savePLdata_only (see GCBO)

x = handles.ntot_fibrils;
contour_lengths = handles.struc.contour_lengths;
intervals_means = handles.struc.intervals_means;

matwo = handles.struc.mat_wormlike;
new = [];
for i = 1:x
    contourL_vec(1,i) = contour_lengths(1,i).contourL;
    intervalmeans_vec(1,i) = intervals_means(1,i).meanitv_fib;
    new = [new; matwo(1,i).wormfib];
end
totaldatapoints = length(new);

contourL_mean = round(mean(contourL_vec));
contourL_sd = round(std(contourL_vec));
contourL_sum = sum(contourL_vec);
contourL_sum = round2(contourL_sum/1000,0.1);
intervals_mean = round(mean(intervalmeans_vec));
intervals_sd = round(std(intervalmeans_vec));

try
    contourL_segmentsep = handles.kurtosisSegL;
catch err
    msgbox(sprintf('It is strongly advised to take a 2D equilibration test before saving your data.  Please hit the Kurtosis button in the 2D Equilibration Tests panel.'),...
        'Warning','Warn') ;
end

contourL_intvmean = intervals_mean;
contourL_intvsd = intervals_sd;
contourL_totpoints = totaldatapoints;
contourL_vec = contourL_vec.' ;

wormPersL_mean = str2double(get(handles.ans_pworm_mean,'String'));
wormPersL_sd = str2double(get(handles.ans_pworm_std,'String'));
wormCB_mean = str2double(get(handles.ans_CBworm_mean,'String'));
wormCB_sd = str2double(get(handles.ans_CBworm_std,'String'));
wormRsquare = str2double(get(handles.ans_rsqworm_mean,'String'));
wolow = str2double(get(handles.set_nb_lowfit_worm,'String'));
wotop = str2double(get(handles.set_nb_upfit_worm,'String'));
wormEdges = [wolow wotop];

corelPersL_mean = str2double(get(handles.ans_pcorel_mean,'String'));
corelPersL_sd = str2double(get(handles.ans_pcorel_std,'String'));
corelCB_mean = str2double(get(handles.ans_CBcorel_mean,'String'));
corelCB_sd = str2double(get(handles.ans_CBcorel_std,'String'));
corelRsquare = str2double(get(handles.ans_rsqcorel_mean,'String'));
colow = str2double(get(handles.set_nb_lowfit_corel,'String'));
cotop = str2double(get(handles.set_nb_upfit_corel,'String'));
corelEdges = [colow cotop];

deviaPersL_mean = str2double(get(handles.ans_pdevia_mean,'String'));
deviaPersL_sd = str2double(get(handles.ans_pdevia_std,'String'));
deviaCB_mean = str2double(get(handles.ans_CBdevia_mean,'String'));
deviaCB_sd = str2double(get(handles.ans_CBdevia_std,'String'));
deviaRsquare = str2double(get(handles.ans_rsqdevia_mean,'String'));
devlow = str2double(get(handles.set_nb_lowfit_devia,'String'));
devtop = str2double(get(handles.set_nb_upfit_devia,'String'));
deviaEdges = [devlow devtop];

try
    worm_length = handles.worm_length;
    worm_binsqrmean = handles.worm_binMean;
    A = [worm_length; worm_binsqrmean];
    worm_graphCoord = A';
    
    corel_length = handles.corel_length;
    corel_binmean = handles.corel_binMean;
    B = [corel_length; corel_binmean];
    corel_graphCoord = B';
    
    devia_length = handles.devia_length;
    devia_binsqrmean = handles.devia_binMean;
    C = [devia_length; devia_binsqrmean];
    devia_graphCoord = C';
catch err
    msgbox(sprintf('Are you sure you do not want to collect data for all fibrils? Whatever your analysis was, you should press the Sum & Bin All Chains button before saving.'),...
        'Warning','Warn') ;
end


try
    test2D_length = handles.test2D_length;
    test2D_binmean = handles.test2D_binmean;
    D = [test2D_length; test2D_binmean];
    surfint_test2D = D';
catch err
    ;
end

try
    test2Dfractal_contour = handles.fractal_contour;
    test2Dfractal_E2E = handles.fractal_E2E; 
    E = [test2Dfractal_contour; test2Dfractal_E2E];
    Test2D_fractal = E'; 
catch err
    ;
end

contourL_max = str2double(get(handles.ans_longestfibL,'String'));

%surfint_param = str2double(get(handles.set_surfintparam,'String'));
        toggle23 = get(handles.toggle_2D3Dswitcher,'String');
        TF = strcmp(toggle23,'Switch to Noneq.');
        if TF == 1
            surfint_param = 'Yes';
        else
            surfint_param = 'No';
        end

% [filename, pathname] = uiputfile('*.mat','save workspace variables as');
% newfilename = fullfile(pathname, filename);
% save(newfilename, 'surfint_param','surfint_test2D',...
%    'contourL_mean','contourL_sd','contourL_sum',...
%    'contourL_intvmean','contourL_intvsd',...
%    'contourL_segmentsep','contourL_vec','contourL_totpoints','contourL_max',...
%    'wormPersL_mean','wormPersL_sd','wormRsquare','wormEdges',...
%    'corelPersL_mean','corelPersL_sd','corelRsquare','corelEdges',...
%    'deviaPersL_mean','deviaPersL_sd','deviaRsquare','deviaEdges',...
%    'corel_graphCoord','devia_graphCoord','worm_graphCoord');

fileID1 = fopen('PLoutput_main.txt','w');
fprintf(fileID1,'No of chains: %6u\r\n \r\n \r\n',size(contourL_vec,1));
fprintf(fileID1,'Contour length: %6u +/- %u nm\r\n',contourL_mean,contourL_sd);
fprintf(fileID1,'Contour length, max value: %6u nm\r\n',contourL_max);
fprintf(fileID1,'Knots interval: %6u +/- %u nm\r\n',contourL_intvmean,contourL_intvsd);
fprintf(fileID1,'Contour, sum of all chains: %6.1f um\r\n',contourL_sum);
fprintf(fileID1,'Contour, Ntot points: %6u\r\n \r\n \r\n',contourL_totpoints);
fprintf(fileID1,'Full equilibration in 2D: %s\r\n \r\n',surfint_param);
fprintf(fileID1,'PL, Contour vs End2End (R-square): %6u +/- %u nm (%.4f)\r\n',wormPersL_mean,wormPersL_sd,wormRsquare);
fprintf(fileID1,'fit, Contour vs End2End, boundaries: %6u ---> %u nm\r\n \r\n',wormEdges(1),wormEdges(2));
fprintf(fileID1,'PL, Correlation function (R-square): %6u +/- %u nm (%.4f)\r\n',corelPersL_mean,corelPersL_sd,corelRsquare);
fprintf(fileID1,'fit, Correlation function, boundaries: %6u ---> %u nm\r\n \r\n',corelEdges(1),corelEdges(2));
fprintf(fileID1,'PL, Deviations to secant midpoint (R-square): %6u +/- %u nm (%.4f)\r\n',deviaPersL_mean,deviaPersL_sd,deviaRsquare);
fprintf(fileID1,'fit, Deviations to secant midpoint, boundaries: %6u ---> %u nm\r\n',deviaEdges(1),deviaEdges(2));
fclose(fileID1);

fileID2 = fopen('PLoutput_contour.txt','w');
fprintf(fileID2,'%6s\r\n' ,'ContourL_nm');
fprintf(fileID2,'%10.1f\r\n',contourL_vec);
fclose(fileID2);

JGVJHGBHKG = contourL_segmentsep; % this line avoids the creation of .txt file that won't close when error is generated
% fileID0 = fopen('PLoutput_segmentsep.txt','w');
% fprintf(fileID0,'%6s\r\n' ,'Segmentsep_nm');
% fprintf(fileID0,'%10.1f\r\n',contourL_segmentsep);
% fclose(fileID0);


A = transpose(worm_graphCoord);
fileID3 = fopen('PLoutput_graphs_CLvsE2E.txt','w');
fprintf(fileID3,'%16s \t %16s\r\n','CLvsE2E_x_nm','CLvsE2E_y_nm^2');
fprintf(fileID3,'%16.3f \t %16.3f\r\n',A);
fclose(fileID3);

B = transpose(corel_graphCoord);
fileID4 = fopen('PLoutput_graphs_correlations.txt','w');
fprintf(fileID4,'%16s \t %16s\r\n','correlations_x_nm','correlations_y');
fprintf(fileID4,'%16.3f \t %16.8f\r\n',B);
fclose(fileID4);

C = transpose(devia_graphCoord);
fileID5 = fopen('PLoutput_graphs_deviations.txt','w');
fprintf(fileID5,'%16s \t %16s\r\n','deviations_x_nm','deviations_y_nm^2');
fprintf(fileID5,'%16.3f \t %16.3f\r\n',C);
fclose(fileID5);

D = transpose(surfint_test2D);
fileID5 = fopen('PLoutput_graphs_Test2D_kurtosis.txt','w');
fprintf(fileID5,'%16s \t %6s\r\n','test2Dkurt_x_nm','test2Dkurt_y');
fprintf(fileID5,'%16.3f \t %8.3f\r\n',D);
fclose(fileID5);

try
E = transpose(Test2D_fractal);
fileID5 = fopen('PLoutput_graphs_Test2D_fractalexp.txt','w');
fprintf(fileID5,'%16s \t %6s\r\n','test2Dfrac_x_nm','test2Dfrac_y_nm');
fprintf(fileID5,'%16.3f \t %8.3f\r\n',E);
fclose(fileID5);
catch err
    ;
end


% message confirmation
msgbox(sprintf('Output files have been successfully generated in the folder that contains the Easyworm executables.'),...
       'Confirmation Note','Help')
   

guidata(hObject,handles);



% --- Executes on button press in btn_saveAll.
function btn_saveAll_Callback(hObject, eventdata, handles)
% hObject    handle to btn_saveAll (see GCBO)
% hObject    handle to btn_savePLdata_only (see GCBO)

x = handles.ntot_fibrils;
contour_lengths = handles.struc.contour_lengths;
intervals_means = handles.struc.intervals_means;
% intervals_sd = handles.struc.intervals_sd;

matwo = handles.struc.mat_wormlike;
new = [];

for i = 1:x
    contourL_vec(1,i) = contour_lengths(1,i).contourL;
    intervalmeans_vec(1,i) = intervals_means(1,i).meanitv_fib;
    % intervalstd_vec(1,i) = intervals_sd(1,i).meansd_fib;
    new = [new; matwo(1,i).wormfib];
end
totaldatapoints = length(new);

contourL_mean = round(mean(contourL_vec));
contourL_sd = round(std(contourL_vec));
contourL_sum = sum(contourL_vec);
contourL_sum = round2(contourL_sum/1000,0.1);
intervals_mean = round(mean(intervalmeans_vec));
intervals_sd = round(std(intervalmeans_vec));

try
    contourL_segmentsep = handles.kurtosisSegL;
catch err
    msgbox(sprintf('It is strongly advised to take an equilibration-in-2D test before saving your data.  Please hit the Kurtosis button in the 2D Equilibration Tests panel.'),...
        'Warning','Warn') ;
end

contourL_intvmean = intervals_mean;
contourL_intvsd = intervals_sd;
contourL_totpoints = totaldatapoints;
contourL_vec = contourL_vec.' ;

wormPL = str2double(get(handles.ans_pworm_mean,'String'));
wormPL_sd = str2double(get(handles.ans_pworm_std,'String'));
wormCB = str2double(get(handles.ans_CBworm_mean,'String'));
wormCB_sd = str2double(get(handles.ans_CBworm_std,'String'));
wormRsquare = str2double(get(handles.ans_rsqworm_mean,'String'));
wolow = str2double(get(handles.set_nb_lowfit_worm,'String'));
wotop = str2double(get(handles.set_nb_upfit_worm,'String'));
wormEdges = [wolow wotop];

corelPL = str2double(get(handles.ans_pcorel_mean,'String'));
corelPL_sd = str2double(get(handles.ans_pcorel_std,'String'));
corelCB = str2double(get(handles.ans_CBcorel_mean,'String'));
corelCB_sd = str2double(get(handles.ans_CBcorel_std,'String'));
corelRsquare = str2double(get(handles.ans_rsqcorel_mean,'String'));
colow = str2double(get(handles.set_nb_lowfit_corel,'String'));
cotop = str2double(get(handles.set_nb_upfit_corel,'String'));
corelEdges = [colow cotop];

deviaPL = str2double(get(handles.ans_pdevia_mean,'String'));
deviaPL_sd = str2double(get(handles.ans_pdevia_std,'String'));
deviaCB = str2double(get(handles.ans_CBdevia_mean,'String'));
deviaCB_sd = str2double(get(handles.ans_CBdevia_std,'String'));
deviaRsquare = str2double(get(handles.ans_rsqdevia_mean,'String'));
devlow = str2double(get(handles.set_nb_lowfit_devia,'String'));
devtop = str2double(get(handles.set_nb_upfit_devia,'String'));
deviaEdges = [devlow devtop];

try
    worm_length = handles.worm_length;
    worm_binsqrmean = handles.worm_binMean;
    A = [worm_length; worm_binsqrmean];
    worm_graphCoord = A';
    
    corel_length = handles.corel_length;
    corel_binmean = handles.corel_binMean;
    B = [corel_length; corel_binmean];
    corel_graphCoord = B';
    
    devia_length = handles.devia_length;
    devia_binsqrmean = handles.devia_binMean;
    C = [devia_length; devia_binsqrmean];
    devia_graphCoord = C';
catch err
    msgbox(sprintf('Are you sure you do not want to collect data for all fibrils? Whatever your analysis was, you should press the Sum & Bin All Chains button before saving.'),...
        'Warning','Warn') ;
end

try
    test2D_length = handles.test2D_length;
    test2D_binmean = handles.test2D_binmean;
    D = [test2D_length; test2D_binmean];
    surfint_test2D = D';
catch err
    ;
end

try
    test2Dfractal_contour = handles.fractal_contour;
    test2Dfractal_E2E = handles.fractal_E2E; 
    E = [test2Dfractal_contour; test2Dfractal_E2E];
    Test2D_fractal = E'; 
catch err
    ;
end


contourL_max = str2double(get(handles.ans_longestfibL,'String'));

% Below is all the stuff saved in addition to PLength data
% (i.e. Young's modulus etc.)

final_PL = str2double(get(handles.set_PL_mean,'String'));
final_PL_sd = str2double(get(handles.set_PL_SD,'String'));
final_BendR = str2double(get(handles.bendr_mean,'String'));
final_BendR_sd = str2double(get(handles.bendr_SD,'String'));

helic_mI = str2double(get(handles.ans_helic_inertia_mean,'String'));
helic_mI_sdplus = str2double(get(handles.ans_helicIplus,'String'));
helic_mI_sdminus = str2double(get(handles.ans_helicIminus,'String'));

tape_mI = str2double(get(handles.ans_ribbonI,'String'));
tape_mI_sdplus = str2double(get(handles.ans_ribbonIplus,'String'));
tape_mI_sdminus = str2double(get(handles.ans_ribbonIminus,'String'));

final_ElasMod = str2double(get(handles.youngm_mean,'String'));
final_ElasMod_sd = str2double(get(handles.youngm_SD,'String'));

Helic_EM = str2double(get(handles.ans_helic_youngm,'String'));
Helic_EMplus = str2double(get(handles.ans_helicYplus,'String'));
Helic_EMminus = str2double(get(handles.ans_helicYminus,'String'));

Ribbon_Nb_protofil = str2double(get(handles.set_number_protofil,'String'));
Ribbon_EM = str2double(get(handles.ans_ribbon_youngm,'String'));
Ribbon_EMplus = str2double(get(handles.ans_ribbon_Eplus,'String'));
Ribbon_EMminus = str2double(get(handles.ans_ribbon_Eminus,'String'));

Height = str2double(get(handles.set_ribfibH,'String'));
Height_sd = str2double(get(handles.set_ribfibH_SD,'String'));
Width = str2double(get(handles.set_width,'String'));
Width_sd = str2double(get(handles.set_width_sd,'String'));

%surfint_param = str2double(get(handles.set_surfintparam,'String'));
        toggle23 = get(handles.toggle_2D3Dswitcher,'String');
        TF = strcmp(toggle23,'Switch to Noneq.');
        if TF == 1
            surfint_param = 'Yes';
        else
            surfint_param = 'No';
        end

% [filename, pathname] = uiputfile('*.mat','Save Workspace Variables As');
% newfilename = fullfile(pathname, filename);
% save(newfilename, 'surfint_param','surfint_test2D',...
%    'contourL_mean','contourL_sd','contourL_sum',...
%    'contourL_intvmean','contourL_intvsd',...
%    'contourL_segmentsep','contourL_vec','contourL_totpoints','contourL_max',...
%    'wormPL','wormPL_sd','wormRsquare','wormEdges',...
%    'corelPL','corelPL_sd','corelRsquare','corelEdges',...
%    'deviaPL','deviaPL_sd','deviaRsquare','deviaEdges',...
%    'corel_graphCoord','devia_graphCoord','worm_graphCoord',...
%    'final_PL','final_PL_sd','final_BendR','final_BendR_sd',...
%    'helic_mI','helic_mI_sdplus','helic_mI_sdminus',...
%    'final_ElasMod','final_ElasMod_sd',...
%    'Helic_EM','Helic_EMplus','Helic_EMminus',...
%    'Ribbon_Nb_protofil','Ribbon_EM','Ribbon_EMplus','Ribbon_EMminus',...
%    'Height','Height_sd','Width','Width_sd');

fileID1 = fopen('output_main.txt','w');
fprintf(fileID1,'No of chains = %6u\r\n \r\n \r\n',size(contourL_vec,1));
fprintf(fileID1,'Contour length = %6u +/- %u nm\r\n',contourL_mean,contourL_sd);
fprintf(fileID1,'Contour length, max value = %6u nm\r\n',contourL_max);
fprintf(fileID1,'Knots interval = %6u +/- %u nm\r\n',contourL_intvmean,contourL_intvsd);
fprintf(fileID1,'Contour, sum of all chains = %6.1f um\r\n',contourL_sum);
fprintf(fileID1,'Contour, Ntot points = %6u\r\n \r\n \r\n',contourL_totpoints);
fprintf(fileID1,'Full equilibration in 2D: %s\r\n \r\n',surfint_param);
fprintf(fileID1,'PL, Contour vs End2End (R-square) = %6u +/- %u nm (%.4f)\r\n',wormPL,wormPL_sd,wormRsquare);
fprintf(fileID1,'fit, Contour vs End2End, boundaries: %6u ---> %u nm\r\n \r\n',wormEdges(1),wormEdges(2));
fprintf(fileID1,'PL, Correlation function (R-square) = %6u +/- %u nm (%.4f)\r\n',corelPL,corelPL_sd,corelRsquare);
fprintf(fileID1,'fit, Correlation function, boundaries: %6u ---> %u nm\r\n \r\n',corelEdges(1),corelEdges(2));
fprintf(fileID1,'PL, Deviations to secant midpoint (R-square) = %6u +/- %u nm (%.4f)\r\n',deviaPL,deviaPL_sd,deviaRsquare);
fprintf(fileID1,'fit, Deviations to secant midpoint, boundaries: %6u ---> %u nm\r\n \r\n \r\n',deviaEdges(1),deviaEdges(2));



fprintf(fileID1,'Final PL = %6u +/- %u nm\r\n',final_PL,final_PL_sd);
fprintf(fileID1,'Final BR = %6.5u +/- %.5u N.m^2\r\n \r\n',final_BendR,final_BendR_sd);

fprintf(fileID1,'Height = %6.2f +/- %.2f nm\r\n',Height,Height_sd);
fprintf(fileID1,'Width = %6.2f +/- %.2f nm\r\n \r\n \r\n',Width,Width_sd);

fprintf(fileID1,'(if Width = 0, I and E below are calculated using helical model)\r\n\r\nHelical *or* Ellipsoidal I: I_hel = %6.5u m^4\r\nI_hel+ = %6.5u m^4\r\nI_hel- = %6.5u m^4\r\n',helic_mI,helic_mI_sdplus,helic_mI_sdminus);
fprintf(fileID1,'\r\nRibbon, Tape I: I_rib = %6.5u m^4\r\nI_rib+ = %6.5u m^4\r\nI_rib- = %6.5u m^4\r\n \r\n \r\n',tape_mI,tape_mI_sdplus,tape_mI_sdminus);

fprintf(fileID1,'Helical *or* Ellipsoidal modulus: E_hel = %6.1f MPa\r\nE_hel+ = %6.1f MPa\r\nE_hel- = %6.1f MPa\r\n',Helic_EM,Helic_EMplus,Helic_EMminus);
fprintf(fileID1,'\r\nRibbon, Tape modulus (%u protofilaments): E_rib = %6.1f MPa\r\nE_rib+ = %6.1f MPa\r\nE_rib- = %6.1f MPa\r\n',Ribbon_Nb_protofil,Ribbon_EM,Ribbon_EMplus,Ribbon_EMminus);
fclose(fileID1);

fileID2 = fopen('output_contour.txt','w');
fprintf(fileID2,'%6s\r\n' ,'ContourL_nm');
fprintf(fileID2,'%10.1f\r\n',contourL_vec);
fclose(fileID2);


JGVJHGBHKG = contourL_segmentsep; % this line avoids the creation of .txt file that won't close when error is generated
% fileID0 = fopen('output_segmentsep.txt','w');
% fprintf(fileID0,'%6s\r\n' ,'Segmentsep_nm');
% fprintf(fileID0,'%10.1f\r\n',contourL_segmentsep);
% fclose(fileID0);

A = transpose(worm_graphCoord);
fileID3 = fopen('output_graphs_CLvsE2E.txt','w');
fprintf(fileID3,'%16s \t %16s\r\n','CLvsE2E_x_nm','CLvsE2E_y_nm^2');
fprintf(fileID3,'%16.3f \t %16.3f\r\n',A);
fclose(fileID3);

B = transpose(corel_graphCoord);
fileID4 = fopen('output_graphs_correlations.txt','w');
fprintf(fileID4,'%16s \t %16s\r\n','correlations_x_nm','correlations_y');
fprintf(fileID4,'%16.3f \t %16.8f\r\n',B);
fclose(fileID4);

C = transpose(devia_graphCoord);
fileID5 = fopen('output_graphs_deviations.txt','w');
fprintf(fileID5,'%16s \t %16s\r\n','deviations_x_nm','deviations_y_nm^2');
fprintf(fileID5,'%16.3f \t %16.3f\r\n',C);
fclose(fileID5);

D = transpose(surfint_test2D);
fileID5 = fopen('output_graphs_Test2D_kurtosis.txt','w');
fprintf(fileID5,'%16s \t %6s\r\n','test2Dkurt_x_nm','test2Dkurt_y');
fprintf(fileID5,'%16.3f \t %8.3f\r\n',D);
fclose(fileID5);

try
E = transpose(Test2D_fractal);
fileID5 = fopen('output_graphs_Test2D_fractalexp.txt','w');
fprintf(fileID5,'%16s \t %6s\r\n','test2Dfrac_x_nm','test2Dfrac_y_nm');
fprintf(fileID5,'%16.3f \t %8.3f\r\n',E);
fclose(fileID5);
catch err
    ;
end


% message confirmation
msgbox(sprintf('Output files have been successfully generated in the folder that contains the Easyworm executables.'),...
       'Confirmation Note','Help')

guidata(hObject,handles);


function txt_threshold_max_Callback(hObject, eventdata, handles)
% hObject    handle to set_PL_mean (see GCBO)
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_threshold_max_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function txt_threshold_min_Callback(hObject, eventdata, handles)
% hObject    handle to set_PL_mean (see GCBO)
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_threshold_min_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function txt_nthknot_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of txt_nthknot as text
%        str2double(get(hObject,'String')) returns contents of txt_nthknot as a double

% --- Executes during object creation, after setting all properties.
function txt_nthknot_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function txt_linewidth_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of txt_linewidth as text
%        str2double(get(hObject,'String')) returns contents of txt_linewidth as a double

% --- Executes during object creation, after setting all properties.
function txt_linewidth_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function set_PL_mean_Callback(hObject, eventdata, handles)
% hObject    handle to set_PL_mean (see GCBO)
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function set_PL_mean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_PL_mean (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function set_helicfib_dia_mean_Callback(hObject, eventdata, handles)
% hObject    handle to set_helicfib_dia_mean (see GCBO)

% --- Executes during object creation, after setting all properties.
function set_helicfib_dia_mean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_helicfib_dia_mean (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function set_PL_SD_Callback(hObject, eventdata, handles)
% hObject    handle to set_PL_SD (see GCBO)
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function set_PL_SD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_PL_SD (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bendr_mean_Callback(hObject, eventdata, handles)
% hObject    handle to bendr_mean (see GCBO)
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function bendr_mean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bendr_mean (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bendr_SD_Callback(hObject, eventdata, handles)
% hObject    handle to bendr_SD (see GCBO)
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function bendr_SD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bendr_SD (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function inertia_mean_Callback(hObject, eventdata, handles)
% hObject    handle to inertia_mean (see GCBO)
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function inertia_mean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inertia_mean (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function inertia_SDplus_Callback(hObject, eventdata, handles)
% hObject    handle to inertia_SDplus (see GCBO)
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function inertia_SDplus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inertia_SDplus (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function youngm_SD_Callback(hObject, eventdata, handles)
% hObject    handle to youngm_SD (see GCBO)
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function youngm_SD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to youngm_SD (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function inertia_SDminus_Callback(hObject, eventdata, handles)
% hObject    handle to inertia_SDminus (see GCBO)
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function inertia_SDminus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inertia_SDminus (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function set_ribfibH_Callback(hObject, eventdata, handles)
% hObject    handle to set_ribfibH (see GCBO)
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function set_ribfibH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_ribfibH (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function set_ribfibH_SD_Callback(hObject, eventdata, handles)
% hObject    handle to set_ribfibH_SD (see GCBO)
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function set_ribfibH_SD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_ribfibH_SD (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function set_number_protofil_Callback(hObject, eventdata, handles)
% hObject    handle to set_number_protofil (see GCBO)
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function set_number_protofil_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_number_protofil (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function set_protofil_dia_SD_Callback(hObject, eventdata, handles)
% hObject    handle to set_protofil_dia_SD (see GCBO)
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function set_protofil_dia_SD_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function set_protofil_dia_Callback(hObject, eventdata, handles)
% hObject    handle to set_protofil_dia (see GCBO)
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function set_protofil_dia_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_protofil_dia (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function set_helicfib_dia_SD_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of set_helicfib_dia_SD as text
%        str2double(get(hObject,'String')) returns contents of set_helicfib_dia_SD as a double
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function set_helicfib_dia_SD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_helicfib_dia_SD (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function youngm_mean_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of youngm_mean as text
%        str2double(get(hObject,'String')) returns contents of youngm_mean as a double
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function youngm_mean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to youngm_mean (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_calc_helical.
function pushbutton_calc_helical_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_calc_helical (see GCBO)
gocalc_helic(hObject, handles);
fibwidth = str2double(get(handles.set_ellipsoid_width_mean,'String'));
if fibwidth ~= 0
    msgbox(sprintf('Note that the chain width was set to a non-zero value.  You might want to push the ''Ellipsoidal'' button next time, otherwise I will keep resetting the width to zero.'),...
       'Warning','Warn')
   set(handles.set_ellipsoid_width_mean,'string','0')
   set(handles.set_ellipsoid_width_SD,'string','0')
   
else
end



% --- Executes on button press in pushbut_ellipsoidcalc.
function pushbut_ellipsoidcalc_Callback(hObject, eventdata, handles)
% hObject    handle to pushbut_ellipsoidcalc (see GCBO)
gocalc_ellipsoid(hObject, handles);


function gocalc_ellipsoid(hObject, handles)
% Calculate area moment of inertia from AFM height data
% USING helical model of fibrils
% + derive young's modulus value 
persisLmean_nm = str2double(get(handles.set_PL_mean,'String'));
persisLSD_nm = str2double(get(handles.set_PL_SD,'String'));
PLmean = persisLmean_nm ./1E9 ; % persistence length in meters
PLsd = persisLSD_nm ./1E9 ;

kB = 1.3806503e-23 ;  %Boltzmann constant

Temp = str2double(get(handles.set_Temperature,'String'));
T = 273.15 + Temp; % Temperature in deg Kelvin
BendingR_mean = kB * T * PLmean ;  % bending rigidity in N.m^2
BendingR_SD = kB * T * PLsd ;  
% set(handles.bendr_mean,'String',round2(BendingR_mean,1E-33)) ;
set(handles.bendr_mean,'String',BendingR_mean) ;
set(handles.bendr_SD,'String',round2(BendingR_SD,1E-33)) ;

%ellipH means the height of the ellipse.  ellipW = its width
ellipH_mean = str2double(get(handles.set_helicfib_dia_mean,'String'));
ellipH_mean = ellipH_mean ./1E9; % converts in meters (from nm)
ellipHrad = ellipH_mean ./2 ;
ellipW_mean = str2double(get(handles.set_ellipsoid_width_mean,'String'));
ellipW_mean = ellipW_mean ./1E9;
ellipWrad = ellipW_mean ./2 ;
Iellip = (pi.*ellipWrad.*(ellipHrad^3))./4 ;
% set(handles.ans_helic_inertia_mean,'String',round2(IH,10e-39)) ;
% set(handles.inertia_mean,'String',round2(IH,10e-39))
set(handles.ans_helic_inertia_mean,'String',Iellip) ;
set(handles.inertia_mean,'String',Iellip)

ellipH_SD = str2double(get(handles.set_helicfib_dia_SD,'String')) ;
ellipH_SD = ellipH_SD ./1E9; % converts in meters (from nm)
ellipW_SD = str2double(get(handles.set_ellipsoid_width_SD,'String')) ;
ellipW_SD = ellipW_SD ./1E9; % converts in meters (from nm)
ellipH_radplus = (ellipH_mean + ellipH_SD) ./2 ;
ellipW_radplus = (ellipW_mean + ellipW_SD) ./2 ;
Iellip_plus = (pi.*ellipW_radplus.*(ellipH_radplus^3))./4 ;
% set(handles.ans_helicIplus,'String',round2(IHplus,10e-39)) ;
% set(handles.inertia_SDplus,'String',round2(abs(IH - IHplus),10e-39)) ;
set(handles.ans_helicIplus,'String',Iellip_plus) ;
set(handles.inertia_SDplus,'String',abs(Iellip - Iellip_plus)) ;

ellipH_radminus = (ellipH_mean - ellipH_SD) ./2 ;
ellipW_radminus = (ellipW_mean - ellipW_SD) ./2 ;
% the 2.66 factor accounts for the fact that centers of mass 
Iellip_minus = (pi.*ellipW_radminus.*(ellipH_radminus^3))./4/2.66 ;
% set(handles.ans_helicIminus,'String',round2(IHminus,10e-39)) ;
% set(handles.inertia_SDminus,'String',round2(abs(IH - IHminus),10e-39)) ;
set(handles.ans_helicIminus,'String',Iellip_minus) ;
set(handles.inertia_SDminus,'String',abs(Iellip - Iellip_minus)) ;


ellipY = BendingR_mean ./ Iellip ;
%converts to MPa and round to 0.1MPa
set(handles.ans_helic_youngm,'String',round2(ellipY ./1E6,0.1)) ;

ellipYplus = (BendingR_mean + BendingR_SD) ./ Iellip_minus ;
%helicYplus = abs(helicY - helicYplus);
set(handles.ans_helicYplus,'String',round2(ellipYplus ./1E6,0.1)) ;
ellipYminus = (BendingR_mean - BendingR_SD) ./ Iellip_plus ;
%helicYminus = abs(helicY - helicYminus);
set(handles.ans_helicYminus,'String',round2(ellipYminus ./1E6,0.1)) ;

% Here we set values to the ribbon / tape Fibrils panel
set(handles.set_ribfibH,'String',ellipH_mean*1E9) ;
set(handles.set_ribfibH_SD,'String',ellipH_SD*1E9) ;
set(handles.set_width,'String',ellipW_mean*1E9) ;
set(handles.set_width_sd,'String',ellipW_SD*1E9) ;

guidata(hObject, handles);


function gocalc_helic(hObject, handles)
% Calculate area moment of inertia from AFM height data
% USING helical model of fibrils
% + derive young's modulus value 
persisLmean_nm = str2double(get(handles.set_PL_mean,'String'));
persisLSD_nm = str2double(get(handles.set_PL_SD,'String'));
PLmean = persisLmean_nm ./1E9 ; % persistence length in meters
PLsd = persisLSD_nm ./1E9 ;

kB = 1.3806503e-23 ;  %Boltzmann constant

Temp = str2double(get(handles.set_Temperature,'String'));
T = 273.15 + Temp; % Temperature in deg Kelvin
BendingR_mean = kB * T * PLmean ;  % bending rigidity in N.m^2
BendingR_SD = kB * T * PLsd ;  

set(handles.bendr_mean,'String',BendingR_mean) ;
set(handles.bendr_SD,'String',round2(BendingR_SD,1E-33)) ;

helicdia_mean = str2double(get(handles.set_helicfib_dia_mean,'String'));
helicdia_mean = helicdia_mean ./1E9; % converts in meters (from nm)
helicrad = helicdia_mean ./2 ;
IH = (pi.*helicrad^4)./4 ;
set(handles.ans_helic_inertia_mean,'String',IH) ;
set(handles.inertia_mean,'String',IH)

helicdia_SD = str2double(get(handles.set_helicfib_dia_SD,'String')) ;
helicdia_SD = helicdia_SD ./1E9; % converts in meters (from nm)

helic_radplus = (helicdia_mean + helicdia_SD) ./2 ;
IHplus = (pi.*helic_radplus^4)./4 ;
set(handles.ans_helicIplus,'String',IHplus) ;
set(handles.inertia_SDplus,'String',abs(IH - IHplus)) ;

helic_radminus = (helicdia_mean - helicdia_SD) ./2 ;
IHminus = (pi.*helic_radminus^4)./4/2.66 ;
set(handles.ans_helicIminus,'String',IHminus) ;
set(handles.inertia_SDminus,'String',abs(IH - IHminus)) ;

helicY = BendingR_mean ./ IH ;
%converts to MPa and round to 0.1MPa
set(handles.ans_helic_youngm,'String',round2(helicY ./1E6,0.1)) ;
helicYplus = (BendingR_mean + BendingR_SD) ./ IHminus ;
set(handles.ans_helicYplus,'String',round2(helicYplus ./1E6,0.1)) ;
helicYminus = (BendingR_mean - BendingR_SD) ./ IHplus ;
set(handles.ans_helicYminus,'String',round2(helicYminus ./1E6,0.1)) ;
set(handles.set_ribfibH,'String',helicdia_mean*1E9) ;
set(handles.set_ribfibH_SD,'String',helicdia_SD*1E9) ;

guidata(hObject, handles);


% --- Executes on button press in pushbutton_quickcalc.
function pushbutton_quickcalc_Callback(hObject, eventdata, handles) %#ok<*INUSL>
quickcalc(hObject,handles);

function quickcalc(hObject,handles)
% calculates E = CB * I and vice-versa
% QC = quick calc ; BR = bending rigidity
QC_BR_mean = str2double(get(handles.bendr_mean,'String'));
QC_BR_SD = str2double(get(handles.bendr_SD,'String'));

QC_minertia_mean = str2double(get(handles.inertia_mean,'String'));
QC_minertia_SDplus = str2double(get(handles.inertia_SDplus,'String'));
QC_minertia_SDminus = str2double(get(handles.inertia_SDminus,'String'));

QC_youngmod_mean = str2double(get(handles.youngm_mean,'String'));
QC_youngmod_SD = str2double(get(handles.youngm_SD,'String'));
QC_youngmod_mean = QC_youngmod_mean .*1E6;  % converts in Pa
QC_youngmod_SD = QC_youngmod_SD .*1E6;

if QC_minertia_mean ~=0
    % Below is condition to use the quickcalc one way or the other, (if
    % Bending R is zero, then that's what we're gonna calculate)
    if QC_BR_mean == 0
        QC_BR_mean = QC_minertia_mean .* QC_youngmod_mean ;
        set(handles.bendr_mean,'String',round2(QC_BR_mean,1E-33)) ;
        
        % Below is error calculation according to simple propagation of errors
        % [ formula is: dA/A = sqrt[(dB/B)^2 + (dC/c)^2] considering independant
        % random errors of dB and dC on B and C, looking for dA]
        QC_BR_SDplus = QC_BR_mean .* ( sqrt ( ...
            ( (QC_youngmod_SD./QC_youngmod_mean).^2) + ...
            ( ( (QC_minertia_mean - QC_minertia_SDplus) ./ QC_minertia_mean ).^2 )...
            )  );
        %Note: we use QC_minertia_mean - QC_minertia_SDplus above and
        %      QC_minertia_mean - QC_minertia_SDminus below because the
        %      error on minertia is not symetric
        QC_BR_SDminus = QC_BR_mean .* ( sqrt ( ...
            ( (QC_youngmod_SD./QC_youngmod_mean).^2) + ...
            ( ( (QC_minertia_mean - QC_minertia_SDminus) ./ QC_minertia_mean ).^2 )...
            )  );
        if QC_BR_SDplus >= QC_BR_SDminus %#ok<BDSCI>
            QC_BR_SD = QC_BR_SDplus;
        else QC_BR_SD = QC_BR_SDminus;  %Keep the largest of the error to print
        end
        set(handles.bendr_SD,'String',round2(QC_BR_SD,1E-33)) ;
        
        kB = 1.3806503e-23 ;  %Boltzmann constant
        
        Temp = str2double(get(handles.set_Temperature,'String'));
        T = 273.15 + Temp; % Temperature in deg Kelvin
        PL_SDminus_PropError = QC_BR_SDminus ./ (kB * T);
        PL_SDplus_PropError = QC_BR_SDplus ./ (kB * T);
        set(handles.set_PL_SD,'String',...
            round2(PL_SDplus_PropError*1E9,1)) ;
        set(handles.ans_PL_SDminus_PropError,'String',...
            round2(PL_SDminus_PropError*1E9,1)) ;
        set(handles.txt_minus,'String','-') ;
        set(handles.txt_plus,'String','+') ;
        
                      
        %Below is Persistence length error calculated with maximum
        %amplitude from both Young's modulus and second moment of area combined
        % (useful to check because the method above tend to lower down
        % the error)
        BRplus_alt = (QC_youngmod_mean + QC_youngmod_SD).* ...
            (QC_minertia_mean + QC_minertia_SDplus);
        BRminus_alt = (QC_youngmod_mean - QC_youngmod_SD).* ...
            (QC_minertia_mean - QC_minertia_SDminus);
        
        kB = 1.3806503e-23 ;  %Boltzmann constant  
        Temp = str2double(get(handles.set_Temperature,'String'));
        T = 273.15 + Temp; % Temperature in deg Kelvin
        PLplus_alt = BRplus_alt ./ (kB * T);
        PLminus_alt = BRminus_alt ./ (kB * T);
        
        set(handles.ans_PLplus,'String',round2(PLplus_alt*1E9,1)) ;
        set(handles.ans_PLminus,'String',round2(PLminus_alt*1E9,1)) ;
        set(handles.sign_plus,'String','PL+ = ') ;
        set(handles.sign_minus,'String','PL- =') ;
        
        
    elseif QC_youngmod_mean == 0 %#ok<BDSCI>
        QC_youngmod_mean = QC_BR_mean ./ QC_minertia_mean ;
        set(handles.youngm_mean,'String',round2(QC_youngmod_mean./1E6,0.1)) ;
        % below is error calculation
        QC_youngmod_SDplus = QC_youngmod_mean .* ( sqrt ( ...
            ( (QC_BR_SD ./ QC_BR_mean).^2) + ...
            ( ( (QC_minertia_mean - QC_minertia_SDplus) ./ QC_minertia_mean ).^2 )...
            )  );
        QC_youngmod_SDminus = QC_youngmod_mean .* ( sqrt ( ...
            ( (QC_BR_SD ./ QC_BR_mean).^2) + ...
            ( ( (QC_minertia_mean - QC_minertia_SDminus) ./ QC_minertia_mean ).^2 )...
            )  );
        if QC_youngmod_SDplus >= QC_youngmod_SDminus %#ok<BDSCI>
            QC_youngmod_SD = QC_youngmod_SDplus;
        else QC_youngmod_SD = QC_youngmod_SDminus;
        end
        % converts in MPa and round to 0.1 MPa
        set(handles.youngm_SD,'String',round2(QC_youngmod_SD ./1E6,0.1)) ;
    else
    end
else
end

% Set PL value in Elsatic Modulus panel according to CB
kB = 1.3806503e-23 ;  %Boltzmann constant
Temp = str2double(get(handles.set_Temperature,'String'));
T = 273.15 + Temp; % Temperature in deg Kelvin
PLmean = QC_BR_mean ./ (kB * T);  % bending rigidity in N.m^2
set(handles.set_PL_mean,'String',round2(PLmean*1E9,1)) ;

%Don't need this anymore because we have separate + and - error (see above)
% on PL value (also on BR but does not appear on the graph)
% kB = 1.3806503e-23 ;
% PLsd = QC_BR_SD ./ (kB * 298.15);
% set(handles.set_PL_SD,'String',round2(PLsd*1E9,1)) ;
      


guidata(hObject,handles);


% --- Executes on button press in btn_zeromod.
function btn_zeromod_Callback(hObject, eventdata, handles)
% hObject    handle to btn_zeromod (see GCBO)
set(handles.youngm_mean,'String',0) ;
set(handles.youngm_SD,'String',0) ;


% --- Executes on button press in btn_zerobendingr.
function btn_zerobendingr_Callback(hObject, eventdata, handles)
% hObject    handle to btn_zerobendingr (see GCBO)
set(handles.bendr_mean,'String',0) ;
set(handles.bendr_SD,'String',0) ;

set(handles.ans_PLplus,'String','') ;
set(handles.ans_PLminus,'String','') ;
set(handles.sign_plus,'String','') ;
set(handles.sign_minus,'String','') ;

set(handles.ans_PL_SDminus_PropError,'String','') ;
set(handles.txt_minus,'String','nm') ;
set(handles.txt_plus,'String','±') ;



% --- Executes on button press in pushbutton_calc_ribbon_height.
function pushbutton_calc_ribbon_height_Callback(hObject, eventdata, handles)
% Calculate area moment of inertia from AFM height data
% USING RIBBON model of fibrils (USing HEIGHT value determined by AFM
% and number of filaments by SEM
% + derive young's modulus value 
persisLmean_nm = str2double(get(handles.set_PL_mean,'String'));
persisLSD_nm = str2double(get(handles.set_PL_SD,'String'));
PLmean = persisLmean_nm ./1E9 ; % persistence length in meters
PLsd = persisLSD_nm ./1E9 ;

kB = 1.3806503e-23 ;  %Boltzmann constant
Temp = str2double(get(handles.set_Temperature,'String'));
T = 273.15 + Temp; % Temperature in deg Kelvin
BendingR_mean = kB * T * PLmean ;  % bending rigidity in N.m^2
BendingR_SD = kB * T * PLsd ;  
% set(handles.bendr_mean,'String',round2(BendingR_mean,1E-33)) ;
set(handles.bendr_mean,'String',BendingR_mean) ;
set(handles.bendr_SD,'String',round2(BendingR_SD,1E-33)) ;


fibH = str2double(get(handles.set_ribfibH,'String'));
fibH_sd = str2double(get(handles.set_ribfibH_SD,'String'));
fibH = fibH ./1E9 ; % Height in meters
fibH_sd = fibH_sd ./1E9 ;

nfil = str2double(get(handles.set_number_protofil,'String'));
nfil_sd = 1;

protofil_dia = str2double(get(handles.set_protofil_dia,'String'));
protofil_dia_sd = str2double(get(handles.set_protofil_dia_SD,'String'));
protofil_dia = protofil_dia ./1E9 ; % converts in meters
protofil_dia_sd = protofil_dia_sd ./1E9 ;

% calculation of second moment of area
Irib = fibH^3 .* (nfil .* protofil_dia) ./ 12 ;
set(handles.ans_ribbonI,'String',Irib) ;
set(handles.inertia_mean,'String',Irib) ;

% calculates upper bound I+
Iribplus = (fibH + fibH_sd)^3 .* ((nfil + nfil_sd) .* ...
    (protofil_dia + protofil_dia_sd) ) ./ 12 ;
set(handles.ans_ribbonIplus,'String',Iribplus) ; 
set(handles.inertia_SDplus,'String',abs(Irib - Iribplus)) ;


% calculates lower bound I-
a = nfil - nfil_sd;
if a ~=0
    Iribminus = (fibH - fibH_sd)^3 .* ((nfil - nfil_sd) .* ...
        (protofil_dia - protofil_dia_sd) ) ./ 12/2.66 ;
else
    Iribminus = (fibH - fibH_sd)^3 .* ( 1 .* ...
        (protofil_dia - protofil_dia_sd) ) ./ 12/2.66 ;
end

set(handles.ans_ribbonIminus,'String',Iribminus) ;
set(handles.inertia_SDminus,'String',abs(Irib - Iribminus)) ;

% calculate Young's modulus
ribbonY = BendingR_mean ./ Irib ;
%converts to MPa and round to 0.1MPa
set(handles.ans_ribbon_youngm,'String',round2(ribbonY ./1E6,0.1)) ;

ribbonEplus = (BendingR_mean + BendingR_SD) ./ Iribminus ;
set(handles.ans_ribbon_Eplus,'String',round2(ribbonEplus ./1E6,0.1)) ;
ribbonEminus = (BendingR_mean - BendingR_SD) ./ Iribplus ;
set(handles.ans_ribbon_Eminus,'String',round2(ribbonEminus ./1E6,0.1)) ;
set(handles.set_helicfib_dia_mean,'String',fibH*1E9) ;
set(handles.set_helicfib_dia_SD,'String',fibH_sd*1E9) ;

guidata(hObject, handles);



% --- Executes on button press in pushbutton_calc_rib_protofildia.
function pushbutton_calc_rib_protofildia_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_calc_rib_protofildia (see GCBO)

persisLmean_nm = str2double(get(handles.set_PL_mean,'String'));
persisLSD_nm = str2double(get(handles.set_PL_SD,'String'));
PLmean = persisLmean_nm ./1E9 ; % persistence length in meters
PLsd = persisLSD_nm ./1E9 ;

kB = 1.3806503e-23 ;  %Boltzmann constant
Temp = str2double(get(handles.set_Temperature,'String'));
T = 273.15 + Temp; % Temperature in deg Kelvin
BendingR_mean = kB * T * PLmean ;  % bending rigidity in N.m^2
BendingR_SD = kB * T * PLsd ;  
set(handles.bendr_mean,'String',BendingR_mean) ;
set(handles.bendr_SD,'String',round2(BendingR_SD,1E-33)) ;

nfil = str2double(get(handles.set_number_protofil,'String'));
nfil_sd = 1;

protofil_dia = str2double(get(handles.set_protofil_dia,'String'));
protofil_dia_sd = str2double(get(handles.set_protofil_dia_SD,'String'));
protofil_dia = protofil_dia ./1E9 ; % converts in meters
protofil_dia_sd = protofil_dia_sd ./1E9 ;

% calculation of second moment of area
Irib = protofil_dia^3 .* (nfil .* protofil_dia) ./ 12 ;
set(handles.ans_ribbonI,'String',Irib) ;
set(handles.inertia_mean,'String',Irib) ;

% calculates upper bound I+
Iribplus = (protofil_dia + protofil_dia_sd)^3 .* ((nfil + nfil_sd) .* ...
    (protofil_dia + protofil_dia_sd) ) ./ 12 ;
set(handles.ans_ribbonIplus,'String',Iribplus) ; 
set(handles.inertia_SDplus,'String',abs(Irib - Iribplus)) ;


% calculates lower bound I-
a = nfil - nfil_sd;
if a ~=0
    Iribminus = (protofil_dia - protofil_dia_sd)^3 .* ((nfil - nfil_sd) .* ...
    (protofil_dia - protofil_dia_sd) ) ./ 12/2.66 ;
else
    Iribminus = (protofil_dia - protofil_dia_sd)^3 .* ( 1 .* ...
    (protofil_dia - protofil_dia_sd) ) ./ 12/2.66 ;
end
set(handles.ans_ribbonIminus,'String',Iribminus) ;
set(handles.inertia_SDminus,'String',abs(Irib - Iribminus)) ;

% calculate Young's modulus
ribbonY = BendingR_mean ./ Irib ;
%converts to MPa and round to 0.1MPa
set(handles.ans_ribbon_youngm,'String',round2(ribbonY ./1E6,0.1)) ;

ribbonEplus = (BendingR_mean + BendingR_SD) ./ Iribminus ;
set(handles.ans_ribbon_Eplus,'String',round2(ribbonEplus ./1E6,0.1)) ;
ribbonEminus = (BendingR_mean - BendingR_SD) ./ Iribplus ;
set(handles.ans_ribbon_Eminus,'String',round2(ribbonEminus ./1E6,0.1)) ;

guidata(hObject,handles);


function set_width_Callback(hObject, eventdata, handles)
% hObject    handle to set_width (see GCBO)

% --- Executes during object creation, after setting all properties.
function set_width_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function set_width_sd_Callback(hObject, eventdata, handles)
% hObject    handle to set_width_sd (see GCBO)

% --- Executes during object creation, after setting all properties.
function set_width_sd_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in btn_calcRib_HandWidth.
function btn_calcRib_HandWidth_Callback(hObject, eventdata, handles)
% calculates I only with FIB HEIGHT and FIB WIDTH info
%                   (NEITHER using Nb of protofils NOR their sizes)

persisLmean_nm = str2double(get(handles.set_PL_mean,'String'));
persisLSD_nm = str2double(get(handles.set_PL_SD,'String'));
PLmean = persisLmean_nm ./1E9 ; % persistence length in meters
PLsd = persisLSD_nm ./1E9 ;

kB = 1.3806503e-23 ;  %Boltzmann constant
Temp = str2double(get(handles.set_Temperature,'String'));
T = 273.15 + Temp; % Temperature in deg Kelvin
BendingR_mean = kB * T * PLmean ;  % bending rigidity in N.m^2
BendingR_SD = kB * T * PLsd ;  
set(handles.bendr_mean,'String',BendingR_mean) ;
set(handles.bendr_SD,'String',round2(BendingR_SD,1E-33)) ;

fibH = str2double(get(handles.set_ribfibH,'String'));
fibH_sd = str2double(get(handles.set_ribfibH_SD,'String'));
fibH = fibH ./1E9 ; % Height in meters
fibH_sd = fibH_sd ./1E9 ;

nfil = str2double(get(handles.set_number_protofil,'String'));
nfil_sd = 1;

fibW = str2double(get(handles.set_width,'String'));
fibW_sd = str2double(get(handles.set_width_sd,'String'));
fibW = fibW ./1E9 ; % converts in meters
fibW_sd = fibW_sd ./1E9 ;

% calculation of second moment of area
Irib = fibH^3 .* fibW ./ 12 ;
set(handles.ans_ribbonI,'String',Irib) ;
set(handles.inertia_mean,'String',Irib) ;

% calculates upper bound I+
Iribplus = (fibH + fibH_sd)^3 .* ...
    (fibW + fibW_sd)  ./ 12 ;
set(handles.ans_ribbonIplus,'String',Iribplus) ; 
set(handles.inertia_SDplus,'String',abs(Irib - Iribplus)) ;

% calculates lower bound I-
Iribminus = (fibH - fibH_sd)^3 .* ...
    (fibW - fibW_sd)  ./ 12/2.66 ;
set(handles.ans_ribbonIminus,'String',Iribminus) ;
set(handles.inertia_SDminus,'String',abs(Irib - Iribminus)) ;

% calculate Young's modulus
ribbonY = BendingR_mean ./ Irib ;
%converts to MPa and round to 0.1MPa
set(handles.ans_ribbon_youngm,'String',round2(ribbonY ./1E6,0.1)) ;

ribbonEplus = (BendingR_mean + BendingR_SD) ./ Iribminus ;
set(handles.ans_ribbon_Eplus,'String',round2(ribbonEplus ./1E6,0.1)) ;
ribbonEminus = (BendingR_mean - BendingR_SD) ./ Iribplus ;
set(handles.ans_ribbon_Eminus,'String',round2(ribbonEminus ./1E6,0.1)) ;
set(handles.set_helicfib_dia_mean,'String',fibH*1E9) ;
set(handles.set_helicfib_dia_SD,'String',fibH_sd*1E9) ;

guidata(hObject, handles);


% --- Executes on button press in btn_save_quickcalc.
function btn_save_quickcalc_Callback(hObject, eventdata, handles)
% hObject    handle to btn_save_quickcalc (see GCBO)

PersisL = str2double(get(handles.set_PL_mean,'String'));
PersisL_sd = str2double(get(handles.set_PL_SD,'String'));
BendR = str2double(get(handles.bendr_mean,'String'));
BendR_sd = str2double(get(handles.bendr_SD,'String'));
momentI = str2double(get(handles.inertia_mean,'String'));
momentI_sdplus = str2double(get(handles.inertia_SDplus,'String'));
momentI_sdminus = str2double(get(handles.inertia_SDminus,'String'));
YoungMod = str2double(get(handles.youngm_mean,'String'));
YoungMod_sd = str2double(get(handles.youngm_SD,'String'));


% [filename, pathname] = uiputfile('*.mat','Save Workspace Variables As');
% newfilename = fullfile(pathname, filename);
% save(newfilename, 'PersisL','PersisL_sd','BendR','BendR_sd',...
%    'momentI','momentI_sdplus','momentI_sdminus',...
%    'YoungMod','YoungMod_sd');

fileID = fopen('quickcalc_output.txt','w');
fprintf(fileID,'Persistence length: PL = %6.0f +/- %.0f nm\r\n',PersisL,PersisL_sd);
fprintf(fileID,'Bending rigidity: BR = %6.6u +/- %.6u N.m^2\r\n',BendR,BendR_sd);
fprintf(fileID,'Second moment of area: I = %6u m^4\r\n',momentI);
fprintf(fileID,'delta_I+: %6u m^4\r\ndelta_I-: %6u m^4\r\n',momentI_sdplus,momentI_sdminus);
fprintf(fileID,'Youngs modulus: E = %6.1f +/- %.1f MPa',YoungMod,YoungMod_sd);
fclose(fileID);

% message confirmation
msgbox(sprintf('The output file has been successfully written and saved.'),...
       'Confirmation Note','Help')

guidata(hObject,handles);




% --- Executes on button press in checkbox_gaussian.
function checkbox_gaussian_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_gaussian (see GCBO)
if (get(hObject,'Value') == get(hObject,'Max'))
   % Checkbox is checked-take appropriate action
   handles.checkbox_gaussian = 1;
else
   % Checkbox is not checked-take appropriate action
   handles.checkbox_gaussian = 0;
end
guidata(hObject,handles);





% --- Executes on button press in pushbut_makeHist.
function pushbut_makeHist_Callback(hObject, eventdata, handles)
% hObject    handle to pushbut_makeHist (see GCBO)

checkb = handles.checkbox_gaussian;
x = handles.ntot_fibrils;
contour_lengths = handles.struc.contour_lengths;

if checkb == 1
    
    for i = 1:x
        contourL_vec(1,i) = contour_lengths(1,i).contourL;
    end
    
    A = contourL_vec ;
    nbins = str2double(get(handles.set_histoNbins,'String'));
    uplim = str2double(get(handles.txt_hgfit_uplim,'String'));
    lowlim = str2double(get(handles.txt_hgfit_lowlim,'String'));
    
    if uplim == 0
        figure('Color',[1 1 1],'Name', 'Chain length distribution')
        hfitg(A,nbins)
        xlabel('Contour length (nm)'); %Write label for x-axis
        ylabel('Number of chains'); %Write label for y-axis
    else
        figure('Color',[1 1 1], 'Name', 'Chain length distribution')
        hfitg(A,nbins,lowlim,uplim)
        xlabel('Contour length (nm)'); %Write label for x-axis
        ylabel('Number of chains'); %Write label for y-axis
    end
    
else
    
    for i = 1:x
        contourL_vec(1,i) = contour_lengths(1,i).contourL;
    end
    
    A = contourL_vec ;
    nbins = str2double(get(handles.set_histoNbins,'String'));
    
    figure('Color',[1 1 1], 'Name', 'Chain length distribution')
    hist(A,nbins)
    xlabel('Contour length (nm)'); %Write label for x-axis
    ylabel('Number of chains'); %Write label for y-axis
end

guidata(hObject,handles);


function set_histoNbins_Callback(hObject, eventdata, handles)
% hObject    handle to set_histoNbins (see GCBO)


% --- Executes during object creation, after setting all properties.
function set_histoNbins_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




%function to emulate the h/pl in paw
function [n1,x1]=hpl(data,npts,miny,maxy)
if (nargin == 2)
  miny=min(min(data));
  maxy=max(max(data));
end
[n1,x1]=hist2(data,npts,miny,maxy);
stairs(x1,n1)
line([x1(1) x1(1)],[0 n1(1)])

function [no,xo] = hist2go(y,x,miny,maxy)
%HIST	Plot histograms.
%	HIST(Y) plots a histogram with 10 equally spaced bins between
%	the minimum and maximum values in Y, showing the distribution
%	of the elements in vector Y.
%	HIST(Y,N), where N is a scalar, uses N bins.
%	HIST(Y,X), where X is a vector, draws a histogram using the
%	bins specified in X.
%	[N,X] = HIST(...) does not draw a graph, but returns vectors
%	X and N such that BAR(X,N) is the histogram.
%
%	See also BAR.

%	J.N. Little 2-06-86
%	Revised 10-29-87, 12-29-88 LS
%	Revised 8-13-91 by cmt, 2-3-92 by ls.
%	Copyright (c) 1984-94 by The MathWorks, Inc.

if nargin == 0
	error('Requires one or two input arguments.')
end
if nargin == 1
    x = 10;
end
if min(size(y))==1, y = y(:); end
if isstr(x) | isstr(y)
	error('Input arguments must be numeric.')
end
[m,n] = size(y);
if max(size(x)) == 1
%    miny = min(min(y));
%    maxy = max(max(y));
    binwidth = (maxy - miny) ./ x;
    xx = miny + binwidth*[0:x];
    xx(length(xx)) = maxy;
%    x = xx(1:length(xx)-1) + binwidth/2
    x = xx(1:length(xx)-1);
else
	xx = x(:)';
    miny = min(min(y));
    maxy = max(max(y));
    binwidth = [diff(xx) 0];
    xx = [xx(1)-binwidth(1)/2 xx+binwidth/2];
    xx(1) = miny;
    xx(length(xx)) = maxy;
end
nbin = max(size(xx));
nn = zeros(nbin,n);
for i=2:nbin
    nn(i,:) = sum(y <= xx(i));
end
nn = nn(2:nbin,:) - nn(1:nbin-1,:);
if nargout == 0
    bar(x,nn);
else
  if min(size(y))==1, % Return row vectors if possible.
    no = nn';
    xo = x;
  else
    no = nn;
    xo = x';
  end
end

function [parout,chisq]=hfitg(x,n,miny,maxy)
global ndf;

hold off
if (nargin == 2)
  miny=min(min(x));
  maxy=max(max(x));
end
[ny,nx]=hpl(x,n,miny,maxy);
%par=[1,1,1000];
ind=find (x >= miny & x <= maxy);
par(1)=median(x(ind));
par(2)=std(x(ind));
par(3)=max(ny);
[parout,chisq]=chisq_min(par,nx,ny);
hold on
plot(nx,parout(3)*g(nx,parout(1),parout(2)),'r','LineWidth',2)
hold off
v=axis;


mean_value=['Mean: ',num2str(parout(1))];
std_value=['Sigma: ',num2str(parout(2))];
max_value=['Max: ',num2str(parout(3))];
chi2_value=['Chi2/ndf: ',num2str(chisq/(length(x)-3))];

xper=.27;yper=.05;
text(v(1)*xper+v(2)*(1-xper),v(3)*yper+v(4)*(1-yper),mean_value)
yper=yper+.05;
text(v(1)*xper+v(2)*(1-xper),v(3)*yper+v(4)*(1-yper),std_value)
yper=yper+.05;
text(v(1)*xper+v(2)*(1-xper),v(3)*yper+v(4)*(1-yper),max_value)
yper=yper+.05;
text(v(1)*xper+v(2)*(1-xper),v(3)*yper+v(4)*(1-yper),chi2_value)


xline=v(1)*(xper-.1)+v(2)*(.9-xper);
yline=v(3)*(yper-.05)+v(4)*(.95-yper);
%line([xline,xline],[v(4),yline])
%line([xline,v(2)],[yline,yline])

function [chisq] = gx(par)
global x;
global y;
global sigma;

% this calculates a gaussian (normal) distribution
%    x     - input value or array
%    mean  - mean of the gaussian
%    stdev - standard deviation of the gaussian
%    value - value of the the gaussian at the input value(s)
%            if input is array, output is an array

mean  = par(1);
stdev = par(2);
amp   = par(3);
fx = g(x,mean,stdev);
fx = amp * fx;

ind=find(y>0);
chisq = sum((y(ind) - fx(ind)).^2 ./ (sigma(ind)).^2);


function [value] = g(x,mean,stdev)

% this calculates a gaussian (normal) distribution
%    x     - input value or array
%    mean  - mean of the gaussian
%    stdev - standard deviation of the gaussian
%    value - value of the the gaussian at the input value(s)
%            if input is array, output is an array

value = (x - mean);
value = value.^2;
value = value/(stdev*stdev);
value = exp(-value/2);

function [parout,chisq]=chisq_min(par, xin, yin, sigmain)

global x;
global y;
global sigma;

clear parout;

x = xin;
y = yin;
if nargin == 4
  sigma = sigmain;
else
%  sigma = ones(size(x));
%  sigma = max(ones(size(x)), sqrt(y));
  sigma = sqrt(y);
end;

% options=foptions;
% options(14)= 500*prod(size(x));
% parout=fmins('gx',par,options);
options.MaxIter = 500 * prod(size(x));
parout = fminsearch( 'gx', par, options );
chisq = gx(parout);



function txt_hgfit_uplim_Callback(hObject, eventdata, handles)
% hObject    handle to txt_hgfit_uplim (see GCBO)

% --- Executes during object creation, after setting all properties.
function txt_hgfit_uplim_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function txt_hgfit_lowlim_Callback(hObject, eventdata, handles)
% hObject    handle to txt_hgfit_lowlim (see GCBO)

% --- Executes during object creation, after setting all properties.
function txt_hgfit_lowlim_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_surfintparam_Callback(hObject, eventdata, handles)
% hObject    handle to set_surfintparam (see GCBO)

% --- Executes during object creation, after setting all properties.
function set_surfintparam_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function set_surfintparam_devia_Callback(hObject, eventdata, handles)
% hObject    handle to set_surfintparam_devia (see GCBO)

% --- Executes during object creation, after setting all properties.
function set_surfintparam_devia_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function set_surfintparam_corel_Callback(hObject, eventdata, handles)
% hObject    handle to set_surfintparam_corel (see GCBO)

% --- Executes during object creation, after setting all properties.
function set_surfintparam_corel_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_margins_Callback(hObject, eventdata, handles)
% hObject    handle to set_margins (see GCBO)

% --- Executes during object creation, after setting all properties.
function set_margins_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_margin.
function checkbox_margin_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_margin (see GCBO)
if (get(hObject,'Value') == get(hObject,'Max'))
   % Checkbox is checked-take appropriate action
   handles.checkbox_margin = 1;
else
   % Checkbox is not checked-take appropriate action
   handles.checkbox_margin = 0;
end
guidata(hObject,handles);



function set_ellipsoid_width_mean_Callback(hObject, eventdata, handles)
% hObject    handle to set_ellipsoid_width_mean (see GCBO)

% --- Executes during object creation, after setting all properties.
function set_ellipsoid_width_mean_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function set_ellipsoid_width_SD_Callback(hObject, eventdata, handles)
% hObject    handle to set_ellipsoid_width_SD (see GCBO)

% --- Executes during object creation, after setting all properties.
function set_ellipsoid_width_SD_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbut_aboutEasy.
function pushbut_aboutEasy_Callback(hObject, eventdata, handles)
% hObject    handle to pushbut_aboutEasy (see GCBO)
msgbox(sprintf('Guillaume Lamour \n*******************\nUniversity of British Columbia, Canada\n*******************\nlamour99@hotmail.com \n*******************\nPlease cite: Lamour et al. Source Code for Biology and Medicine 2014;9:16\n*******************\n2014 \n*******************\nEnjoy'),...
       'About Easyworm','Help')



function set_Temperature_Callback(hObject, eventdata, handles)
% hObject    handle to set_Temperature (see GCBO)


% --- Executes during object creation, after setting all properties.
function set_Temperature_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_Temperature (see GCBO)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in toggle_showall.
function toggle_showall_Callback(hObject, eventdata, handles)
a=get(hObject,'Value');
if a ==1
    set(handles.uipanel_contourE2E,'Visible','On') 
    set(handles.uipanel_corel,'Visible','On')
    set(handles.uipanel_deviations,'Visible','On')
    set(handles.uipanel_inputdata,'Visible','On')
    set(handles.uipanel_contourE2E,'Visible','On')  
    set(handles.uipanel_plot,'Visible','On')
    set(handles.uipanel_quickcalc,'Visible','On') 
    set(handles.uipanel_surfparam,'Visible','On') 
    set(handles.uipanel_elasticmod,'Visible','On')  
    set(handles.uipanel_lengths,'Visible','On')
    set(handles.uipanel_saveoutputs,'Visible','On')
else
    set(handles.uipanel_contourE2E,'Visible','Off') 
    set(handles.uipanel_corel,'Visible','Off')
    set(handles.uipanel_deviations,'Visible','Off')
    set(handles.uipanel_inputdata,'Visible','On')
    set(handles.uipanel_contourE2E,'Visible','Off')  
    set(handles.uipanel_plot,'Visible','Off')
    set(handles.uipanel_quickcalc,'Visible','Off') 
    set(handles.uipanel_surfparam,'Visible','Off') 
    set(handles.uipanel_elasticmod,'Visible','Off')  
    set(handles.uipanel_lengths,'Visible','Off')
    set(handles.uipanel_saveoutputs,'Visible','Off')
end
guidata(hObject, handles)


% --- Executes on button press in toggle_output.
function toggle_output_Callback(hObject, eventdata, handles)
a=get(hObject,'Value');
if a ==1
    set(handles.uipanel_saveoutputs,'Visible','On')  
else
    set(handles.uipanel_saveoutputs,'Visible','Off')
end



% --- Executes on button press in pushbut_help.
function pushbut_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbut_help (see GCBO)
msgbox(sprintf('GETTING STARTED\n*********************\n1- In Input Data panel, press the Load Data button and load a .mat file generated by Easyworm1 (experimental chains) or by Synchains (synthetic ones).\n2- In Contour/End-to-end, Tangent Correlations, and Deviations/Secant Midpoint panels, click on the Launch Fit buttons.\n3- Click on Kurtosis in the 2D Equilibration Tests panel to determine whether your chains have equilibrated in 2D.  If they did, the kurtosis of the angle distribution is very close to three up to the persistence length, then it starts decreasing.  If indeed it is very close to three, keep on analyzing your chains with the 2D fit.  Otherwise, click on Switch to Noneq. (stands for "nonequilibrated 2D conformations").   \n4- Optimize the fitting parameters and press the Launch Fit buttons to fit the data again.  Making figures is useful to determine the optimal boundaries and bin numbers for fitting.  \n*********************\nIMPORTANT\n*********************\nUse the Deviations/Secant Midpoint panel only if the persistence length is much higher than the contour length of the chains.\n*********************\n5- In the Save Outputs panel, click on Sum & Bin All Chains.  It will compile the results and generate data for all chains.  \n6- Click on Save PL Data, which will create suitable .txt outputs from the data you just generated ("PL" means Persistence Length).\n*********************\nOTHER FUNCTIONS\n*********************\n- To plot the fibrils, go to the Plot Chains panel.  Clicking on Raw Trace will plot fibrils as they were fitted using Easyworm1.  Clicking on Smoothed Trace will plot them considering a fitting parameter (see Easyworm1 interface) multiplied by 4.  \n- Use the Chain Lengths panel to display the contour length distribution.\n- Use the Elastic Modulus panel to derive an axial elastic modulus for the chain, using the values determined for persistence length and size (height and width) of the chain.  Note that three distinct models are available for the geometry of the cross-section (CS): helical (circular CS), ellipsoidal (ellipsoidal CS) and tape (rectangular CS).   \n- Note that the values inside the boxes of the Quick Calculator panel will be updated as soon as buttons will be pressed in the Elastic Modulus panel.\n'),...
       'Help','Help')

 

% --- Executes on button press in btn_fractal_exp.
function btn_fractal_exp_Callback(hObject, eventdata, handles)
% hObject    handle to btn_fractal_exp (see GCBO)
test2Dfractal(hObject, handles)


function test2Dfractal(hObject, handles)

PLworm_mean = str2double(get(handles.ans_pworm_mean,'String'));

if PLworm_mean == 0
    msgbox(sprintf('Before pressing this button, you must collect some data regarding the end-to-end distances as a function of contour lengths. This can be accomplished by pressing the "Launch Fit" button in the "Contour / End-to-end" panel.'),...
        'Attention','Warn')
    return
else 
end

z = handles.ntot_fibrils;
mod1 = str2double(get(handles.ans_nfibrilsworm,'String'));
mod2 = str2double(get(handles.ans_nfibrilscorel,'String'));
mod3 = str2double(get(handles.ans_nfibrilsdevia,'String'));

if mod1 < z
    z = mod1;
elseif mod2 < z
    z = mod2;
elseif mod3 < z
    z = mod3;
end

numBins_wo = str2double(get(handles.bins_test2D,'String'));
maxval_wo = str2double(get(handles.set_nb_upmax_worm,'String')); 
upfit = str2double(get(handles.set_nb_upfit_worm,'String'));
if maxval_wo < upfit
   upfit = maxval_wo;
   set(handles.set_nb_upfit_worm,'String',upfit) ;
 else
 end

minval_wo = 0;

    matwo = handles.struc.mat_wormlike;
    new1 = [];
     for i = 1:z
        new1 = [new1; matwo(1,i).wormfib];
     end
     vecworm = new1;  % vecworm(,2) is the contour length of each fib
                       % vecworm(,1) is the end-to-end distance

%BIN DATA
    a = vecworm(:,2); %split into x and y
    b = vecworm(:,1);

%%%---------------------------------------------------------------------    
    botEdge = minval_wo; % define limits
    topEdge = maxval_wo; % define limits
    binEdges = linspace(botEdge, topEdge, numBins_wo+1);
    [~,whichBin] = histc(a, binEdges);
    for i = 1:numBins_wo
        flagBinMembers = (whichBin == i);
        binMembers     = b(flagBinMembers);
        %---MODIFIED EXPRESSION HERE
         binMean_E2E(i) = mean(binMembers);
        %---setup the length as the center of two separate bin edges
        length_E2E(i) = binEdges(i) + (((topEdge - botEdge)/numBins_wo)/2 );
    end
    test2Dfractal_length = length_E2E;
    test2Dfractal_binMean = binMean_E2E;

    f = figure('Color',[1 1 1], 'Name', 'Test 2D - Fractal Exponent' );
    set(f, 'Position', [800 100 700 500])
    xmax = topEdge; ymax = max(binMean_E2E);
    axis([1 xmax 1 ymax])
    loglog(test2Dfractal_length,test2Dfractal_binMean, '-s') 
    hold on
      
% Fit: linear regression to determine scaling  exponent for contour length
% >> persistence length
   datapoints_x = test2Dfractal_length;
   datapoints_y = test2Dfractal_binMean;
[xData, yData] = prepareCurveData( datapoints_x, datapoints_y );

% Set up fittype and options.
ft = fittype( 'a*x + b', 'independent', 'x', 'dependent', 'y' ); 
opts = fitoptions( ft ); opts.Display = 'Off'; opts.Lower = -Inf; opts.Upper = Inf;
upfit = str2double(get(handles.set_nb_upfit_worm,'String'));
lowfit = str2double(get(handles.set_nb_lowfit_worm,'String'));
ex = excludedata( xData, yData, 'domain', [lowfit  upfit] );
opts.Exclude = ex;

% Fit model to data.
[~, gof] = fit( xData, yData, ft, opts );
c = fit( xData, yData, ft, opts );
fractal_fit = coeffvalues(c);
fractal_gof = gof;

slope = fractal_fit(1,1);
intersect = fractal_fit(1,2);
Rsquare = fractal_gof.rsquare;
     x = lowfit:.1:upfit;
     y = slope*x + intersect;
     p = plot(x,y);
     set(p,'color','red','linewidth',2);

xlabel('Contour length (nm)')
ylabel('Mean end-to-end distance (nm)')
%axis equal

pre = ['Linear regression analysis:'];
lowfit = num2str(lowfit); upfit = num2str(upfit); 
slope = num2str(round2(slope, 0.01)); Rsquare = num2str(round2(Rsquare,0.001));
Fit_slope = ['Slope  =   ',slope, '   on interval: [', lowfit, '  ', upfit, '] (set in Outliers of Contour panel)'];
CoefDet = ['Coeff. of determination   =   ',Rsquare];
%note = ['(set in Outliers of Contour panel)']

        mTextBox = text;
        set(mTextBox,'units','normalized');
        set(mTextBox,'Position',[0.05 0.9 0]);
        set(mTextBox,'HorizontalAlignment', 'left');
        set(mTextBox,'fontsize',9);
        set(mTextBox,'color',[0 0 0 ]);
        set(mTextBox,'edgecolor',[1 1 1]);
        set(mTextBox,'String',{pre, Fit_slope, CoefDet});

Note1 = ['Note:'];
Note2 = ['For contour length >> persistence length, '];
Note3 = ['Self-avoiding random walk in 2D -> Slope  = 3/4  '];
Note4 = ['           "                    "          "    in 3D -> Slope  = 3/5  '];
            
        mTextBox = text;
        set(mTextBox,'units','normalized');
        set(mTextBox,'Position',[0.45 0.14 0]);
        set(mTextBox,'HorizontalAlignment', 'left');
        set(mTextBox,'fontsize',9);
        set(mTextBox,'color',[0 0 0 ]);
        set(mTextBox,'edgecolor',[1 1 1]);
        set(mTextBox,'String',{Note1, Note2, Note3, Note4 });
                
handles.fractal_contour = test2Dfractal_length;
handles.fractal_E2E = test2Dfractal_binMean;        
        
guidata(hObject, handles);


% --- Executes on button press in toggle_2D3Dswitcher.
function toggle_2D3Dswitcher_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of toggle_2D3Dswitcher
a=get(hObject,'Value');
if a ==1
   set(handles.toggle_2D3Dswitcher,'String','Switch to 2D eq.');
   set(handles.Info_2D3D,'String','Noneq. fit','FontWeight', 'bold'); 
else
   set(handles.toggle_2D3Dswitcher,'String','Switch to Noneq.');
   set(handles.Info_2D3D,'String','2D eq. fit', 'FontWeight', 'bold'); 
end
