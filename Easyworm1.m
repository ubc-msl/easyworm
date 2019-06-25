function varargout = Easyworm1(varargin)
% EASYWORM1 MATLAB code for Easyworm1.fig 
%      EASYWORM1, by itself, creates a new EASYWORM1 or raises the existing
%      singleton*.
%
%      H = EASYWORM1 returns the handle to a new EASYWORM1 or the handle to
%      the existing singleton*.
%
%      EASYWORM1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EASYWORM1.M with the given input arguments.
%
%      EASYWORM1('Property','Value',...) creates a new EASYWORM1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Easyworm1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Easyworm1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Easyworm1

% Last Modified by GUIDE v2.5 27-Aug-2013 11:36:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Easyworm1_OpeningFcn, ...
                   'gui_OutputFcn',  @Easyworm1_OutputFcn, ...
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


% --- Executes just before Easyworm1 is made visible.
function Easyworm1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Easyworm1 (see VARARGIN)

% setup part of the closing confirmation pop-up window (see last function)
set(handles.figure1,'CloseRequestFcn',@closeGUI);

% masks the axis of the image
axis off
% Choose default command line output for Easyworm1
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Easyworm1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Easyworm1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure

% Get default command line output from handles structure
varargout{1} = handles.output;





% --- Executes on button press in btn_browse_input.
function btn_browse_input_Callback(hObject, eventdata, handles)

% grab the file containing the height maps
% input_filename = uigetfile('*.txt')
[input_filename,PathName] = uigetfile('*','Select the height map');
set(handles.ans_filename,'String',input_filename);

% fid = fopen(input_filename)
% [input_filename,rt,n] = fopen(fid)
[~, ~, ext] = fileparts(input_filename);
handles.ext = ext;

if(input_filename ~= 0)
       
    % read the heights map 
    switch ext
        case '.txt'
            A = dlmread(input_filename, '\t', 1, 1);
        case {'.bmp','.png','.jpg','.ppm','.tif'}
            A = rgb2gray(imread(input_filename));
        otherwise
            msgbox(sprintf('This software accepts only files with following extensions: .bmp, .png, .jpg, .ppm, .tif, and .txt.  Please convert your image in one of these associated formats.'),...
                'Wrong Image File Format','Warn')
    end
    
    % get resolution and store it
    image_resol1 = size(A,1);
    handles.image_resol = image_resol1;
    set(handles.txt_resol1,'String',image_resol1);
    image_resol2 = size(A,2);
    set(handles.txt_resol2,'String',image_resol2); 
    Z = [image_resol1 image_resol2];
    tolerance = max(Z)/100;
    if abs(image_resol1 - image_resol2) > tolerance
        msgbox(sprintf('Please crop your height map so that both image sides share the same number of pixels (or at least that they do not differ by more than 1%%), then reload it.'),...
       'Tolerance threshold exceeded','Warn')
        %return
    else
    % normalize heights to make the image readable
    C = round((A-min(min(A)))./(max(max(A))-min(min(A)))*255-min(min(A)));
    
    switch ext
        case '.txt'
        image(C);
    otherwise
        image(A);
    end
    colormap(bone(255));
    axis off
    % Store filename
    handles.input_filename = input_filename;
    handles.A = A;
    handles.C = C;
    end
end
guidata(hObject,handles);






% --------------------------------------------------------------------
% Select fibril
function btn_select_Callback(h, eventdata, handles)
getpoints(h,handles);


% --------------------------------------------------------------------
% Gets points on the fibril from user input, traces the
% fibril that corresponds to these points, and fits a spline
% to the resulting points.
function getpoints(h,handles)

ext = handles.ext;

text1 = 'Add next';
set(handles.txt_confirm_add_chain,'String',text1,...
    'foregroundcolor',[1 0.5 0],'FontWeight','bold');

C = handles.C;
A = handles.A;

switch ext
    case '.txt'
        image(C);
    otherwise
        image(A);
end
colormap(bone(255));
axis off;
hold on;

% Collect points on a fibril
[X, Y] = ginput;

% Show points
plot(X, Y, '-ob');

% Vectors for storing the traced points 
count = 0;
XT = 0;
YT = 0;

% Loop through all the segments
for i=1:length(X)-1
  
   % Get the inputed segment beginning
   P = [X(i) Y(i)];
   % Segment vector to next segment point
   V = [X(i+1)-X(i), Y(i+1)-Y(i)];
   l = sqrt(V(1)^2+V(2)^2);
   
   % Go along the segment and adjust to nearest maximum
   u = l;   
   for k=0:u
      % Point on the segment
      P1= P + k/u* V;
            
      max = -1e10;
      
      % The range where we look for the maximum.
      % If the points the user selected are dense then
      % reduce the range.
      if(l < 20)
          range = 2;
      else
          range = 10;
      end
      
      % Now refine by looking perpendicularly to the line 
      for t=-range:range
          try
              P2 = P1 + [V(2), -V(1)]/l * t;
              n = round(P2(2));
              m = round(P2(1));
              
              B = (double(A(n,m)) ...
                  + double(A(n+1,m)) ...
                  + double(A(n-1,m)) ...
                  + double(A(n,m+1)) ...
                  + double(A(n,m-1)));
              
              if(B > max)
                  max = B;
                  bestn = n;
                  bestm = m;
              end
          catch err
              ;
          end
      end
      
 	    count = count + 1;
        YT(count) = bestn;
	    XT(count) = bestm;
   end
end

plot(XT, YT, 'y-');

X1=XT+1e-5*rand(size(XT));
Y1=YT+1e-5*rand(size(XT));

xy = [X1;Y1]; 
df = diff(xy.').'; 

t = cumsum([0, sqrt([1 1]*(df.*df))]); 

% get some user input that influences the sizes of the knots recovered by
% the least-square fit written below
fitting_param = str2double(get(handles.txt_fitting_param,'String'));
fp = fitting_param;

% Least square fit a spline
seg = round(t(size(t, 2))/fp);  % number of pieces
if(seg == 0)
    seg = 1;
end
cv = spap2(seg, 3, t, xy);

% Least square fit a second spline with higher distance between intervals
% for drawing purposes (bspline_norm)
seg2 = round(t(size(t, 2))/fp/4);  % number of pieces
if(seg2 == 0)
    seg2 = 1;
end
cv2 = spap2(seg2, 3, t, xy);

fnplt(cv2, 'b-')
hold on
fnplt(cv, 'r-')  % plot the fitted spline in red
hold off;

handles.cv = cv;
handles.cv2 = cv2;
handles.X = X;
handles.Y = Y;
handles.XT = XT;
handles.YT = YT;
handles.seg = seg;



% ---------------------------------------------------------------
% --- Compile spline fitting and image data (--> normalize lengths vs. pixels)
%---- Executes on button press in btn_compile_and_check.
image_width = str2double(get(handles.txt_input_width,'String'));
handles.image_width = image_width;
image_resol = handles.image_resol;

output_base = get(handles.txt_output_path,'String');
output_no = get(handles.txt_no,'String');
file_name=[output_base, output_no];

L_initial_fit = length(XT);
handles.L_initial_fit = L_initial_fit;

% normalization of pixel length 
x = cv.coefs(1,:); %pass cv.coefs values (first line) to the 'x' series
y = cv.coefs(2,:); %pass cv.coefs values (second line) to the 'y' series
x_norm = round(x*(image_width*1000)/image_resol);  % lengths are converted from microns to nanometers (*1000)
y_norm = round(y*(image_width*1000)/image_resol);

% normalization of pixel length for bspline with BS parameter SET BY USER 
a = cv2.coefs(1,:); %pass cv2.coefs values (first line) to the 'a' series
b = cv2.coefs(2,:); %pass cv2.coefs values (second line) to the 'b' series
a_norm = round(a*(image_width*1000)/image_resol);  % lengths are converted from microns to nanometers (*1000)
b_norm = round(b*(image_width*1000)/image_resol);

% pass bspline with BS parameter SET BY USER to the handles structure
bspline_norm_coord(1,:) = x_norm;
bspline_norm_coord(2,:) = y_norm;
handles.bspline_norm_coord = bspline_norm_coord;

% pass bspline with BS parameter SET BY USER to the handles structure
bspline_norm_coord2(1,:) = a_norm;
bspline_norm_coord2(2,:) = b_norm;
handles.bspline_norm_coord2 = bspline_norm_coord2;


handles.file_name = file_name;
handles.output_base = output_base;
handles.output_no = output_no;

handles.x_norm = x_norm; % used in interval/check increment function
handles.y_norm = y_norm; 


% ---------------------------------------------------------------
% Check knots spline increments over the fibril
%INTERVAL FUNCTION
x_norm = handles.x_norm;
y_norm = handles.y_norm;
seg = handles.seg;

for i=1:(seg+1)
    dist = ( x_norm(i) - x_norm(i+1) ).^2 + ( y_norm(i) - y_norm(i+1) ).^2;
    intervals(i) = sqrt(dist);
end

% average_interval between each knot over all the spline
mean_itv = mean (intervals(1,:));
mean_itv = round (mean_itv);
% standard deviation
sd_itv = std(intervals(1,:));
sd_itv = round (sd_itv);

%-------- Send the values in the box of the GUI dedicated to it
set(handles.ans_mean_itv,'String',mean_itv);
set(handles.ans_sd_itv,'String',sd_itv);
set(handles.txt_nb_knots,'string',seg+2);

handles.mean_itv = mean_itv;
handles.sd_itv = sd_itv;
handles.intervals = intervals;

guidata(h,handles);


function txt_fitting_param_Callback(hObject, eventdata, handles)
%store the contents of input1_editText as a string. if the string
%is not a number then input will be empty
input = str2num(get(hObject,'String'));
%checks to see if input is empty. if so, default input1_editText to twenty
if (isempty(input))
    set(hObject,'String','15')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_fitting_param_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --------------------------------------------------------------------
% --- user input - image side size in microns 
function txt_input_width_Callback(hObject, eventdata, handles)
% hObject    handle to txt_input_width (see GCBO)

function txt_input_width_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Setup part of filename by user input
function txt_output_path_Callback(hObject, eventdata, handles)
% hObject    handle to txt_output_path (see GCBO)

% --- Executes during object creation, after setting all properties.
function txt_output_path_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Setup fibril # to be saved by user input
function txt_no_Callback(hObject, eventdata, handles)
% hObject    handle to txt_no (see GCBO)

% --- Executes during object creation, after setting all properties.
function txt_no_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





%---- Executes on button press in btn_combocalc.
%-------Executes all three procedures at the same time
function btn_combocalc_Callback(h, eventdata, handles)
% hObject    handle to btn_combocalc (see GCBO)
combo_calc(h, handles);

function combo_calc(h, handles)

%MIDPOINT-FLUCT
%--------------------------------------------------------------------
%-------Calculate the deviation of the fibril from secant midpoints

%recover data from previous functions and pass them through this one
x_norm = handles.x_norm;
y_norm = handles.y_norm;

midpointX = 0;
midpointY = 0;
sec_length = 0;
lastpoint = length(x_norm);


for n = 2:lastpoint - 1
    for j = 1:lastpoint - n
        % coordinates of secant midpoint
        midpointX = (x_norm(j) + x_norm(j+n))/2;
        midpointY = (y_norm(j) + y_norm(j+n))/2;
        %secant length associated to previous midpoint
        sec_length = sqrt (( y_norm(j+n) - y_norm(j) ).^2 +...
            (( x_norm(j+n) - x_norm(j) ).^2 ));
        shortest_dist = 1e+10;
        for i = 1:lastpoint
            % dist is the distance between secant midpoint defined above and
            % each knot/point of the fibril spline
            dist = (( midpointY - y_norm(i) ).^2 +...
                (( midpointX - x_norm(i) ).^2));
            dist = sqrt (dist);
            % Select the closest knot of the fibril spline to the secant
            % midpoint and keeps the distance that separate them
            if (dist <= shortest_dist)
                shortest_dist = dist;
                closest_i = i;
            end
        end
        % refine shortest dist
        if closest_i ~= 1 & closest_i <= lastpoint -1
            ci = closest_i;
            Avec = sqrt ( (x_norm(ci) - midpointX).^2 + (y_norm(ci) - midpointY).^2 );
            Bvec = sqrt ( (x_norm(ci+1) - midpointX).^2 + (y_norm(ci+1) - midpointY).^2 );
            Cvec = sqrt ( (x_norm(ci+1) - x_norm(ci)).^2 + (y_norm(ci+1) - y_norm(ci)).^2 );
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
        % stores the result of each loop (one loop corresponding to one value
        % of j; ie the first of the two knots that belong to the secant, and
        % the loop being reproduced for each knot over the fibril spline)
        delta(j)= shortest_dist;
        secant(j) = sec_length;
    end
    deltas(:,n-1) = delta;
    secants(:,n-1) = secant;
end
handles.deltas = deltas;
handles.secants = secants;

% pick the 2 matrices and make one in 2 dimensions
A = deltas;
B = secants;
p = length(A);
j = 0;
for n = 0:(p - 1)
    for k = 1:(p - n)
        u = k + (n*p - j);  % jump to next column (-j in terms of rows)
        C(u, 1) = A(k, n + 1 );
        C(u, 2) = B(k, n + 1 );
    end
    j = j + n;
end
deviat_single = C;
handles.deviat_single = deviat_single;


%TANTAN-COREL
%--------------------------------------------------------------------
%-------Calculate the decay of tangent-tangent correlations over the fibril

intervals = handles.intervals;

for j = 1:lastpoint - 2
    
    for i = 1:lastpoint - j - 1
        % defines coordinates of each tangent vector along the spline
        spX = x_norm(i+1) - x_norm(i);
        spY = y_norm(i+1) - y_norm(i);
        spnextX = x_norm(i+j+1) - x_norm(i+j);
        spnextY = y_norm(i+j+1) - y_norm(i+j);
        % define tangent vectors to the spline at the i^th point
        tan_i = [ spX, spY ];
        tan_k = [ spnextX, spnextY ]; %k is the (i+j)th point of the fibril
        % scalar product of the normalized vectors
        scal_prod = dot(tan_i, tan_k) / ( (norm(tan_i)) * (norm(tan_k)) );
        % Return the tangent-tangent correl = cos of normalized vector for each
        % knot and the corresponding length over the fibril between the 2 knots
        tantan_corel(i) = scal_prod;
        % increment the value of the contour length by adding the value of one
        % more segment at each loop
        contour_length(i) = sum(intervals(1,i:i+j-1));
        
    end
    cosines(:,j) = tantan_corel;
    contours(:,j) = contour_length;
end
handles.cosines = cosines;
handles.contours = contours;

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
handles.corel_single = corel_single;


%CONTOUR-ENDEND
%--------------------------------------------------------------------
%-------Calculate the contour length vs. the end-to-end length

for j = 1:lastpoint - 2
    for i = 1:lastpoint - j - 1
        %secant length associated to previous midpoint
        end2end_length(i) = sqrt (( y_norm(i+j) - y_norm(i) ).^2 +...
            (( x_norm(i+j) - x_norm(i) ).^2 ));
        % increment the value of the contour length by adding the value of 
        % one more segment at each loop
        contour_length(i) = sum(intervals(1,i:i+j-1));
    end
    end2end(:,j) = end2end_length;
    contour(:,j) = contour_length;
end
handles.end2end = end2end;
handles.contour = contour;

% pick the 2 matrices and make one in 2 dimensions
G = end2end;
H = contour;
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
handles.end2cont_single = end2cont_single;

% extract the maximum contour length calculated, ie. the fibril length
goodcontour = max(I(:,2));
handles.goodcontour = goodcontour ;


%--------------------------------------------------------------------
%BELOW is the code that fill in the structure with all the elements, pertaining 
% to one fibril, calculated above in this very same function

% IMPORTANT to keep this
output_no = handles.output_no;
n=str2num(output_no);

handles.mat_length_struct(n).contourL = handles.goodcontour;
handles.mat_intervals_struct(n).meanitv_fib = handles.mean_itv;
handles.mat_sd_struct(n).meansd_fib = handles.sd_itv;

handles.bsplines_struct(n).bspline = handles.bspline_norm_coord;
handles.bsplines2_struct(n).bspline = handles.bspline_norm_coord2;

handles.deviations_struct(n).deviafib = handles.deviat_single;
handles.correlations_struct(n).corelfib = handles.corel_single;
handles.wormlike_struct(n).wormfib = handles.end2cont_single;

% Increment the filename
n=str2num(output_no);
set(handles.txt_no, 'String', num2str(n+1));

% Set the string and color of text to confirm chain has been added
text1 = 'Chain added!';
set(handles.txt_confirm_add_chain,'String',text1,...
    'foregroundcolor',[0 0.5 0],'FontWeight','bold');

guidata(h, handles);


% --- Executes on button press in loadmat.
function loadmat_Callback(hObject, eventdata, handles)

struc = uiimport

handles.mat_length_struct = struc.contour_lengths;
handles.mat_intervals_struct = struc.intervals_means;
handles.mat_sd_struct = struc.intervals_sd;

handles.bsplines_struct = struc.bsplines_norm;
handles.bsplines2_struct = struc.bsplines_norm2pl;

handles.deviations_struct = struc.mat_deviations;
handles.correlations_struct = struc.mat_correlations;
handles.wormlike_struct = struc.mat_wormlike;

s = size(handles.mat_length_struct,2);
set(handles.txt_no, 'String', num2str(s+1));

guidata(hObject,handles);




% --- Executes on button press in increment_fibril_nb.
function increment_fibril_nb_Callback(hObject, eventdata, handles)
n = str2num(get(handles.txt_no,'String'));
set(handles.txt_no, 'String', num2str(n+1));

guidata(hObject, handles);

% --- Executes on button press in pushminus.
function pushminus_Callback(hObject, eventdata, handles)
n = str2num(get(handles.txt_no,'String'));
set(handles.txt_no, 'String', num2str(n-1));

guidata(hObject, handles);



% --- Executes on button press in btn_savedata.
function btn_savedata_Callback(h, eventdata, handles)

% Important Stuff
contour_lengths = handles.mat_length_struct;
intervals_means = handles.mat_intervals_struct;
intervals_sd = handles.mat_sd_struct;

bsplines_norm = handles.bsplines_struct;
bsplines_norm2pl = handles.bsplines2_struct;

mat_deviations = handles.deviations_struct;
mat_correlations = handles.correlations_struct;
mat_wormlike = handles.wormlike_struct;

%output_base = handles.output_base;
output_base = get(handles.txt_output_path,'String');
output_sup = '_ew1';
updated_filename=[output_base, output_sup];

%[updated_filename, pathname] = uiputfile('*.mat','Save Easyworm1 data As');
%newfilename = fullfile(pathname, updated_filename);
save(updated_filename,'contour_lengths','intervals_means',...
    'intervals_sd','bsplines_norm','bsplines_norm2pl','mat_deviations',...
    'mat_correlations','mat_wormlike');

% message confirmation
msgbox(sprintf('The data has been successfully saved. You can now use Easyworm2 to load the .mat file.'),...
       'Confirmation Note','Help')


guidata(h, handles);



function closeGUI(src,evnt)
%this function is called when the user attempts to close the GUI window

%this command brings up the close window dialog
selection = questdlg('Do you want to close Easyworm1? Note that you can save a .mat file and reload it later, using the ''Reload Data'' button.',...
                     'Close Request Function',...
                     'Yes','No','No');

%if user clicks yes, the GUI window is closed, otherwise
%nothing happens
switch selection,
   case 'Yes',
    delete(gcf)
   case 'No'
     return
end


% --- Executes during object creation, after setting all properties.
function txt_nb_knots_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_nb_knots (see GCBO)


% --- Executes on button press in pushbut_aboutEasy.
function pushbut_aboutEasy_Callback(hObject, eventdata, handles)
% hObject    handle to pushbut_aboutEasy (see GCBO)
msgbox(sprintf('Guillaume Lamour \n*******************\nUniversity of British Columbia, Canada\n*******************\nlamour99@hotmail.com \n*******************\nPlease cite: Lamour et al. Source Code for Biology and Medicine 2014;9:16\n*******************\n2014 \n*******************\nEnjoy'),...
       'About Easyworm','Help')


% --- Executes on button press in help.
function help_Callback(hObject, eventdata, handles)
msgbox(sprintf('Step-by-step instructions on how to use Easyworm are available online at: http://www.scfbm.org/content/9/1/16. \n\n*********************\nIMPORTANT\n*********************\n\b The height maps must be in the same folder than the Easyworm1 executable in order to be loaded. \n\b Be very careful if you load an image file.  Loaded images will be converted to black and white scale. The height map is generated from the brightness of the pixels, with the brighest/darkest pixels corresponding to the highest/lowest points, respectively.\n\b In addition, it should be symmetric, that is, with both sides sharing the same number of pixels (a maximum difference of 1%% is tolerated).\n*********************\nGETTING STARTED\n*********************\n1- Load an image file or text file by clicking on Select Height Map.  \n2- Enter the size of the image in the appropriate textbox.\n3- Click on Select Chain and use your mouse to select points along one chain. Press Enter on your computer keyboard when you have enough points.  \n4- Modify the fitting parameter and/or repeat the previous step if you are not satisfied with the spline fit (red line). \n5- Once you are ok with the fit, click Add Chain.  You can now add other chains of the same image and/or load other images.  \n6- When all the chains of your sample have been added, click on Save Data to generate a .mat file that compiles and records the results. It is readable by Easyworm2.'),...
    'Help','Help')
