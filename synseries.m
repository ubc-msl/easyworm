%Script that generates Monte-Carlo polymers
clear
P = 50; %set (in nm)
%contour = 1856;  %set (in nm)
contour = 1850;  %set (in nm)
contour_SD = 0;  %set (in nm)
checkboxMC = 1; %set to 1 only if you want normal distribution of contour length
%set to any other number will make the distribution comprised between
%intervals +/ contour_SD
segment = 5;  %set (in nm)
segment_SD = 0;  %set (in nm)
linewidth = 3;  %set (in px)
Nfib = 3000; %set
margin = 100; %set (in nm) to provide blank space in matlab figure

AAA = normrnd(contour,contour_SD,[1 200]); %don't change this
BBB = normrnd(segment,segment_SD,[1 200]);  %don't change this
fnames = {'file1', 'file2', 'file3','file4','file5'};
% fnames = {'file1', 'file2', 'file3','file4','file5','file6','file7','file8','file9','file10',...
%     'file11', 'file12', 'file13','file14','file15','file16','file17','file18','file19','file20',...
%     'file21', 'file22', 'file23','file24','file25','file26','file27','file28','file29','file30',...
%     'file31', 'file32', 'file33','file34','file35','file36','file37','file38','file39','file40',...
%     'file41', 'file42', 'file43','file44','file45','file46','file47','file48','file49','file50',...
%     'file51', 'file52', 'file53','file54','file55','file56','file57','file58','file59','file60',...
%     'file61', 'file62', 'file63','file64','file65','file66','file67','file68','file69','file70'};

    bsplines_struct = struct([]);
    bsplines_norm2pl = struct([]);
    mat_length_struct = struct([]);
    mat_intervals_struct = struct([]);
    mat_sd_struct = struct([]);
    deviations_struct = struct([]);
    correlations_struct = struct([]);
    wormlike_struct = struct([]);
    
    angle_struct = struct([]);
    
kappa= 1;
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
        
%         angledata = [];
%         angledata(1) = R(1,1);
        ibis = 4;
        for ibis = 4:(n+2)
            seg = BBB(randi(200));
             spline_x(ibis)= seg * cos(sum(R(1:(ibis-2)))) + spline_x(ibis-1);
             spline_y(ibis)= seg * sin(sum(R(1:(ibis-2)))) + spline_y(ibis-1);
             %angledata = [angledata; R(ibis-2)];
             
             %%HERE WE NEED TO WRITE A MATRIX THAT REPORTS THE ANGLE VALUE
             %%ALONG WITH ITS CORRESPONDING CONTOUR LENGTH
             
             %%AND WRITE THE SAME IN GUICALC (ALSO NEEDS TO GENERATE THE
             %%angles somehow)
             
        end
    else
    end
    % THE FOLLOWING LINES can be commented out if we dont need to print
    % figure files as outputs
     plot(spline_x,spline_y,'LineWidth',linewidth)
     xlabel('nm'); %Write label for x-axis
     ylabel('nm'); %Write label for y-axis
     axis equal
     %Provide blankappaspace around fibril extremities for plotting
     minx = min(spline_x)-margin;
     minx = round2(minx,margin);
     maxx = max(spline_x)+margin;
     maxx = round2(maxx,margin);
     miny = min(spline_y)-margin;
     miny = round2(miny,margin);
     maxy = max(spline_y)+margin;
     maxy = round2(maxy,margin);
     axis([ minx maxx miny maxy])
     %Export procedure
     print('-dtiff','-r300',fnames{kappa});
    
    
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
     
%      for n = 2:lastpoint - 1
%              for j = 1:lastpoint - n
%              midpointX = (x_norm(j) + x_norm(j+n))/2;
%              midpointY = (y_norm(j) + y_norm(j+n))/2;
%              sec_length = sqrt (( y_norm(j+n) - y_norm(j) ).^2 +...
%                  (( x_norm(j+n) - x_norm(j) ).^2 ));
%              shortest_dist = 1e+10;
%              for i = 1:lastpoint
%                  dist = (( midpointY - y_norm(i) ).^2 +...
%                      (( midpointX - x_norm(i) ).^2));
%                  dist = sqrt (dist);
%                  if (dist <= shortest_dist)
%                      shortest_dist = dist;
%                  end
%              end
%              delta(j)= shortest_dist;
%              secant(j) = sec_length;
%          end
%          deltas(:,n-1) = delta;
%          secants(:,n-1) = secant;
%      end
%      % pick the 2 matrices and make one in 2 dimensions
%      A = deltas;
%      B = secants;
%      p = length(A);
%      j = 0;
%      for n = 0:(p - 1)
%          for k = 1:(p - n)
%              u = k + (n*p - j);
%              C(u, 1) = A(k, n + 1 );
%              C(u, 2) = B(k, n + 1 );
%          end
%          j = j + n;
%      end
%      deviat_single = C;
    
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
    
    hello = kappa
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
    goodcontour = max(F(:,2));
    
    %--------------------------------------------------------------------
    %BELOW is the code that fill in the structure with all the elements, pertaining
    % to one fibril, calculated above in this very same function
    bspline_coord = [];
    bspline_coord(1,:) = x_norm;
    bspline_coord(2,:) = y_norm;
    
    bsplines_struct(kappa).bspline = bspline_coord;
    % angle_struct(kappa).ang = angledata; 
    
    mat_length_struct(kappa).contourL = goodcontour;
    mat_intervals_struct(kappa).meanitv_fib = mean_itv;
    mat_sd_struct(kappa).meansd_fib = sd_itv;
   % deviations_struct(kappa).deviafib = deviat_single;
    correlations_struct(kappa).corelfib = corel_single;
    wormlike_struct(kappa).wormfib = end2cont_single;
    
   angle_struct(kappa).ang = 0;
   deviations_struct(kappa).deviafib = 0;
   %correlations_struct(kappa).corelfib = 0;
   
    %Go to next fibril
    kappa= kappa+1;
end

contour_lengths = mat_length_struct;
intervals_means = mat_intervals_struct;
intervals_sd = mat_sd_struct;

bsplines_norm = bsplines_struct;
bsplines_norm2pl = 0;
angles = angle_struct;

mat_deviations = deviations_struct;
mat_correlations = correlations_struct;
mat_wormlike = wormlike_struct;

output_base = 'SampleXXX';
output_sup = '_synchains';
updated_filename=[output_base, output_sup];
save(updated_filename,'angles','contour_lengths','intervals_means',...
    'intervals_sd','bsplines_norm','bsplines_norm2pl','mat_deviations',...
    'mat_correlations','mat_wormlike');

clear