function varargout = Taste_GUI(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Taste_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Taste_GUI_OutputFcn, ...
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

% --- Executes just before Taste_GUI is made visible.
function Taste_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
handles.savedplot = []; % for plotting in selected
handles.t_step = 0.44;

handles.colorvalue = cell(17,1); % Color for plot
handles.colorvalue{1,1} =  [0.5,0,0];
handles.colorvalue{2,1} =  [1,0,0.5];
handles.colorvalue{3,1} =  [1,0.5,0.5];
handles.colorvalue{4,1} =  [1,0.5,1];
handles.colorvalue{5,1} =  [1,0,0];
handles.colorvalue{6,1} =  [1,1,0];
handles.colorvalue{7,1} =  [1,0.5,0];
handles.colorvalue{8,1} =  [0.5,0,1];
handles.colorvalue{9,1} =  [0.5,1,0];
handles.colorvalue{10,1} = [0,0.5,0];
handles.colorvalue{11,1} = [0,0,1];
handles.colorvalue{12,1} = [0,0.5,0.5];
handles.colorvalue{13,1} = [0,1,1];
handles.colorvalue{14,1} = [1,0,1];
handles.colorvalue{15,1} = [0.5,0.5,0];
handles.colorvalue{16,1} = [0,0,0.5];
handles.colorvalue{17,1} = [0.5,0,0.5];

guidata(hObject, handles);% Update handles structure

% UIWAIT makes main wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = Taste_GUI_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

% --- Executes on button press in CIRCLE.
function CIRCLE_Callback(hObject, eventdata, handles)

handles.ratio = [];
handles.gcamp = [];
handles.tdtomato = [];
handles.time = [0:handles.t_step:handles.t_step*(handles.N_img-1)];

axes(handles.axes2)
imshow(handles.ref_img./511); %imshow(evalin('base','ref_img./511'));
h = imellipse; % Select the elliptical ROI
handles.mask = createMask(h); % Binary mask for ROI
guidata(hObject, handles);

N_img = handles.N_img;

for i=1:N_img
    tmp_gcamp = mean(mean(handles.data_reg(:,:,2,i).*handles.mask));
    tmp_tdtomato = mean(mean(handles.data_reg(:,:,3,i).*handles.mask));;
    handles.gcamp = [handles.gcamp tmp_gcamp];
    handles.tdtomato = [handles.tdtomato tmp_tdtomato]; 
end

handles.ratio = handles.gcamp./handles.tdtomato;
tmp = handles.ratio;
tmp = mean(tmp(1:100));
handles.ratio = handles.ratio./tmp-1;

axes(handles.axes3)
cla;
plot(handles.time,handles.gcamp,'g');
hold on;
plot(handles.time,handles.tdtomato,'r');
xlabel('Time (s)','FontSize',11);
ylabel('Intensity (a.u.)','FontSize',11);
title('GCaMP(green) / tdTomato(red)','FontSize',11); 
axis tight;
line([20 20],get(gca,'ylim'));
line([26 26],get(gca,'ylim'));
grid on;

axes(handles.axes4)
cla;
plot(handles.time,handles.ratio*100,'m');
hold on
line([20 20],get(gca,'ylim'));
line([26 26],get(gca,'ylim'));
ylim([min(handles.ratio) max(handles.ratio)]);
xlabel('Time (s)','FontSize',11);
ylabel('dF/F (%)','FontSize',11);
title('Calcium (ratiometric)','FontSize',11); 
axis tight;
grid on;

guidata(hObject,handles);

% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.savedplot = [handles.savedplot; handles.ratio];
guidata(hObject,handles);

axes(handles.selected)
cla;
hold on;
for i=1:size(handles.savedplot,1)
    plot(handles.time,15-i+10*handles.savedplot(i,:),'Color',handles.colorvalue{i,1});
    axis([min(handles.time) max(handles.time) 0 16])
    ylabel('dF/F (%)')
end

axes(handles.axes2)
imshow(rgb2gray(handles.ref_img./511));
handles.mask = bwperim(handles.mask);

ind = size(handles.savedplot,1);
[yy xx zz] = size(handles.ref_img);

for i=1:yy
    for j=1:xx
        if handles.mask(i,j)~= 0
            handles.colormask(i,j,:) = handles.colorvalue{ind,1};
        end
    end
end

hold on
h = imshow(handles.colormask);
hold off
set(h,'AlphaData',mean(handles.colormask,3)>0);


guidata(hObject,handles);




% --- Executes on button press in Delete.
function Delete_Callback(hObject, eventdata, handles)
% hObject    handle to Delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if size(handles.savedplot,1)==0
    axes(handles.selected)
    plot([]);
else
    tmp = handles.savedplot(1:(size(handles.savedplot,1)-1),:);
    handles.savedplot = tmp;
    guidata(hObject,handles);
end

axes(handles.selected)
cla;
hold on;
for i=1:size(handles.savedplot,1)
    plot(handles.time,15-i+10*handles.savedplot(i,:),'Color',handles.colorvalue{i,1});
    axis([min(handles.time) max(handles.time) 0 16])
    ylabel('dF/F (%)')
end

axes(handles.axes2)
imshow(rgb2gray(handles.ref_img./511));

ind = size(handles.savedplot,1)+1; % index to be deleted in colormask
[yy xx zz] = size(handles.ref_img);

tmp = handles.colorvalue;

for i=1:yy
    for j=1:xx
        if squeeze(handles.colormask(i,j,:))' == tmp{ind,1}
            handles.colormask(i,j,:) = [0 0 0];
        end
    end
end

hold on
h = imshow(handles.colormask);
hold off
set(h,'AlphaData',mean(handles.colormask,3)>0);


guidata(hObject,handles);



function Filename_Callback(hObject, eventdata, handles)
% hObject    handle to Filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Filename as text
%        str2double(get(hObject,'String')) returns contents of Filename as a double

handles.fname = get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function Filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Finish.
function Finish_Callback(hObject, eventdata, handles)
disp(['Saving to  ' handles.fname '.xls'])
xlswrite([handles.fname '.xls'],[handles.time; handles.savedplot]');

axes(handles.axes2);
set(gcf,'PaperPositionMode','auto') 
saveas(gcf,handles.fname,'bmp');


% --- Executes on button press in NewAnal.
function NewAnal_Callback(hObject, eventdata, handles)

handles.savedplot = [];
handles.mask = [];
handles.colormask = [];

[name path] = uigetfile('*.TIF','Select the stack');
cd(path);
N_img = length(imfinfo(name)); % Number of images
N_pixel = 512; % 512 by 512 pixels
data = zeros(N_pixel,N_pixel,3,N_img);
for i = 1:N_img
    tmp1 = imread(name,'TIF',i); 
    data(:,:,:,i) = tmp1(:,:,:); % Load the RGB images
end
clear tmp;

% Crop the CIRCLE
disp('Crop the region of interest')
axes(handles.axes2)
imshow(data(:,:,:,1)./511);
h = imrect; % Select the rectangular ROI
position = wait(h);
data_crop = zeros([size(imcrop(data(:,:,:,1),position)) N_img]);
for i=1:N_img
    data_crop(:,:,:,i) = imcrop(data(:,:,:,i),position);
end
axes(handles.axes2)
imshow(data_crop(:,:,:,1)./511);
title('Cropped ROI')

% Image registration
% Blurring the image with 2D Gaussian filtering to reduce speckle noise
h = fspecial('gaussian',[4 4], 1.5); % Gaussian filter
data_reg = zeros(size(data_crop));
for i=1:N_img
    data_reg(:,:,:,i) = imfilter(imfilter(data_crop(:,:,:,i),h),h); % image blurred
end
clear h;

% Set parameters for image registration
[optimizer metric] = imregconfig('monomodal');
optimizer.GradientMagnitudeTolerance = 5.0e-04;
optimizer.MinimumStepLength = 5.0e-05;
optimizer.MaximumStepLength = 6.25e-02;
optimizer.MaximumIterations = 100;
optimizer.RelaxationFactor = 0.5;

% Registration
% Use blue channel to get the transform (tform)
% and apply the 'tform' to both green and blue channel (imwarp)
ref_img=zeros(size(data_reg(:,:,3,i))); % Reference image
for i = 1:30
    ref_img = ref_img + data_reg(:,:,3,i);
end
ref_img = ref_img/30; % Averaged image (frame 1-30)

for i = 2:N_img
    tform = imregtform(data_reg(:,:,3,i), ref_img, 'translation', optimizer, metric); % Get the transform
    data_reg(:,:,3,i) = imwarp(data_reg(:,:,3,i), tform, 'OutputView', imref2d(size(ref_img)));
    data_reg(:,:,2,i) = imwarp(data_reg(:,:,2,i), tform, 'OutputView', imref2d(size(ref_img)));
    %data_reg(:,:,1,i) = imwarp(data_reg(:,:,1,i), tform, 'OutputView',
    %imref2d(size(ref_img))); % Red alignment (optional)
    if i == 2
        disp('Registration in progress');
    elseif i == N_img
        disp('Registration completed');
    end
end

% Save the registered images (Optional)
imwrite(data_reg(:,:,:,1)./511,[name((1:length(name)-4)) '_reg.tif']);
for i=2:N_img
    imwrite(data_reg(:,:,:,i)./511,[name((1:length(name)-4)) '_reg.tif'],'writemode','append');
end

handles.data_reg = data_reg;
handles.ref_img = cat(3, ref_img, ref_img, ref_img);
handles.N_img = N_img;
handles.time = [0:handles.t_step:handles.t_step*(handles.N_img-1)];
handles.colormask = zeros(size(handles.ref_img));
guidata(hObject,handles);

axes(handles.axes2)
imshow(handles.ref_img./511);
axes(handles.axes3)
cla;
axes(handles.axes4)
cla;
axes(handles.selected)
cla;

% --- Executes on button press in Load.
function Load_Callback(hObject, eventdata, handles)
handles.savedplot = [];
handles.mask = [];
handles.colormask = [];

[name path] = uigetfile('*_reg.TIF','Select the registered stack');
cd(path);
tmp = imfinfo(name);
ypixel = tmp.Height;
xpixel = tmp.Width;
handles.N_img = length(tmp); % Number of images
handles.data_reg = zeros(ypixel,xpixel,3,handles.N_img);

for i=1:handles.N_img
    handles.data_reg(:,:,:,i) = imread(name,'TIF',i);
end
tmp = squeeze(handles.data_reg(:,:,2,:));
handles.ref_img = mean(tmp(:,:,[1:100]),3);
handles.ref_img = cat(3, handles.ref_img,handles.ref_img,handles.ref_img);
handles.time = [0:handles.t_step:handles.t_step*(handles.N_img-1)];

axes(handles.axes2)
imshow(handles.ref_img./511);
axes(handles.axes3)
cla;
axes(handles.axes4)
cla;
axes(handles.selected)
cla;

handles.colormask = zeros(ypixel,xpixel,3);
guidata(hObject,handles);


function framerate_Callback(hObject, eventdata, handles)
handles.t_step = str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function framerate_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in POLYGON.
function POLYGON_Callback(hObject, eventdata, handles)
% hObject    handle to POLYGON (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.ratio = [];
handles.gcamp = [];
handles.tdtomato = [];
handles.time = [0:handles.t_step:handles.t_step*(handles.N_img-1)];

axes(handles.axes2)
imshow(handles.ref_img./511); %imshow(evalin('base','ref_img./511'));
h = impoly; % Select the elliptical ROI
handles.mask = createMask(h); % Binary mask for ROI
guidata(hObject, handles);

N_img = handles.N_img;

for i=1:N_img
    tmp_gcamp = mean(mean(handles.data_reg(:,:,2,i).*handles.mask));
    tmp_tdtomato = mean(mean(handles.data_reg(:,:,3,i).*handles.mask));;
    handles.gcamp = [handles.gcamp tmp_gcamp];
    handles.tdtomato = [handles.tdtomato tmp_tdtomato]; 
end

handles.ratio = handles.gcamp./handles.tdtomato;
tmp = handles.ratio;
tmp = mean(tmp(1:100));
handles.ratio = handles.ratio./tmp-1;

axes(handles.axes3)
cla;
plot(handles.time,handles.gcamp,'g');
hold on;
plot(handles.time,handles.tdtomato,'r');
xlabel('Time (s)','FontSize',11);
ylabel('Intensity (a.u.)','FontSize',11);
title('GCaMP(green) / tdTomato(red)','FontSize',11); 
axis tight;
line([20 20],get(gca,'ylim'));
line([26 26],get(gca,'ylim'));
grid on;

axes(handles.axes4)
cla;
plot(handles.time,handles.ratio*100,'m');
hold on
line([20 20],get(gca,'ylim'));
line([26 26],get(gca,'ylim'));
ylim([min(handles.ratio) max(handles.ratio)]);
xlabel('Time (s)','FontSize',11);
ylabel('dF/F (%)','FontSize',11);
title('Calcium (ratiometric)','FontSize',11); 
axis tight;
grid on;

guidata(hObject,handles);
