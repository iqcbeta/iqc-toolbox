function varargout = guiblk_ver_update(varargin)
% GUIBLK_VER_UPDATE MATLAB code for guiblk_ver_update.fig
%      GUIBLK_VER_UPDATE, by itself, creates a new GUIBLK_VER_UPDATE or
%      raises the existing singleton*.
%
%      H = GUIBLK_VER_UPDATE returns the handle to a new gui_blk_ver_update or
%      the handle to the existing singleton*.
%
%      GUIBLK_VER_UPDATE('CALLBACK',hObject,eventData,handles,...) calls
%      the local function named CALLBACK in GUIBLK_VER_UPDATE.M with the
%      given input arguments.
%
%      GUIBLK_VER_UPDATE('Property','Value',...) creates a new gui_blk_ver_update
%      or raises the existing singleton*.  Starting from the left, property
%      value pairs are applied to the GUI before gui_blk_ver_update_OpeningFcn gets
%      called.  An unrecognized property name or invalid value makes
%      property application stop.
%      All inputs are passed to gui_blk_ver_update_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_blk_ver_update

% Last Modified by GUIDE v2.5 21-Aug-2012 15:47:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gui_blk_ver_update_OpeningFcn, ...
    'gui_OutputFcn',  @gui_blk_ver_update_OutputFcn, ...
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


% --- Executes just before gui_blk_ver_update is made visible.
function gui_blk_ver_update_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_blk_ver_update (see VARARGIN)


if ~isempty(varargin)
    switch length(varargin)
        case 1
            o_file=which(varargin{1});
            o_file_v=exist(o_file);
            n_file=o_file;
            if o_file_v~=4
                error('the file class not mdl')
            end
            if size(o_file,1)>1
                error('the file name not only one')
            end
            temp_file=o_file;
            temp_file=[temp_file(1:length(temp_file)-4),'_temp.mdl'];
            copyfile(o_file,temp_file);
            temp_filename=temp_file;
        case 2
            o_file=which(varargin{1});
            o_file_v=exist(o_file);
            if o_file_v~=4
                error('the file class not mdl')
            end
            if size(o_file,1)>1
                error('the file name not only one')
            end
            n_file=varargin{2};
            copyfile(o_file,n_file)
            temp_filename = n_file;
        otherwise
            error('error input number')
    end
    
    
    computer_info=computer;
    if ~isempty(strfind(computer_info,'PCWIN'))
        dot_pos=strfind(o_file,'\');
    else
        dot_pos=strfind(o_file,'/');
    end
    % dot_pos = strfind(temp_filename,'\');
    o_file_name=o_file(dot_pos(length(dot_pos))+1:...
        length(o_file)-4);
    
    load_system(o_file);
    l = which('iqc_lib');
    if isempty(l)
        error('don''t find iqc_lib')
    end
    re_blk = 1;
    load_system(l);
    allblk = find_system('iqc_lib');
    for i1=2:length(allblk)
        if length(strfind(allblk{i1},'/'))<2
            re_blk = [re_blk i1];
        end
    end
    allblk(re_blk)='';
    clear re_blk
    
    AllSblk = find_system(['',o_file_name,'']);
    AllSblk(1)='';
    listbox1dpy='';
    for i1=1:length(AllSblk)
        mask_type = get_param(AllSblk{i1},'MaskType');
        if ~isempty(strfind(mask_type,'IQC'))
            iqc_name = get_param(AllSblk{i1},'name');
            iqc_ver =  mask_type(strfind(mask_type,'(ver:')+5:...
                strfind(mask_type,')')-1);
            crtver = getcrtver(allblk,...
                mask_type(1:strfind(mask_type,'(ver:')-1));
            if isempty(crtver)
                crtver=0;
            end
            fprintf('name: %s ,version: %s ,current version: %s \n',iqc_name,iqc_ver,crtver);
            listbox1dpy=[listbox1dpy;{iqc_name,iqc_ver,crtver}];
        end
    end
    
    
    % computer_info=computer;
    if ~isempty(strfind(computer_info,'PCWIN'))
        dot_pos=strfind(temp_filename,'\');
    else
        dot_pos=strfind(temp_filename,'/');
    end
    % dot_pos = strfind(temp_filename,'\');
    temp_filename = temp_filename(dot_pos(length(dot_pos))+1:...
        length(temp_filename)-4);
    
    new_filename = n_file;
    if ~isempty(strfind(computer_info,'PCWIN'))
        dot_pos=strfind(new_filename,'\');
    else
        dot_pos=strfind(new_filename,'/');
    end
    % dot_pos = strfind(new_filename,'\');
    new_filename = new_filename(dot_pos(length(dot_pos))+1:...
        length(new_filename)-4);
    
    disp('Copy File...OK!');
    disp('Replace Block...');
    
    
    load_system(o_file);
    all_blk_data = listbox1dpy;
    
    load_system(temp_filename);
    load_system('iqc_lib');
    
    for i1=1:size(all_blk_data,1)
        if strcmp(all_blk_data{i1,3},'don''t find')
            %         handles.str{handles.sc} = ['#',num2str(i1),...
            %             ' block error...Ignore'];
            %         handles.sc = handles.sc+1;
            %         set(handles.edit4,'String',handles.str);
            disp(['#',num2str(i1),' block error...Ignore']);
            
            continue
        end
        
%         if str2num(all_blk_data{i1,2})<str2num(all_blk_data{i1,3})
%             % 版本是否較舊，配合下方else 
            blk_path = [temp_filename,'/',all_blk_data{i1,1}];
            par = get_param(blk_path,'DialogParameters');
            par = fieldnames(par);
            for i2=1:size(par,1)
                par_val = get_param(blk_path,par{i2,1});
                par{i2,2} = par_val;
            end
            par{size(par,1)+1,1}='name';
            par{size(par,1),2}=get_param(blk_path,'name');
            o_mask=get_param(blk_path,'Masktype');
            n_mask=o_mask;
            n_mask=[n_mask(1:strfind(n_mask,'(ver:')+4),...
                num2str(all_blk_data{i1,3}),')'];
            n_path = find_system('iqc_lib','Masktype',n_mask);
            replace_block(new_filename,'Masktype',o_mask,n_path{1},'noprompt');
            %         new_path=find_system(new_filename,'Masktype',n_mask);
            new_path= [new_filename,'/',all_blk_data{i1,1}];
            for i2=1:size(par,1)
                set_param(new_path,par{i2,1},par{i2,2});
            end
            %         handles.str{handles.sc} = ['#',num2str(i1),...
            %             ' block is replace...OK'];
            %         handles.sc = handles.sc+1;
            %         set(handles.edit4,'String',handles.str);
            disp(['#',num2str(i1),' block is replace...OK']);
            
%         else
%             %         handles.str{handles.sc} = ['#',num2str(i1),...
%             %             ' block is version is new...Ignore'];
%             %         handles.sc = handles.sc+1;
%             %         set(handles.edit4,'String',handles.str);
%            disp(['#',num2str(i1),' block is version is new...Ignore']);
%         end
    end

   disp('Replace Block...OK!');
    disp('Done.');
    
   
    
    if ~isempty(temp_file)
        close_system(temp_file);
        delete(temp_file);
    end
    
    close_system('iqc_lib');
    save_system(new_filename,[],'OverwriteIfChangedOnDisk',true);
    close_system(new_filename);
    exit
else
    % Choose default command line output for gui_blk_ver_update
    handles.output = hObject;
    
    % Update handles structure
    guidata(hObject, handles);
    
    
    % UIWAIT makes gui_blk_ver_update wait for user response (see UIRESUME)
    % uiwait(handles.figure1);
    
end
% --- Outputs from this function are returned to the command line.
function varargout = gui_blk_ver_update_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% handles.str = '';
% handles.str{1} = 'Run Step 1.';
% handles.sc = 2;
% guidata(hObject,handles)
% set(handles.edit4,'String',handles.str);

handles.str = '';
handles.sc = 1;
[hObject,handles]=disp_statu('Run Step 1.',handles,hObject);

[FileName,PathName] = uigetfile('*.mdl','Select the Simulink file');
set(handles.edit2,'String',[PathName,FileName]);
listbox1dpy='';

% chk file path
filepath=which(FileName);
if isempty(filepath)
    %     handles.str{handles.sc} = 'The MATLAB search path can''t find file!';
    %     handles.sc = handles.sc+1;
    %     set(handles.edit4,'String',handles.str);
    [hObject,handles]=disp_statu('The MATLAB search path can''t find file!',handles,hObject);
else
    mdlfullpath=[PathName,FileName];
    load_system(mdlfullpath);
    l = which('iqc_lib');
    if isempty(l)
        error('don''t find iqc_lib')
    end
    re_blk = 1;
    load_system(l);
    allblk = find_system('iqc_lib');
    for i1=2:length(allblk)
        if length(strfind(allblk{i1},'/'))<2
            re_blk = [re_blk i1];
        end
    end
    allblk(re_blk)='';
    clear re_blk
    
    AllSblk = find_system(['',FileName(1:length(FileName)-4),'']);
    AllSblk(1)='';
    for i1=1:length(AllSblk)
        mask_type = get_param(AllSblk{i1},'MaskType');
        if ~isempty(strfind(mask_type,'IQC'))
            iqc_name = get_param(AllSblk{i1},'name');
            iqc_ver =  mask_type(strfind(mask_type,'(ver:')+5:...
                strfind(mask_type,')')-1);
            crtver = getcrtver(allblk,...
                mask_type(1:strfind(mask_type,'(ver:')-1));
            if isempty(crtver)
                crtver=0;
            end
            listbox1dpy=[listbox1dpy;{iqc_name,iqc_ver,crtver}];
        end
    end
    
    if isempty(listbox1dpy)
        %         handles.str{handles.sc} = 'Don''t have need to update a block!';
        %         handles.sc = handles.sc+1;
        %         set(handles.edit4,'String',handles.str);
        [hObject,handles]=disp_statu('Don''t have need to update a block!',handles,hObject);
        set(handles.pushbutton2,'Enable','off');
        set(handles.pushbutton1,'Enable','on');
        set(handles.checkbox1,'Enable','on');
        return
    end
    
    set(handles.uitable1,'Data',listbox1dpy);
    
    set(handles.pushbutton1,'Enable','off');
    %     handles.str{handles.sc} = 'Step 1...OK!';
    %     handles.str{handles.sc+1} = 'Run Step 2.';
    %     handles.sc = handles.sc+2;
    %     set(handles.edit4,'String',handles.str);
    [hObject,handles]=disp_statu('Step 1...OK!',handles,hObject);
    [hObject,handles]=disp_statu('Run Step 2.',handles,hObject);
    set(handles.pushbutton2,'Enable','on');
    set(handles.checkbox1,'Enable','on');
end
guidata(hObject,handles)
close_system(mdlfullpath);
close_system(l);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uiputfile('*.mdl');
set(handles.edit3,'String',[PathName,FileName]);
set(handles.pushbutton2,'Enable','off');

% handles.str{handles.sc} = 'Step 2...OK!';
% handles.str{handles.sc+1} = 'Push Start Transform.';
% handles.sc = handles.sc+2;
% set(handles.edit4,'String',handles.str);

[hObject,handles]=disp_statu('Step 2...OK!',handles,hObject);
[hObject,handles]=disp_statu('Push Start Transform.',handles,hObject);

set(handles.pushbutton4,'Enable','on');
guidata(hObject,handles)
% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% handles.str{handles.sc} = 'Copy File...';
% handles.sc = handles.sc+1;
% set(handles.edit4,'String',handles.str);

set(handles.checkbox1,'Enable','off');
set(handles.pushbutton2,'Enable','off');

[hObject,handles]=disp_statu('Copy File...',handles,hObject);

if ~strcmp(get(handles.edit2,'String'),get(handles.edit3,'String'))
    copyfile(get(handles.edit2,'String'),get(handles.edit3,'String'))
    temp_filename = get(handles.edit3,'String');
else
    temp_file=get(handles.edit2,'String');
    temp_file=[temp_file(1:length(temp_file)-4),'_temp.mdl'];
    copyfile(get(handles.edit2,'String'),temp_file);
    temp_filename=temp_file;
end
computer_info=computer;
if ~isempty(strfind(computer_info,'PCWIN'))
    dot_pos=strfind(temp_filename,'\');
else
    dot_pos=strfind(temp_filename,'/');
end
% dot_pos = strfind(temp_filename,'\');
temp_filename = temp_filename(dot_pos(length(dot_pos))+1:...
    length(temp_filename)-4);

new_filename = get(handles.edit3,'String');
if ~isempty(strfind(computer_info,'PCWIN'))
    dot_pos=strfind(new_filename,'\');
else
    dot_pos=strfind(new_filename,'/');
end
% dot_pos = strfind(new_filename,'\');
new_filename = new_filename(dot_pos(length(dot_pos))+1:...
    length(new_filename)-4);


% handles.str{handles.sc} = 'Copy File...OK!';
% handles.sc = handles.sc+1;
% set(handles.edit4,'String',handles.str);
[hObject,handles]=disp_statu('Copy File...OK!',handles,hObject);


% handles.str{handles.sc} = 'Replace Block...';
% handles.sc = handles.sc+1;
% set(handles.edit4,'String',handles.str);
[hObject,handles]=disp_statu('Replace Block...',handles,hObject);


load_system(get(handles.edit3,'String'));
all_blk_data = get(handles.uitable1,'Data');

load_system(temp_filename);
load_system('iqc_lib');

for i1=1:size(all_blk_data,1)
    if strcmp(all_blk_data{i1,3},'don''t find')
        %         handles.str{handles.sc} = ['#',num2str(i1),...
        %             ' block error...Ignore'];
        %         handles.sc = handles.sc+1;
        %         set(handles.edit4,'String',handles.str);
        [hObject,handles]=disp_statu(['#',num2str(i1),...
            ' block error...Ignore'],handles,hObject);
        
        continue
    end
    if str2num(all_blk_data{i1,2})<str2num(all_blk_data{i1,3})
        blk_path = [temp_filename,'/',all_blk_data{i1,1}];
        par = get_param(blk_path,'DialogParameters');
        par = fieldnames(par);
        for i2=1:size(par,1)
            par_val = get_param(blk_path,par{i2,1});
            par{i2,2} = par_val;
        end
        par{size(par,1)+1,1}='name';
        par{size(par,1),2}=get_param(blk_path,'name');
        o_mask=get_param(blk_path,'Masktype');
        n_mask=o_mask;
        n_mask=[n_mask(1:strfind(n_mask,'(ver:')+4),...
            num2str(all_blk_data{i1,3}),')'];
        n_path = find_system('iqc_lib','Masktype',n_mask);
        replace_block(new_filename,'Masktype',o_mask,n_path{1},'noprompt');
        %         new_path=find_system(new_filename,'Masktype',n_mask);
        new_path= [new_filename,'/',all_blk_data{i1,1}];
        for i2=1:size(par,1)
            set_param(new_path,par{i2,1},par{i2,2});
        end
        %         handles.str{handles.sc} = ['#',num2str(i1),...
        %             ' block is replace...OK'];
        %         handles.sc = handles.sc+1;
        %         set(handles.edit4,'String',handles.str);
        [hObject,handles]=disp_statu(['#',num2str(i1),...
            ' block is replace...OK'],handles,hObject);
        
    else
        %         handles.str{handles.sc} = ['#',num2str(i1),...
        %             ' block is version is new...Ignore'];
        %         handles.sc = handles.sc+1;
        %         set(handles.edit4,'String',handles.str);
        [hObject,handles]=disp_statu(['#',num2str(i1),...
            ' block is version is new...Ignore'],handles,hObject);
    end
end
% handles.str{handles.sc} = 'Replace Block...OK!';
% handles.str{handles.sc+1} = 'Done.';
% handles.sc = handles.sc+2;
% set(handles.edit4,'String',handles.str);
[hObject,handles]=disp_statu('Replace Block...OK!',handles,hObject);
[hObject,handles]=disp_statu('Done.',handles,hObject);

guidata(hObject,handles)

if ~isempty(temp_file)
    close_system(temp_file);
    delete(temp_file);
end

close_system('iqc_lib');
save_system(new_filename,[],'OverwriteIfChangedOnDisk',true);
close_system(new_filename);

set(handles.pushbutton4,'Enable','off');
set(handles.pushbutton1,'Enable','on');
set(handles.edit2,'String','');
set(handles.edit3,'String','');
set(handles.checkbox1,'value',0);


function crtver = getcrtver(all_blk,mask_name)

all_crt_mask = get_param(all_blk,'MaskType');

for i1=1:length(all_crt_mask)
    crt_mask = all_crt_mask{i1};
    if strcmp(crt_mask(1:strfind(crt_mask,'(ver:')-1),mask_name)
        crtver = crt_mask(strfind(crt_mask,'(ver:')+5:...
            strfind(crt_mask,')')-1);
        return
    end
end

crtver = 'don''t find';


function [hObject,handles]=disp_statu(dis_str,handles,hObject)

if handles.sc == 10
    handles.str = '';
    handles.sc = 1;
end

handles.str{handles.sc} = dis_str;
handles.sc = handles.sc+1;
set(handles.edit4,'String',handles.str);
guidata(hObject,handles)


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
switch get(handles.checkbox1,'value')
    case 1
        set(handles.pushbutton2,'Enable','off');
        set(handles.edit3,'String',get(handles.edit2,'String'));
        set(handles.pushbutton4,'Enable','on');
        [hObject,handles]=disp_statu('Step 2...OK!',handles,hObject);
        [hObject,handles]=disp_statu('Push Start Transform.',handles,hObject);
        
    case 0
        [hObject,handles]=disp_statu('Run Step 2.',handles,hObject);
        set(handles.pushbutton2,'Enable','on');
        set(handles.edit3,'String','');
        set(handles.pushbutton4,'Enable','off');
end
