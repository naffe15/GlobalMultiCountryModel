global xGUI ;

%% initialize
xGUI.left_path = pwd ;
if ispc, xGUI.right_path = 'Z:\Global_Estimated_Model\GM' ; %PC
else  xGUI.right_path = 'Z:/Global_Estimated_Model/GM' ; %MAC
end
xGUI.current_path = pwd ;
xGUI.search_pattern = '*.fig' ;

% temp get list of pictures
yourcell = get_file_names(xGUI.current_path, xGUI.search_pattern) ;

% yourcell = {} ;
% a = dir([xGUI.left_path '\*.fig']) ;
% for c = 1:size(a, 1)
%     yourcell{c} = a(c).name
% end


%% Produce the window
xGUI.main = dialog('units','normalized','outerposition',[.4 .0 .2 .97],'Name','<<HERE                    OTHER>>') ;

xGUI.pattern = uicontrol(xGUI.main ,'Style', 'edit','Units','normalized','Position',[.05 .19 .9 .03],...
    'string','*.fig','Callback',@F_pattern);

xGUI.listbox = uicontrol(xGUI.main ,'Style', 'listbox','Units','normalized','Position',[.05 .24 .9 .75],...
    'string',yourcell,'Callback',@F_listbox);

xGUI.left_button = uicontrol(xGUI.main,'Style', 'push', ...
    'String', ['<html> <b>(((-- LEFT (HERE) ------------------------------------------</b> <br />' make_html(xGUI.left_path) '</html>'], ...
    'Units','normalized',...
    'Position', [.05 .12 .9 .07], ...
    'BackgroundColor',[.6 .7 .9], ...
    'CallBack', @F_left_button);

xGUI.right_button = uicontrol(xGUI.main,'Style', 'push', ...
    'String', ['<html> <b>                  ------------------------------------------RIGHT (OTHER) --)))</b> <br />' make_html(xGUI.right_path) '</html>'], ...
    'Units','normalized',...
    'Position', [.05 .05 .9 .07] , ...
    'BackgroundColor',[.6 .7 .9], ...
    'CallBack', @F_right_button);


%---------------add your function---------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F_window(hObj, event)
event
if event == 'CloseRequestFcn'
    try 
            close( xGUI.pix_left); close( xGUI.pix_right)
    catch; end
end
end

function F_pattern(hObj,event)
global xGUI ;
listbox = xGUI.listbox ; % handle to other contorl (the list)
xGUI.search_pattern = get(hObj, 'String');

% get new names depending on pattern selected
yourcell = get_file_names(xGUI.current_path, xGUI.search_pattern) ;

% update the other control (list of pix)
set(listbox, 'Value' , 1) ;
set(listbox, 'string' , yourcell) ;
xGUI.current_file_names = get(listbox, 'String') ;
%     drawnow

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = F_listbox(hObj,event)
global xGUI ;
value=get(hObj,'value') ;
% close previous pix ONLY, if any
try
    close( xGUI.pix_left)
catch; end
try
    close( xGUI.pix_right)
catch
end

% load left
try
    pix_left = openfig([xGUI.left_path filesep xGUI.current_file_names{value}]) ;
    set(pix_left,'Units','normalized')
    set(pix_left,'Position',[.0 .04 .40 .91])
    xGUI.pix_left = pix_left ;
    set(xGUI.left_button, 'BackgroundColor',[.6 .7 .9])
catch  set(xGUI.left_button, 'BackgroundColor',[.9 .7 .6])
end
% load right
try
    pix_right = openfig([xGUI.right_path filesep xGUI.current_file_names{value}]) ;
    set(pix_right,'Units','normalized')
    set(pix_right,'Position',[.60 .04 .4 .91])
    xGUI.pix_right = pix_right ;
    set(xGUI.right_button, 'BackgroundColor',[.6 .7 .9])
catch set(xGUI.right_button, 'BackgroundColor',[.9 .7 .6])
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outString = make_html(inString)
try
    outString = [inString(1:31) strrep(inString(32:end), filesep, [filesep '<br>'])] ;
catch outString = inString ; end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F_left_button(hObj,~)
global xGUI ;
temp = uigetdir(xGUI.current_path) ;
if exist(string(temp))
    xGUI.left_path = temp ;
    xGUI.current_path = xGUI.left_path ;
    set(hObj,     'String', ['<html> <b>(((-- LEFT (HERE) ------------------------------------------</b> <br />' ...
        make_html(xGUI.left_path) '</html>']) ;
    F_pattern(xGUI.pattern, 1) ;
    
    set(xGUI.right_button, 'BackgroundColor',[.6 .7 .9])
    set(xGUI.left_button, 'BackgroundColor',[.8 .95 .95])
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F_right_button(hObj,event)
global xGUI ;
temp = uigetdir(xGUI.current_path) ;
if exist(string(temp))
    xGUI.right_path = temp ;
    xGUI.current_path = xGUI.right_path ;
    set(hObj,     'String', ['<html> <b>------------------------------------------  RIGHT (OTHER) --)))</b> <br />' ...
        make_html(xGUI.right_path) '</html>']) ;
    F_pattern(xGUI.pattern, 1) ;   %update list of files
    
    set(xGUI.left_button, 'BackgroundColor',[.6 .7 .9])
    set(xGUI.right_button, 'BackgroundColor',[.8 .95 .95])
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function files_list = get_file_names(inDir, inPattern)
% temp get list of pictures
files_list = {} ;
a = dir([inDir filesep inPattern]) ;
for c = 1:size(a, 1)
    files_list{c} = a(c).name ;
end
end
