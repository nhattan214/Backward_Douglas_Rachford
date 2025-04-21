
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NAME: magnifyOnFigure 
%  
% AUTHOR: David Fernandez Prim (david.fernandez.prim@gmail.com) 
% 
% PURPOSE: Shows a functional zoom tool, suitable for publishing of zoomed 
% images and 2D plots 
%  
% INPUT ARGUMENTS: 
%                   figureHandle [double 1x1]: graphic handle of the target figure 
%                   axesHandle [double 1x1]: graphic handle of the target axes. 
%                    
% OUTPUT ARGUMENTS:  
%                   none 
%                    
% SINTAX: 
%           1)  magnifyOnFigure; 
%               $ Adds magnifier on the first axes of the current figure, with 
%               default behavior. 
% 
%           2)  magnifyOnFigure( figureHandle ); 
%               $ Adds magnifier on the first axes of the figure with handle 
%               'figureHandle', with default behavior.  
% 
%           3)  magnifyOnFigure( figureHandle, 'property1', value1,... ); 
%               $ Adds magnifier on the first axes of the figure with handle 
%               'figureHandle', with modified behavior.  
% 
%           4)  magnifyOnFigure( axesHandle ); 
%               $ Adds magnifier on the axes with handle 'axesHandle', with 
%               default behavior.  
% 
%           5)  magnifyOnFigure( axesHandle, 'property1', value1,... ); 
%               $ Adds magnifier on the axes with handle 'axesHandle', with 
%               modified behavior.  
% 
%  
% USAGE EXAMPLES: see script 'magnifyOnFigure_examples.m' 
%                  
% PROPERTIES:  
%        'magnifierShape':  'Shape of th magnifier ('rectangle' or 'ellipse' allowed, 'rectangle' as default) 
%        'secondaryAxesFaceColor':  ColorSpec 
%        'edgeWidth'                Color of the box surrounding the secondary  
%                                   axes, magnifier and link. Default 1 
%        'edgeColor':               Color of the box surrounding the secondary  
%                                   axes, magnifier and link. Default 'black' 
%        'displayLinkStyle':        Style of the link. 'none', 'straight' or 
%                                   'edges', with 'straight' as default. 
%        'mode':                    'manual' or 'interactive' (allowing 
%                                   adjustments through mouse/keyboard). Default 
%                                   'interactive'.                                 
%        'units'                    Units in which the position vectors are 
%                                   given. Only 'pixels' currently supported 
%        'initialPositionSecondaryAxes':    Initial position vector ([left bottom width height]) 
%                                           of secondary axes, in pixels  
%        'initialPositionMagnifier':        Initial position vector ([left bottom width height]) 
%                                           of magnifier, in pixels  
%        'secondaryAxesXLim':       Initial XLim value of the secondary axes 
%        'secondaryAxesYLim':       Initial YLim value of the secondary axes 
%        'frozenZoomAspectRatio':   Specially usefull for images, forces the use of the same zoom  
%                                   factor on both X and Y axes, in order to keep the aspect ratio  
%                                   ('on' or 'off' allowed, 'off' by default  
%  
% HOT KEYS (active if 'mode' set to 'interactive') 
%  
%       'up arrow':             Moves magnifier 1 pixel upwards 
%       'down arrow':           Moves magnifier 1 pixel downwards 
%       'left arrow':           Moves magnifier 1 pixel to the left 
%       'right arrow':          Moves magnifier 1 pixel to the right 
%       'Shift+up arrow':       Expands magnifier 10% on the Y-axis 
%       'Shift+down arrow':     Compress magnifier 10% on the Y-axis 
%       'Shift+left arrow':     Compress magnifier 10% on the X-axis 
%       'Shift+right arrow':    Expands magnifier 10% on the X-axis 
%       'Control+up arrow':     Moves secondary axes 1 pixel upwards 
%       'Control+down arrow':   Moves secondary axes 1 pixel downwards 
%       'Control+left arrow':   Moves secondary axes 1 pixel to the left 
%       'Control+right arrow':  Moves secondary axes 1 pixel to the right 
%       'Alt+up arrow':         Expands secondary axes 10% on the Y-axis 
%       'Alt+down arrow':       Compress secondary axes 10% on the Y-axis 
%       'Alt+left arrow':       Compress secondary axes 10% on the X-axis 
%       'Alt+right arrow':      Expands secondary axes 10% on the X-axis 
%       'PageUp':               Increase additional zooming factor on X-axis 
%       'PageDown':             Decrease additional zooming factor on X-axis 
%       'Shift+PageUp':         Increase additional zooming factor on Y-axis 
%       'Shift+PageDown':       Decrease additional zooming factor on Y-axis 
%       'Control+Q':            Resets the additional zooming factors to 0 
%       'Control+A':            Displays position of secondary axes and 
%                               magnifier in the command window 
%       'Mouse pointer on magnifier+left click'         Drag magnifier to any 
%                                                       direction 
%       'Mouse pointer on secondary axes+left click'    Drag secondary axes in any 
%                                                       direction 
% 
% TODO: 
%   - Use another axes copy as magnifier instead of rectangle (no ticks). 
%   - Adapt to work on 3D plots. 
%   - Add tip tool with interface description?. 
%   - Support several active instances on a given axes (object-oriented?). 
%   - Support 'delete' property to eliminate a given magnifier. 
% 
% 
% CHANGE HISTORY: 
%  
%   Version     |       Date    |   Author          |   Description 
%---------------|---------------|-------------------|--------------------------------------- 
%   1.0         |   28/11/2009  |   D. Fernandez    |   First version    
%   1.1         |   29/11/2009  |   D. Fernandez    |   Added link from magnifier to secondary axes    
%   1.2         |   30/11/2009  |   D. Fernandez    |   Keyboard support added    
%   1.3         |   01/11/2009  |   D. Fernandez    |   Properties added 
%   1.4         |   02/11/2009  |   D. Fernandez    |   Manual mode supported 
%   1.5         |   03/11/2009  |   D. Fernandez    |   New link style added ('edges') 
%   1.6         |   03/11/2009  |   D. Fernandez    |   Bug solved in display of link style 'edges' 
%   1.7         |   04/11/2009  |   D. Fernandez    |   Target axes selection added 
%   1.8         |   07/11/2009  |   D. Fernandez    |   Solved bug when any of the axes are reversed. 
%   1.9         |   08/11/2009  |   D. Fernandez    |   Adapted to work under all axes modes (tight, square, image, ...) 
%   1.10        |   08/11/2009  |   D. Fernandez    |   Added frozenZoomAspectRatio zoom mode, useful for images 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function magnifyOnFigure( varargin ) 
 
clear global appDataStruct 
global appDataStruct 
 
%Initialize 'appDataStruct' with default values 
setDefaultProperties(); 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%CHECK OUTPUT ARGUMENTS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
switch nargout 
     
    case 0 
        %Correct 
         
    otherwise 
        error('Number of output arguments not supported.'); 
end 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%CHECK INPUT ARGUMENTS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if nargin == 0  
    appDataStruct.figure.handle = gcf; 
    % Get number of axes in the same figure 
    childHandle = get(appDataStruct.figure.handle, 'Children'); 
    iAxes = find(strcmpi(get(childHandle, 'Type'), 'axes')); 
    % If no target axes specified, select the first found as mainAxes 
    appDataStruct.mainAxes.handle = childHandle( iAxes(1) ); 
elseif nargin > 0 
    if ishandle(varargin{1}) && strcmpi(get(varargin{1}, 'Type'), 'figure') 
        %Set figure handle 
        appDataStruct.figure.handle = varargin{1}; 
        % Get number of axes in the same figure 
        childHandle = get(appDataStruct.figure.handle, 'Children'); 
        iAxes = find(strcmpi(get(childHandle, 'Type'), 'axes')); 
        % If no target axes specified, select the first found as mainAxes 
        appDataStruct.mainAxes.handle = childHandle( iAxes(1) ); 
         
    elseif ishandle(varargin{1}) && strcmpi(get(varargin{1}, 'Type'), 'axes') 
        appDataStruct.mainAxes.handle = varargin{1};         
        % Get figure handle 
        parentHandle = get(appDataStruct.mainAxes.handle, 'Parent'); 
        iHandle = find(strcmpi(get(parentHandle, 'Type'), 'figure')); 
        % Figure is the parent of the axes 
        appDataStruct.figure.handle = parentHandle( iHandle(1) );         
         
    else 
        warning('Wrong figure/axes handle specified. The magnifier will be applied on the current figure.'); 
    end 
 
    if mod(nargin-1, 2) == 0 
                        
        %Check input properties 
        for i = 2:2:nargin  
            if  ~strcmpi( varargin{i}, 'frozenZoomAspectRatio' ) &&... 
                ~strcmpi( varargin{i}, 'magnifierShape' ) &&... 
                ~strcmpi( varargin{i}, 'secondaryaxesxlim' ) &&... 
                ~strcmpi( varargin{i}, 'secondaryaxesylim' ) &&...     
                ~strcmpi( varargin{i}, 'secondaryaxesfacecolor' ) &&... 
                ~strcmpi( varargin{i}, 'edgewidth' ) &&... 
                ~strcmpi( varargin{i}, 'edgecolor' ) &&... 
                ~strcmpi( varargin{i}, 'displayLinkStyle' ) &&... 
                ~strcmpi( varargin{i}, 'mode' ) &&... 
                ~strcmpi( varargin{i}, 'units' ) &&... 
                ~strcmpi( varargin{i}, 'initialpositionsecondaryaxes' ) &&... 
                ~strcmpi( varargin{i}, 'initialpositionmagnifier' ) 
                error('Illegal property specified. Please check.'); 
            end 
            if strcmpi( varargin{i}, 'frozenZoomAspectratio' ) 
                if ischar(varargin{i+1}) &&... 
                   ( strcmpi(varargin{i+1}, 'on') || strcmpi(varargin{i+1}, 'off') ) 
                    appDataStruct.global.zoomMode = lower(varargin{i+1}); 
                else 
                    warning(sprintf('Specified zoom mode not supported. Default values will be applied [%s].', appDataStruct.global.zoomMode)); 
                end 
            end 
            if strcmpi( varargin{i}, 'mode' ) 
                if ischar(varargin{i+1}) &&... 
                   ( strcmpi(varargin{i+1}, 'manual') || strcmpi(varargin{i+1}, 'interactive') ) 
                    appDataStruct.global.mode = lower(varargin{i+1}); 
                else 
                    warning(sprintf('Specified mode descriptor not supported. Default values will be applied [%s].', appDataStruct.global.mode)); 
                end 
            end 
            if strcmpi( varargin{i}, 'magnifierShape' ) 
                if ischar(varargin{i+1}) &&... 
                   ( strcmpi(varargin{i+1}, 'rectangle') || strcmpi(varargin{i+1}, 'ellipse') ) 
                    appDataStruct.magnifier.shape = lower(varargin{i+1}); 
                else 
                    warning(sprintf('Specified magnifier shape not supported. Default values will be applied [%s].', appDataStruct.magnifier.shape)); 
                end 
            end 
            if strcmpi( varargin{i}, 'displayLinkStyle' ) 
                if ischar(varargin{i+1}) &&... 
                   ( strcmpi(varargin{i+1}, 'straight') || strcmpi(varargin{i+1}, 'none') || strcmpi(varargin{i+1}, 'edges') ) 
                    if ~strcmpi(appDataStruct.magnifier.shape, 'rectangle') && strcmpi(varargin{i+1}, 'edges') 
                        warning(sprintf('Specified link style not supported. Default values will be applied for ''displayLinkStyle''[%s].', appDataStruct.link.displayLinkStyle)); 
                    else 
                        appDataStruct.link.displayLinkStyle = lower(varargin{i+1}); 
                    end 
                else 
                    warning(sprintf('Specified descriptor not supported. Default values will be applied for ''displayLink''[%s].', appDataStruct.link.displayLinkStyle)); 
                end 
            end 
            if strcmpi( varargin{i}, 'units' ) 
                if ischar(varargin{i+1}) && strcmpi(varargin{i+1}, 'pixels') 
                    appDataStruct.global.units = lower(varargin{i+1}); 
                else 
                    warning(sprintf('Specified units descriptor nbot supported. Default values will be applied [%s].', appDataStruct.global.units)); 
                end 
            end 
            if strcmpi( varargin{i}, 'edgewidth' ) 
                if length(varargin{i+1})==1 && isnumeric(varargin{i+1}) 
                    appDataStruct.global.edgeWidth = varargin{i+1}; 
                else 
                    warning(sprintf('Incorrect edge width value. Default value will be applied [%g].', appDataStruct.global.edgeWidth )) 
                end 
            end 
            if strcmpi( varargin{i}, 'edgecolor' ) 
                if ( length(varargin{i+1})==3 && isnumeric(varargin{i+1}) ) ||... 
                   ( ischar(varargin{i+1}) )      
                    appDataStruct.global.edgeColor = varargin{i+1}; 
                else 
                    warning('Incorrect edge color value. Default black will be applied.'); 
                end 
            end                       
            if strcmpi( varargin{i}, 'secondaryaxesfacecolor' ) 
                if ( length(varargin{i+1})==3 && isnumeric(varargin{i+1}) ) ||... 
                   ( ischar(varargin{i+1}) )      
                    appDataStruct.secondaryAxes.faceColor = varargin{i+1}; 
                else 
                    warning('Incorrect secondary axes face color value. Default white will be applied.'); 
                end 
            end  
             
            if strcmpi( varargin{i}, 'secondaryaxesxlim' ) 
                if ( length(varargin{i+1})==2 && isnumeric(varargin{i+1}) ) 
                    appDataStruct.secondaryAxes.initialXLim = varargin{i+1}; 
                else 
                    warning('Incorrect secondary axes XLim value. Default white will be applied.'); 
                end 
            end  
             
            if strcmpi( varargin{i}, 'secondaryaxesylim' ) 
                if ( length(varargin{i+1})==2 && isnumeric(varargin{i+1}) ) 
                    appDataStruct.secondaryAxes.initialYLim = varargin{i+1}; 
                else 
                    warning('Incorrect secondary axes YLim value. Default white will be applied.'); 
                end 
            end              
             
            if strcmpi( varargin{i}, 'initialpositionsecondaryaxes' ) 
                if length(varargin{i+1})==4 && isnumeric(varargin{i+1}) 
                    appDataStruct.secondaryAxes.initialPosition = varargin{i+1}; 
                else 
                    warning('Incorrect initial position of secondary axes. Default values will be applied.') 
                end 
            end 
            if strcmpi( varargin{i}, 'initialpositionmagnifier' ) 
                if length(varargin{i+1})==4 && isnumeric(varargin{i+1}) 
                    appDataStruct.magnifier.initialPosition = varargin{i+1}; 
                else 
                    warning('Incorrect initial position of magnifier. Default values will be applied.') 
                end                     
            end 
        end 
         
    else 
        error('Number of input arguments not supported.'); 
    end 
end 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%ENTRY POINT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Create secondary axes 
appDataStruct.secondaryAxes.handle = copyobj(appDataStruct.mainAxes.handle,appDataStruct.figure.handle); 
 
%Configure secondary axis 
set(appDataStruct.secondaryAxes.handle, 'Color',get(appDataStruct.mainAxes.handle,'Color'), 'Box','on'); 
set(appDataStruct.secondaryAxes.handle, 'FontWeight', 'bold',... 
                                        'LineWidth', appDataStruct.global.edgeWidth,... 
                                        'XColor', appDataStruct.global.edgeColor,... 
                                        'YColor', appDataStruct.global.edgeColor,... 
                                        'Color', appDataStruct.secondaryAxes.faceColor ); 
set(appDataStruct.figure.handle, 'CurrentAxes',appDataStruct.secondaryAxes.handle ); 
xlabel(''); ylabel(''); zlabel(''); title(''); 
% axis( appDataStruct.secondaryAxes.handle, 'normal'); %Ensure that secondary axes are not resizing 
set(appDataStruct.figure.handle, 'CurrentAxes',appDataStruct.mainAxes.handle ); 
 
%Default magnifier position 
appDataStruct.magnifier.initialPosition = computeMagnifierDefaultPosition(); 
 
%Default secondary axes position 
appDataStruct.secondaryAxes.initialPosition = computeSecondaryAxesDefaultPosition(); 
 
%Set initial position of secondary axes 
setSecondaryAxesPositionInPixels( appDataStruct.secondaryAxes.initialPosition ); 
 
%Set initial position of magnifier 
setMagnifierPositionInPixels( appDataStruct.magnifier.initialPosition );                             
 
%Update view limits on secondary axis 
refreshSecondaryAxisLimits(); 
    
%Update link between secondary axes and magnifier 
refreshMagnifierToSecondaryAxesLink(); 
 

 
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NAME: refreshSecondaryAxisLimits 
%  
% PURPOSE: Updates the view on the secondary axis, based on position and 
% span of magnifier, and extend of secondary axis. 
%  
% INPUT ARGUMENTS: 
%                   appDataStruct [struct 1x1]: global variable 
% OUTPUT ARGUMENTS:  
%                   change 'XLim' and 'YLim' of secondary axis (ACTION) 
%                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function refreshSecondaryAxisLimits() 
         
global appDataStruct; 
 
 
%If limits specified 
if ~(isempty(appDataStruct.secondaryAxes.initialXLim) ||... 
   isempty(appDataStruct.secondaryAxes.initialYLim))      
 
    limitsToSet = [... 
                        appDataStruct.secondaryAxes.initialXLim(1)... 
                        appDataStruct.secondaryAxes.initialXLim(2)... 
                        appDataStruct.secondaryAxes.initialYLim(1)... 
                        appDataStruct.secondaryAxes.initialYLim(2)... 
                    ];         
    axis(appDataStruct.secondaryAxes.handle, limitsToSet); 
     
    appDataStruct.secondaryAxes.initialXLim = []; 
    appDataStruct.secondaryAxes.initialYLim = []; 
else 
    %Get main axes limits, in axes units 
    mainAxesXLim = get( appDataStruct.mainAxes.handle, 'XLim' ); 
    mainAxesYLim = get( appDataStruct.mainAxes.handle, 'YLim' ); 
    mainAxesXDir = get( appDataStruct.mainAxes.handle, 'XDir' ); 
    mainAxesYDir = get( appDataStruct.mainAxes.handle, 'YDir' ); 
     
    %Get position and size of main axes in pixels 
    mainAxesPositionInPixels = getMainAxesPositionInPixels(); 
     
    %Compute Pixels-to-axes units conversion factors 
    xMainAxisPixels2UnitsFactor = determineSpan( mainAxesXLim(1), mainAxesXLim(2) )/mainAxesPositionInPixels(3); 
    yMainAxisPixels2UnitsFactor = determineSpan( mainAxesYLim(1), mainAxesYLim(2) )/mainAxesPositionInPixels(4); 
 
    %Get position and extend of magnifier, in pixels                       
    magnifierPosition = getMagnifierPositionInPixels(); %In pixels 
     
    %Relative to the lower-left corner of the axes 
    magnifierPosition(1) = magnifierPosition(1) - mainAxesPositionInPixels(1); 
    magnifierPosition(2) = magnifierPosition(2) - mainAxesPositionInPixels(2); 
     
    %Compute position and exted of magnifier, in axes units 
    magnifierPosition(3) = magnifierPosition(3) * xMainAxisPixels2UnitsFactor; 
    magnifierPosition(4) = magnifierPosition(4) * yMainAxisPixels2UnitsFactor; 
    if strcmpi(mainAxesXDir, 'normal') && strcmpi(mainAxesYDir, 'normal') 
        magnifierPosition(1) = mainAxesXLim(1) + magnifierPosition(1)*xMainAxisPixels2UnitsFactor; 
        magnifierPosition(2) = mainAxesYLim(1) + magnifierPosition(2)*yMainAxisPixels2UnitsFactor; 
    end 
    if strcmpi(mainAxesXDir, 'normal') && strcmpi(mainAxesYDir, 'reverse') 
        magnifierPosition(1) = mainAxesXLim(1) + magnifierPosition(1)*xMainAxisPixels2UnitsFactor; 
        magnifierPosition(2) = mainAxesYLim(2) - magnifierPosition(2)*yMainAxisPixels2UnitsFactor - magnifierPosition(4);    
    end 
    if strcmpi(mainAxesXDir, 'reverse') && strcmpi(mainAxesYDir, 'normal') 
        magnifierPosition(1) = mainAxesXLim(2) - magnifierPosition(1)*xMainAxisPixels2UnitsFactor - magnifierPosition(3); 
        magnifierPosition(2) = mainAxesYLim(1) + magnifierPosition(2)*yMainAxisPixels2UnitsFactor;  
    end 
    if strcmpi(mainAxesXDir, 'reverse') && strcmpi(mainAxesYDir, 'reverse') 
        magnifierPosition(1) = mainAxesXLim(2) - magnifierPosition(1)*xMainAxisPixels2UnitsFactor - magnifierPosition(3); 
        magnifierPosition(2) = mainAxesYLim(2) - magnifierPosition(2)*yMainAxisPixels2UnitsFactor - magnifierPosition(4); 
    end 
         
    secondaryAxisXlim = [magnifierPosition(1) magnifierPosition(1)+magnifierPosition(3)]; 
    secondaryAxisYlim = [magnifierPosition(2) magnifierPosition(2)+magnifierPosition(4)];  
     
    xZoom = appDataStruct.secondaryAxes.additionalZoomingFactor(1); 
    yZoom = appDataStruct.secondaryAxes.additionalZoomingFactor(2); 
 
    aux_secondaryAxisXlim(1) =  mean(secondaryAxisXlim) -... 
                            determineSpan( secondaryAxisXlim(1), mean(secondaryAxisXlim) )*(1-xZoom); 
    aux_secondaryAxisXlim(2) =  mean(secondaryAxisXlim) +... 
                            determineSpan( secondaryAxisXlim(2), mean(secondaryAxisXlim) )*(1-xZoom); 
    aux_secondaryAxisYlim(1) =  mean(secondaryAxisYlim) -... 
                            determineSpan( secondaryAxisYlim(1), mean(secondaryAxisYlim) )*(1-yZoom); 
    aux_secondaryAxisYlim(2) =  mean(secondaryAxisYlim) +... 
                            determineSpan( secondaryAxisYlim(2), mean(secondaryAxisYlim) )*(1-yZoom); 
 
    if aux_secondaryAxisXlim(1)<aux_secondaryAxisXlim(2) &&... 
       all(isfinite(aux_secondaryAxisXlim))      
        set( appDataStruct.secondaryAxes.handle, 'XLim', aux_secondaryAxisXlim ); 
    end 
    if aux_secondaryAxisYlim(1)<aux_secondaryAxisYlim(2) &&... 
       all(isfinite(aux_secondaryAxisYlim))              
        set( appDataStruct.secondaryAxes.handle, 'YLim', aux_secondaryAxisYlim ); 
    end
    
    %% add another box for secondaryAxes 
    annotation(  'rectangle',...
    'Units', 'pixels',...
    'Position', appDataStruct.secondaryAxes.handle.Position,...
    'LineWidth', appDataStruct.global.edgeWidth,...
    'LineStyle','-',...
    'EdgeColor',appDataStruct.global.edgeColor...
    );
    %% modified
 
end 
 
%Increase line width in plots on secondary axis 
childHandle = get(appDataStruct.secondaryAxes.handle, 'Children'); 
for iChild = 1:length(childHandle) 
    if strcmpi(get(childHandle(iChild), 'Type'), 'line') 
        set(childHandle(iChild), 'LineWidth', 2); 
    end 
    if strcmpi(get(childHandle(iChild), 'Type'), 'image') 
        %Do nothing for now 
    end 
end 
 
%Modify X and Y axes ticks for better display 
%Get main axes limits, in axes units 
% secondaryAxesXLim = get( appDataStruct.secondaryAxes.handle, 'XLim' ); 
% secondaryAxesYLim = get( appDataStruct.secondaryAxes.handle, 'YLim' ); 
% secondaryAxesXTick = get( appDataStruct.secondaryAxes.handle, 'XTick' ); 
% secondaryAxesYTick = get( appDataStruct.secondaryAxes.handle, 'YTick' ); 
%  
% if secondaryAxesXTick(1) ~= secondaryAxesXLim(1) 
%     secondaryAxesXTick = [secondaryAxesXLim(1) secondaryAxesXTick]; 
% end 
% if secondaryAxesXTick(end) ~= secondaryAxesXLim(2) 
%     secondaryAxesXTick = [secondaryAxesXTick secondaryAxesXLim(2)]; 
% end 
% if secondaryAxesYTick(1) ~= secondaryAxesYLim(1) 
%     secondaryAxesYTick = [secondaryAxesYLim(1) secondaryAxesYTick]; 
% end 
% if secondaryAxesYTick(end) ~= secondaryAxesYLim(2) 
%     secondaryAxesYTick = [secondaryAxesYTick secondaryAxesYLim(2)]; 
% end 
%  
% set( appDataStruct.secondaryAxes.handle, 'XTick', secondaryAxesXTick ); 
% set( appDataStruct.secondaryAxes.handle, 'YTick', secondaryAxesYTick ); 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NAME: determineSpan 
%  
% PURPOSE: Computes the distance between two real numbers on a 1D space. 
%  
% INPUT ARGUMENTS: 
%                   v1 [double 1x1]: first number 
%                   v2 [double 1x1]: second number 
% OUTPUT ARGUMENTS:  
%                   span [double 1x1]: computed span  
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function span = determineSpan( v1, v2 ) 
 
if v1>=0 && v2>=0 
    span = max(v1,v2) - min(v1,v2); 
end 
if v1>=0 && v2<0 
    span = v1 - v2; 
end 
if v1<0 && v2>=0 
    span = -v1 + v2; 
end 
if v1<0 && v2<0 
    span = max(-v1,-v2) - min(-v1,-v2); 
end 
    
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NAME: KeyPressCallback 
%  
% PURPOSE: Service routine to KeyPress event. 
%  
% INPUT ARGUMENTS: 
%                    
% OUTPUT ARGUMENTS:  
%                    
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function KeyPressCallback(src,eventdata) 
 
global appDataStruct 
 
currentCaracter = eventdata.Key; 
currentModifier = eventdata.Modifier; 
 
switch(currentCaracter) 
    case {'leftarrow'} % left arrow 
        %Move magnifier to the left 
        if isempty(currentModifier) 
            position = getMagnifierPositionInPixels(); 
            magnifierPosition(1) = position(1)-1; 
            magnifierPosition(2) = position(2); 
            magnifierPosition(3) = position(3); 
            magnifierPosition(4) = position(4); 
            setMagnifierPositionInPixels( magnifierPosition );          
        end 
        %Compress magnifier on the X axis 
        if strcmp(currentModifier, 'shift') 
            position = getMagnifierPositionInPixels(); 
            magnifierPosition(3) = position(3)*(1 - 0.1); 
            if strcmpi( appDataStruct.global.zoomMode, 'off') 
                magnifierPosition(4) = position(4); 
            else 
                %If 'freezeZoomAspectRatio' to 'on', be consistent 
                magnifierPosition(4) = position(4)*(1 - 0.1);  
            end 
            magnifierPosition(1) = position(1)-(-position(3)+magnifierPosition(3))/2; 
            magnifierPosition(2) = position(2)-(-position(4)+magnifierPosition(4))/2; 
            setMagnifierPositionInPixels( magnifierPosition ); 
        end 
        %Move secondary axes to the left         
        if strcmp(currentModifier, 'control') 
            position = getSecondaryAxesPositionInPixels(); 
            secondaryAxesPosition(1) = position(1)-1; 
            secondaryAxesPosition(2) = position(2); 
            secondaryAxesPosition(3) = position(3); 
            secondaryAxesPosition(4) = position(4); 
            setSecondaryAxesPositionInPixels( secondaryAxesPosition );         
        end 
        %Compress secondary axes on the X axis         
        if strcmp(currentModifier, 'alt') 
            position = getSecondaryAxesPositionInPixels(); 
            secondaryAxesPosition(3) = position(3)*(1 - 0.1); 
            if strcmpi( appDataStruct.global.zoomMode, 'off') 
                secondaryAxesPosition(4) = position(4); 
            else 
                %If 'freezeZoomAspectRatio' to 'on', be consistent 
                secondaryAxesPosition(4) = position(4)*(1 - 0.1);  
            end 
            secondaryAxesPosition(1) = position(1)-(-position(3)+secondaryAxesPosition(3))/2; 
            secondaryAxesPosition(2) = position(2)-(-position(4)+secondaryAxesPosition(4))/2; 
            setSecondaryAxesPositionInPixels( secondaryAxesPosition );  
        end 
         
    case {'rightarrow'} % right arrow 
         %Move magnifier to the right 
        if isempty(currentModifier) 
            position = getMagnifierPositionInPixels(); 
            magnifierPosition(1) = position(1)+1; 
            magnifierPosition(2) = position(2); 
            magnifierPosition(3) = position(3); 
            magnifierPosition(4) = position(4); 
            setMagnifierPositionInPixels( magnifierPosition );          
        end 
        %Expand magnifier on the X axis 
        if strcmp(currentModifier, 'shift') 
            position = getMagnifierPositionInPixels(); 
            magnifierPosition(3) = position(3)*(1 + 0.1); 
            if strcmpi( appDataStruct.global.zoomMode, 'off') 
                magnifierPosition(4) = position(4); 
            else 
                %If 'freezeZoomAspectRatio' to 'on', be consistent 
                magnifierPosition(4) = position(4)*(1 + 0.1);  
            end 
            magnifierPosition(1) = position(1)-(-position(3)+magnifierPosition(3))/2; 
            magnifierPosition(2) = position(2)-(-position(4)+magnifierPosition(4))/2; 
            setMagnifierPositionInPixels( magnifierPosition ); 
        end 
        %Move secondary axes to the right         
        if strcmp(currentModifier, 'control') 
            position = getSecondaryAxesPositionInPixels(); 
            secondaryAxesPosition(1) = position(1)+1; 
            secondaryAxesPosition(2) = position(2); 
            secondaryAxesPosition(3) = position(3); 
            secondaryAxesPosition(4) = position(4); 
            setSecondaryAxesPositionInPixels( secondaryAxesPosition );       
        end    
        %Expand secondary axes on the X axis         
        if strcmp(currentModifier, 'alt') 
            position = getSecondaryAxesPositionInPixels(); 
            secondaryAxesPosition(3) = position(3)*(1 + 0.1); 
            if strcmpi( appDataStruct.global.zoomMode, 'off') 
                secondaryAxesPosition(4) = position(4); 
            else 
                %If 'freezeZoomAspectRatio' to 'on', be consistent 
                secondaryAxesPosition(4) = position(4)*(1 + 0.1);  
            end 
            secondaryAxesPosition(1) = position(1)-(-position(3)+secondaryAxesPosition(3))/2; 
            secondaryAxesPosition(2) = position(2)-(-position(4)+secondaryAxesPosition(4))/2; 
            setSecondaryAxesPositionInPixels( secondaryAxesPosition );   
        end 
         
    case {'uparrow'} % up arrow 
        %Move magnifier to the top 
        if isempty(currentModifier) 
            position = getMagnifierPositionInPixels(); 
            magnifierPosition(1) = position(1); 
            magnifierPosition(2) = position(2)+1; 
            magnifierPosition(3) = position(3); 
            magnifierPosition(4) = position(4); 
            setMagnifierPositionInPixels( magnifierPosition );          
        end 
        %Expand magnifier on the Y axis 
        if strcmp(currentModifier, 'shift') 
            position = getMagnifierPositionInPixels(); 
            if strcmpi( appDataStruct.global.zoomMode, 'off') 
                magnifierPosition(3) = position(3); 
            else 
                %If 'freezeZoomAspectRatio' to 'on', be consistent 
                magnifierPosition(3) = position(3)*(1 + 0.1);  
            end 
            magnifierPosition(4) = position(4)*(1 + 0.1); 
            magnifierPosition(1) = position(1)-(-position(3)+magnifierPosition(3))/2; 
            magnifierPosition(2) = position(2)-(-position(4)+magnifierPosition(4))/2; 
            setMagnifierPositionInPixels( magnifierPosition ); 
        end 
        %Move secondary axes to the top         
        if strcmp(currentModifier, 'control') 
            position = getSecondaryAxesPositionInPixels(); 
            secondaryAxesPosition(1) = position(1); 
            secondaryAxesPosition(2) = position(2)+1; 
            secondaryAxesPosition(3) = position(3); 
            secondaryAxesPosition(4) = position(4); 
            setSecondaryAxesPositionInPixels( secondaryAxesPosition );                   
        end  
        %Expand secondary axes on the Y axis         
        if strcmp(currentModifier, 'alt') 
            position = getSecondaryAxesPositionInPixels(); 
            if strcmpi( appDataStruct.global.zoomMode, 'off') 
                secondaryAxesPosition(3) = position(3); 
            else 
                %If 'freezeZoomAspectRatio' to 'on', be consistent 
                secondaryAxesPosition(3) = position(3)*(1 + 0.1);  
            end 
            secondaryAxesPosition(4) = position(4)*(1 + 0.1); 
            secondaryAxesPosition(1) = position(1)-(-position(3)+secondaryAxesPosition(3))/2; 
            secondaryAxesPosition(2) = position(2)-(-position(4)+secondaryAxesPosition(4))/2; 
            setSecondaryAxesPositionInPixels( secondaryAxesPosition );    
        end 
         
    case {'downarrow'} % down arrow 
        %Move magnifier to the bottom 
        if isempty(currentModifier) 
            position = getMagnifierPositionInPixels(); 
            magnifierPosition(1) = position(1); 
            magnifierPosition(2) = position(2)-1; 
            magnifierPosition(3) = position(3); 
            magnifierPosition(4) = position(4); 
            setMagnifierPositionInPixels( magnifierPosition );          
        end 
        %Compress magnifier on the Y axis 
        if strcmp(currentModifier, 'shift') 
            position = getMagnifierPositionInPixels(); 
            if strcmpi( appDataStruct.global.zoomMode, 'off') 
                magnifierPosition(3) = position(3); 
            else 
                %If 'freezeZoomAspectRatio' to 'on', be consistent 
                magnifierPosition(3) = position(3)*(1 - 0.1);  
            end 
            magnifierPosition(4) = position(4)*(1 - 0.1); 
            magnifierPosition(1) = position(1)-(-position(3)+magnifierPosition(3))/2; 
            magnifierPosition(2) = position(2)-(-position(4)+magnifierPosition(4))/2; 
            setMagnifierPositionInPixels( magnifierPosition ); 
        end 
        %Move secondary axes to the bottom         
        if strcmp(currentModifier, 'control') 
            position = getSecondaryAxesPositionInPixels(); 
            secondaryAxesPosition(1) = position(1); 
            secondaryAxesPosition(2) = position(2)-1; 
            secondaryAxesPosition(3) = position(3); 
            secondaryAxesPosition(4) = position(4); 
            setSecondaryAxesPositionInPixels( secondaryAxesPosition );          
        end  
        %Compress secondary axes on the Y axis         
        if strcmp(currentModifier, 'alt') 
            position = getSecondaryAxesPositionInPixels(); 
            if strcmpi( appDataStruct.global.zoomMode, 'off') 
                secondaryAxesPosition(3) = position(3); 
            else 
                %If 'freezeZoomAspectRatio' to 'on', be consistent 
                secondaryAxesPosition(3) = position(3)*(1 - 0.1);  
            end 
            secondaryAxesPosition(4) = position(4)*(1 - 0.1); 
            secondaryAxesPosition(1) = position(1)-(-position(3)+secondaryAxesPosition(3))/2; 
            secondaryAxesPosition(2) = position(2)-(-position(4)+secondaryAxesPosition(4))/2; 
            setSecondaryAxesPositionInPixels( secondaryAxesPosition );         
        end       
 
    case {'a'} % 'a' 
        %Debug info 
        if strcmp(currentModifier, 'control') 
            magnifierPosition = getMagnifierPositionInPixels(); 
            disp(sprintf('Magnifier position: [%g %g %g %g];', magnifierPosition(1), magnifierPosition(2), magnifierPosition(3), magnifierPosition(4)  )); 
            secondaryAxesPosition = getSecondaryAxesPositionInPixels();             
            disp(sprintf('Secondary axes position: [%g %g %g %g];', secondaryAxesPosition(1), secondaryAxesPosition(2), secondaryAxesPosition(3), secondaryAxesPosition(4)  ));             
        end 
 
    case {'q'} % 'a' 
        %additional xooming factors reseted 
        if strcmp(currentModifier, 'control') 
            appDataStruct.secondaryAxes.additionalZoomingFactor = [0 0]; 
        end 
                 
    case {'pageup'} % '+' 
        %Increase additional zooming factor on X-axis 
        if isempty(currentModifier) 
            appDataStruct.secondaryAxes.additionalZoomingFactor(1) = appDataStruct.secondaryAxes.additionalZoomingFactor(1) + 0.1; 
            if strcmpi( appDataStruct.global.zoomMode, 'on') 
                appDataStruct.secondaryAxes.additionalZoomingFactor(2) = appDataStruct.secondaryAxes.additionalZoomingFactor(2) + 0.1; 
            end 
        end 
        %Increase additional zooming factor on Y-axis         
        if strcmp(currentModifier, 'shift') 
            appDataStruct.secondaryAxes.additionalZoomingFactor(2) = appDataStruct.secondaryAxes.additionalZoomingFactor(2) + 0.1; 
            if strcmpi( appDataStruct.global.zoomMode, 'on') 
                appDataStruct.secondaryAxes.additionalZoomingFactor(1) = appDataStruct.secondaryAxes.additionalZoomingFactor(1) + 0.1; 
            end 
        end 
         
    case {'pagedown'} % '-' 
        %Redude additional zooming factor on X-axis 
        if isempty(currentModifier) 
            appDataStruct.secondaryAxes.additionalZoomingFactor(1) = appDataStruct.secondaryAxes.additionalZoomingFactor(1) - 0.1; 
            if strcmpi( appDataStruct.global.zoomMode, 'on') 
                appDataStruct.secondaryAxes.additionalZoomingFactor(2) = appDataStruct.secondaryAxes.additionalZoomingFactor(2) - 0.1; 
            end 
        end 
        %Redude additional zooming factor on Y-axis         
        if strcmp(currentModifier, 'shift') 
            appDataStruct.secondaryAxes.additionalZoomingFactor(2) = appDataStruct.secondaryAxes.additionalZoomingFactor(2) - 0.1; 
            if strcmpi( appDataStruct.global.zoomMode, 'on') 
                appDataStruct.secondaryAxes.additionalZoomingFactor(1) = appDataStruct.secondaryAxes.additionalZoomingFactor(1) - 0.1; 
            end 
        end         
         
    otherwise 
 
         
end 
 
%Update view limits on secondary axis 
refreshSecondaryAxisLimits(); 
 
%Update link between secondary axes and magnifier 
refreshMagnifierToSecondaryAxesLink(); 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NAME: ButtonMotionCallback 
%  
% PURPOSE: Service routine to ButtonMotion event. 
%  
% INPUT ARGUMENTS: 
%                    
% OUTPUT ARGUMENTS:  
%                    
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function ButtonMotionCallback(src,eventdata) 
 
global appDataStruct 
 
% pointerPos = get(appDataStruct.figure.handle, 'CurrentPoint'); 
% disp(sprintf('X: %g ; Y: %g', pointerPos(1), pointerPos(2)) ); 
 
getPointerArea(); 
 
%If Left mouse button not pressed, exit 
if  ~isfield(appDataStruct, 'ButtonDown') ||... 
    appDataStruct.ButtonDown == false 
   return; 
end 
 
%If Left mouse button pressed while the pointer is moving (drag) 
switch appDataStruct.pointerArea 
 
    case 'insideSecondaryAxis' 
        %Get current position of seconday axes (in pixels) 
        currentPositionOfSecondaryAxis = getSecondaryAxesPositionInPixels(); 
         
        %Get pointer position on figure's frame 
        currentPointerPositionOnFigureFrame = getPointerPositionOnFigureFrame(); 
         
        %Modify position 
        secondaryAxisPosition_W = currentPositionOfSecondaryAxis(3); 
        secondaryAxisPosition_H = currentPositionOfSecondaryAxis(4); 
        secondaryAxisPosition_X = currentPositionOfSecondaryAxis(1) + (-appDataStruct.pointerPositionOnButtonDown(1)+currentPointerPositionOnFigureFrame(1)); 
        secondaryAxisPosition_Y = currentPositionOfSecondaryAxis(2) + (-appDataStruct.pointerPositionOnButtonDown(2)+currentPointerPositionOnFigureFrame(2)); 
        appDataStruct.pointerPositionOnButtonDown = currentPointerPositionOnFigureFrame; 
         
        %Set initial position and size of secondary axes 
        setSecondaryAxesPositionInPixels( [... 
                                            secondaryAxisPosition_X,... 
                                            secondaryAxisPosition_Y,... 
                                            secondaryAxisPosition_W,... 
                                            secondaryAxisPosition_H... 
                                            ] ); 
         
    case 'insideMagnifier'                  
        %Get magnifier current position and size 
        currentPositionOfMagnifier = getMagnifierPositionInPixels(); 
         
        %Get pointer position on figure's frame 
        currentPointerPosition = getPointerPositionOnFigureFrame(); 
         
        %Modify magnifier position 
        magnifierPosition_W = currentPositionOfMagnifier(3); 
        magnifierPosition_H = currentPositionOfMagnifier(4); 
        magnifierPosition_X = currentPositionOfMagnifier(1) + (-appDataStruct.pointerPositionOnButtonDown(1)+currentPointerPosition(1)); 
        magnifierPosition_Y = currentPositionOfMagnifier(2) + (-appDataStruct.pointerPositionOnButtonDown(2)+currentPointerPosition(2)); 
        appDataStruct.pointerPositionOnButtonDown = currentPointerPosition; 
         
        %Set initial position and size of magnifying rectangle 
        setMagnifierPositionInPixels( [... 
                                          magnifierPosition_X... 
                                          magnifierPosition_Y... 
                                          magnifierPosition_W... 
                                          magnifierPosition_H... 
                                       ] ); 
         
        %Refresh zooming on secondary axis, based on magnifier position and extend 
        refreshSecondaryAxisLimits();     
            
    otherwise 
%         appDataStruct.pointerArea 
 
end 
 
%Update link between secondary axes and magnifier 
refreshMagnifierToSecondaryAxesLink(); 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NAME: ButtonDownCallback 
%  
% PURPOSE: Service routine to ButtonDown event. 
%  
% INPUT ARGUMENTS: 
%                    
% OUTPUT ARGUMENTS:  
%                    
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function ButtonDownCallback(src,eventdata) 
 
global appDataStruct 
 
if strcmpi( get(appDataStruct.figure.handle, 'SelectionType'), 'normal' ) 
     
    %Respond to left mous button 
    appDataStruct.ButtonDown = true; 
    %Get pointer position on figure's frame 
    appDataStruct.pointerPositionOnButtonDown = getPointerPositionOnFigureFrame(); 
 
elseif strcmpi( get(appDataStruct.figure.handle, 'SelectionType'), 'alt' ) 
     
    %Display contextual menu? 
     
end 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NAME: ButtonUpCallback 
%  
% PURPOSE: Service routine to ButtonUp event. 
%  
% INPUT ARGUMENTS: 
%                    
% OUTPUT ARGUMENTS:  
%                    
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function ButtonUpCallback(src,eventdata) 
 
global appDataStruct 
 
appDataStruct.ButtonDown = false; 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NAME: DeleteCallback 
%  
% PURPOSE: Service routine to Delete event. 
%  
% INPUT ARGUMENTS: 
%                    
% OUTPUT ARGUMENTS:  
%                    
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function DeleteCallback(src,eventdata) 
 
global appDataStruct; 
 
%Recover old callback handles 
set( appDataStruct.figure.handle, 'WindowButtonDownFcn', appDataStruct.figure.oldWindowButtonDownFcn); 
set( appDataStruct.figure.handle, 'WindowButtonUpFcn', appDataStruct.figure.oldWindowButtonUpFcn); 
set( appDataStruct.figure.handle, 'WindowButtonMotionFcn', appDataStruct.figure.oldWindowButtonMotionFcn); 
set( appDataStruct.figure.handle, 'KeyPressFcn', appDataStruct.figure.oldKeyPressFcn); 
set( appDataStruct.figure.handle, 'DeleteFcn', appDataStruct.figure.oldDeleteFcn); 
 
%Clear global variable when figure is closed 
clear global appDataStruct; 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NAME: getPointerPositionOnFigureFrame 
%  
% PURPOSE: determine if the position of the mouse pointer on the figure frame, in pixels. 
%  
% INPUT ARGUMENTS: 
%                   none 
% OUTPUT ARGUMENTS:  
%                   pointerPositionOnFigureFrame [double 1x2]: (X Y) 
%                   position 
%                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function pointerPositionOnFigureFrame = getPointerPositionOnFigureFrame() 
 
global appDataStruct 
 
%Get position of mouse pointer on screen 
defaultUnits = get(appDataStruct.figure.handle,'Units'); 
set(appDataStruct.figure.handle, 'Units', 'pixels'); 
pointerPositionOnFigureFrame = get(appDataStruct.figure.handle,'CurrentPoint'); 
set(appDataStruct.figure.handle,'Units', defaultUnits); 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NAME: getPointerArea 
%  
% PURPOSE: determine if the mouse pointer is on an active area. Change 
% pointer image if this is the case, and communicate the status. 
%  
% INPUT ARGUMENTS: 
%                   appDataStruct [struct 1x1]: global variable 
% OUTPUT ARGUMENTS:  
%                   change image of pointer (ACTION) 
%                   appDataStruct.pointerArea: ID of the active area 
%                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function getPointerArea() 
 
global appDataStruct 
 
%Get current pointer position on figure frame 
pointerPositionOnFigureFrame = getPointerPositionOnFigureFrame(); 
 
%Get current secondaryAxes position 
secondaryAxesPosition = getSecondaryAxesPositionInPixels(); 
 
%Get current magnifier position 
magnifierPosition = getMagnifierPositionInPixels(); 
 
%If mouse pointer on the secondary axis 
if  pointerPositionOnFigureFrame(1)>=secondaryAxesPosition(1) &&... 
    pointerPositionOnFigureFrame(1)<=secondaryAxesPosition(1)+secondaryAxesPosition(3) &&... 
    pointerPositionOnFigureFrame(2)>=secondaryAxesPosition(2) &&... 
    pointerPositionOnFigureFrame(2)<=secondaryAxesPosition(2)+secondaryAxesPosition(4) 
    %Pointer inside secondary axis 
    set(appDataStruct.figure.handle, 'Pointer', 'fleur'); 
     
    appDataStruct.pointerArea = 'insideSecondaryAxis'; 
     
elseif  pointerPositionOnFigureFrame(1)>=magnifierPosition(1) &&... 
        pointerPositionOnFigureFrame(1)<=magnifierPosition(1)+magnifierPosition(3) &&... 
        pointerPositionOnFigureFrame(2)>=magnifierPosition(2) &&... 
        pointerPositionOnFigureFrame(2)<=magnifierPosition(2)+magnifierPosition(4)  
    %Pointer inside magnifier 
    set(appDataStruct.figure.handle, 'Pointer', 'fleur'); 
     
    appDataStruct.pointerArea = 'insideMagnifier'; 
      
else 
    %Otherwise 
    set(appDataStruct.figure.handle, 'Pointer', 'arrow');   
     
    appDataStruct.pointerArea = 'none'; 
end 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NAME: getFigurePositionInPixels 
%  
% PURPOSE: obtain the position and size of the figure, relative to the 
% lower left corner of the screen, in pixels. 
%  
% INPUT ARGUMENTS: 
%                   none 
% OUTPUT ARGUMENTS:  
%                   position [double 1x4]:  
%                               X of lower left corner of the figure frame 
%                               Y of lower left corner of the figure frame 
%                               Width of the figure frame 
%                               Height of the figure frame 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function position = getFigurePositionInPixels() 
 
global appDataStruct 
 
defaultUnits = get(appDataStruct.figure.handle,'Units'); 
set(appDataStruct.figure.handle,'Units', 'pixels'); 
position = get(appDataStruct.figure.handle,'Position'); %pixels [ low bottom width height] 
set(appDataStruct.figure.handle,'Units', defaultUnits); 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NAME: getMainAxesPositionInPixels 
%  
% PURPOSE: obtain the position and size of the main axes, relative to the 
% lower left corner of the figure, in pixels. This fucntion locates the 
% lower-left corner of the displayed axes, accounting for all conditions of 
% DataAspectRatio and PlotBoxAspectRatio. 
%  
% INPUT ARGUMENTS: 
%                   none 
% OUTPUT ARGUMENTS:  
%                   position [double 1x4]:  
%                               X of lower left corner of the axis frame 
%                               Y of lower left corner of the axis frame 
%                               Width of the axis frame 
%                               Height of the axis frame 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function position = getMainAxesPositionInPixels() 
 
global appDataStruct 
 
%Characterize mainAxes in axes units 
mainAxisXLim = get( appDataStruct.mainAxes.handle, 'XLim' ); 
mainAxisYLim = get( appDataStruct.mainAxes.handle, 'YLim' ); 
spanX = determineSpan(mainAxisXLim(1), mainAxisXLim(2) ); 
spanY = determineSpan(mainAxisYLim(1), mainAxisYLim(2)); 
 
%Capture default units of mainAxes, and fix units to 'pixels' 
defaultUnits = get(appDataStruct.mainAxes.handle,'Units'); 
set(appDataStruct.mainAxes.handle,'Units', 'pixels'); 
 
%Obtain values in 'pixels' 
mainAxesPosition = get(appDataStruct.mainAxes.handle, 'Position'); 
dataAspectRatioMode  = get(appDataStruct.mainAxes.handle, 'DataAspectRatioMode'); 
dataAspectRatio = get(appDataStruct.mainAxes.handle, 'DataAspectRatio'); 
plotBoxAspectRatioMode = get(appDataStruct.mainAxes.handle, 'PlotBoxAspectRatioMode'); 
plotBoxAspectRatio = get(appDataStruct.mainAxes.handle, 'PlotBoxAspectRatio'); 
 
%Determine correction values 
dataAspectRatioLimits = (spanX/dataAspectRatio(1))/(spanY/dataAspectRatio(2)); 
plotBoxAspectRatioRelation = plotBoxAspectRatio(1)/plotBoxAspectRatio(2); 
mainAxesRatio = mainAxesPosition(3)/mainAxesPosition(4); 
 
%Id DataAspectRatio to auto and PlotBoxAspectRatio to auto 
if ~strcmpi( dataAspectRatioMode, 'manual') && ~strcmpi( plotBoxAspectRatioMode, 'manual') 
        
    %Recover default units of mainAxes 
    set(appDataStruct.mainAxes.handle,'Units', defaultUnits); 
 
    %Obtain 'real' position from a temporal axes 
    temporalAxes = axes('Visible', 'off'); 
    set(temporalAxes, 'Units', 'pixels'); 
    set(temporalAxes, 'Position', mainAxesPosition); 
    position = get(temporalAxes, 'Position'); 
    delete(temporalAxes); 
     
    return; 
end  
 
%If DataAspectRatio to manual 
if strcmpi( dataAspectRatioMode, 'manual') 
    if dataAspectRatioLimits <= mainAxesRatio        
        position(4) = mainAxesPosition(4); 
        position(3) = mainAxesPosition(4) * dataAspectRatioLimits; 
        position(2) = mainAxesPosition(2); 
        position(1) = mainAxesPosition(1) + (mainAxesPosition(3) - position(3))/2; 
         
    else  
        position(1) = mainAxesPosition(1); 
        position(3) = mainAxesPosition(3); 
        position(4) = mainAxesPosition(3)/dataAspectRatioLimits; 
        position(2) = mainAxesPosition(2) + (mainAxesPosition(4) - position(4))/2;      
         
    end 
elseif strcmpi( plotBoxAspectRatioMode, 'manual') 
    % Or PlotBoxAspectRatio to manual 
    if plotBoxAspectRatioRelation <= mainAxesRatio        
        position(4) = mainAxesPosition(4); 
        position(3) = mainAxesPosition(4) * plotBoxAspectRatioRelation; 
        position(2) = mainAxesPosition(2); 
        position(1) = mainAxesPosition(1) + (mainAxesPosition(3) - position(3))/2; 
         
    else 
        position(1) = mainAxesPosition(1); 
        position(3) = mainAxesPosition(3); 
        position(4) = mainAxesPosition(3)/plotBoxAspectRatioRelation; 
        position(2) = mainAxesPosition(2) + (mainAxesPosition(4) - position(4))/2;   
         
    end 
end 
 
%Recover default units of mainAxes 
set(appDataStruct.mainAxes.handle,'Units', defaultUnits); 
 
%Obtain 'real' position from a temporal axes 
temporalAxes = axes('Visible', 'off'); 
set(temporalAxes, 'Units', 'pixels'); 
set(temporalAxes, 'Position', position ); 
position = get(temporalAxes, 'Position'); 
delete(temporalAxes); 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NAME: getMagnifierPositionInPixels 
%  
% PURPOSE: obtain the position (of the lower left corner) and size of the  
% magnifier, relative to the lower left corner of the figure, in pixels. 
%  
% INPUT ARGUMENTS: 
%                   none 
% OUTPUT ARGUMENTS:  
%                   position [double 1x4]:  
%                               X of lower left corner of the magnifier 
%                               Y of lower left corner of the magnifier 
%                               Width of the magnifier 
%                               Height of the magnifier 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function position = getMagnifierPositionInPixels() 
 
global appDataStruct; 
 
defaultUnits = get(appDataStruct.magnifier.handle, 'Units'); 
set(appDataStruct.magnifier.handle, 'Units', 'pixels'); 
position = get(appDataStruct.magnifier.handle, 'Position'); 
set(appDataStruct.magnifier.handle, 'Units', defaultUnits ); 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NAME: getSecondaryAxesPositionInPixels 
%  
% PURPOSE: obtain the position and size of the secondary axis, relative to the 
% lower left corner of the figure, in pixels. Includes legends and axes 
% numbering 
%  
% INPUT ARGUMENTS: 
%                   none 
% OUTPUT ARGUMENTS:  
%                   position [double 1x4]:  
%                               X of lower left corner of the axis frame 
%                               Y of lower left corner of the axis frame 
%                               Width of the axis frame 
%                               Height of the axis frame 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function position = getSecondaryAxesPositionInPixels() 
 
global appDataStruct 
 
defaultUnits = get(appDataStruct.secondaryAxes.handle,'Units');  
set(appDataStruct.secondaryAxes.handle,'Units', 'pixels');  
position = get(appDataStruct.secondaryAxes.handle,'Position'); %[ left bottom width height] 
set(appDataStruct.secondaryAxes.handle,'Units', defaultUnits);  
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NAME: setSecondaryAxesPositionInPixels 
%  
% PURPOSE: fix the position and size of the secondary axis, relative to the 
% lower left corner of the figure, in pixels.  
%  
% INPUT ARGUMENTS: 
%                   position [double 1x4]:  
%                               X of lower left corner of the axis frame 
%                               Y of lower left corner of the axis frame 
%                               Width of the axis frame 
%                               Height of the axis frame 
% OUTPUT ARGUMENTS:  
%                   none      
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function setSecondaryAxesPositionInPixels( position ) 
 
global appDataStruct 
 
%Get position of secondary axes 
defaultUnits = get(appDataStruct.secondaryAxes.handle,'Units'); 
set(appDataStruct.secondaryAxes.handle,'Units', 'pixels'); 
set(    appDataStruct.secondaryAxes.handle,... 
        'Position', [... 
                    position(1),... 
                    position(2),... 
                    position(3),... 
                    position(4)... 
                    ]... 
    );  
% tightInset = get( appDataStruct.secondaryAxes.handle, 'TightInset' );  
set(appDataStruct.secondaryAxes.handle,'Units', defaultUnits); 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NAME: setMagnifierPositionInPixels 
%  
% PURPOSE: fix the position and size of the magnifier, relative to the 
% lower left corner of the figure, in pixels.  
%  
% INPUT ARGUMENTS: 
%                   position [double 1x4]:  
%                               X of lower left corner of the magnifier 
%                               Y of lower left corner of the magnifier 
%                               Width of the magnifier frame 
%                               Height of the magnifier frame 
% OUTPUT ARGUMENTS:  
%                   none      
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function setMagnifierPositionInPixels( position ) 
 
 
global appDataStruct 
 
%Limit position of magnifier within the main axes 
mainAxesPosition = getMainAxesPositionInPixels(); 
if position(1)<mainAxesPosition(1) 
    position(1) = mainAxesPosition(1); 
end 
if position(1)+position(3)>mainAxesPosition(1)+mainAxesPosition(3) 
    position(1) = mainAxesPosition(1)+mainAxesPosition(3)-position(3); 
end 
if position(2)<mainAxesPosition(2) 
    position(2) = mainAxesPosition(2); 
end 
if position(2)+position(4)>mainAxesPosition(2)+mainAxesPosition(4) 
    position(2) = mainAxesPosition(2)+mainAxesPosition(4)-position(4); 
end 
 
%Create of set magnifier 
if ~isfield(appDataStruct, 'magnifier') ||... 
   (isfield(appDataStruct, 'magnifier') && ~isfield(appDataStruct.magnifier, 'handle')) ||... 
   (isfield(appDataStruct, 'magnifier') && isfield(appDataStruct.magnifier, 'handle') && ~ishandle(appDataStruct.magnifier.handle)) ||... 
   (isfield(appDataStruct, 'magnifier') && isfield(appDataStruct.magnifier, 'handle') && isempty(appDataStruct.magnifier.handle))  
    
    if strcmpi(appDataStruct.magnifier.shape, 'rectangle') 
        appDataStruct.magnifier.handle =... 
                       annotation(  'rectangle',... 
                                    'Units', 'pixels',... 
                                    'Position', position,... 
                                    'LineWidth', appDataStruct.global.edgeWidth,... 
                                    'LineStyle','-',... 
                                    'EdgeColor',appDataStruct.global.edgeColor...                                 
                                    ); 
    end 
    if strcmpi(appDataStruct.magnifier.shape, 'ellipse') 
        appDataStruct.magnifier.handle =... 
                       annotation(  'ellipse',... 
                                    'Units', 'pixels',... 
                                    'Position', position,... 
                                    'LineWidth', appDataStruct.global.edgeWidth,... 
                                    'LineStyle','-',... 
                                    'EdgeColor',appDataStruct.global.edgeColor...                                 
                                    ); 
    end 
     
else 
   set( appDataStruct.magnifier.handle, 'Position', position );  
end                       
         
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NAME: refreshMagnifierToSecondaryAxesLink 
%  
% PURPOSE: Updates the line connection between the magnifier and the secondary axes. 
%  
% INPUT ARGUMENTS: 
% 
% OUTPUT ARGUMENTS:  
%                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function refreshMagnifierToSecondaryAxesLink() 
         
global appDataStruct; 
 
%Don't display link if not requestred 
if strcmpi( appDataStruct.link(1).displayLinkStyle, 'none') 
    return; 
end 
 
%Get position and size of figure in pixels 
figurePosition = getFigurePositionInPixels(); 
 
%Get position and size of secondary axes in pixels 
secondaryAxesPosition = getSecondaryAxesPositionInPixels(); 
 
defaultUnits = get(appDataStruct.secondaryAxes.handle, 'Units'); 
set(appDataStruct.secondaryAxes.handle, 'Units', 'pixels'); 
tightInset = get(appDataStruct.secondaryAxes.handle, 'TightInset'); 
set(appDataStruct.secondaryAxes.handle, 'Units', defaultUnits); 
 
%Get position and size of secondary axes in pixels 
magnifierPosition = getMagnifierPositionInPixels(); 
    
 
if strcmpi( appDataStruct.link(1).displayLinkStyle, 'straight') 
 
    %Magnifier Hot points 
    magnifierHotPoints = [... 
        magnifierPosition(1) + magnifierPosition(3)/2   magnifierPosition(2);... 
        magnifierPosition(1) + magnifierPosition(3)     magnifierPosition(2)+magnifierPosition(4)/2;... 
        magnifierPosition(1) + magnifierPosition(3)/2   magnifierPosition(2)+magnifierPosition(4);... 
        magnifierPosition(1)                            magnifierPosition(2)+magnifierPosition(4)/2;... 
        ]; 
 
    %Secondary axes Hot points 
    secondaryAxesHotPoints = [... 
        secondaryAxesPosition(1) + secondaryAxesPosition(3)/2   secondaryAxesPosition(2) - tightInset(2) - 2;... 
        secondaryAxesPosition(1) + secondaryAxesPosition(3)     secondaryAxesPosition(2)+secondaryAxesPosition(4)/2;... 
        secondaryAxesPosition(1) + secondaryAxesPosition(3)/2   secondaryAxesPosition(2)+secondaryAxesPosition(4);...     
        secondaryAxesPosition(1) - tightInset(1) - 2            secondaryAxesPosition(2)+secondaryAxesPosition(4)/2;...     
        ]; 
     
    %Minimize distance between hot spots 
    L1 = size(magnifierHotPoints, 1); 
    L2 = size(secondaryAxesHotPoints, 1); 
    [iMagnifierHotPoints iSecondaryAxesHotPoints] = meshgrid(1:L1, 1:L2); 
    D2  =   ( magnifierHotPoints(iMagnifierHotPoints(:),1) - secondaryAxesHotPoints(iSecondaryAxesHotPoints(:),1) ).^2 + ... 
            ( magnifierHotPoints(iMagnifierHotPoints(:),2) - secondaryAxesHotPoints(iSecondaryAxesHotPoints(:),2) ).^2; 
 
    [C,I] = sort( D2, 'ascend' ); 
 
    X(1) = magnifierHotPoints(iMagnifierHotPoints(I(1)),1); 
    Y(1) = magnifierHotPoints(iMagnifierHotPoints(I(1)),2); 
    X(2) = secondaryAxesHotPoints(iSecondaryAxesHotPoints(I(1)),1); 
    Y(2) = secondaryAxesHotPoints(iSecondaryAxesHotPoints(I(1)),2); 
     
    %Plot/update line 
    if ~isfield(appDataStruct, 'link') ||... 
       (isfield(appDataStruct, 'link') && ~isfield(appDataStruct.link, 'handle')) ||... 
       (isfield(appDataStruct, 'link') && isfield(appDataStruct.link, 'handle') && ~ishandle(appDataStruct.link.handle)) ||... 
       (isfield(appDataStruct, 'link') && isfield(appDataStruct.link, 'handle') && isempty(appDataStruct.link.handle)) 
 
        appDataStruct.link.handle = annotation( 'line', X/figurePosition(3), Y/figurePosition(4),... 
                                                'LineWidth', appDataStruct.global.edgeWidth,... 
                                                'Color', appDataStruct.global.edgeColor ); 
 
    else 
 
        set(appDataStruct.link.handle, 'X', X/figurePosition(3), 'Y', Y/figurePosition(4));  
 
    end 
end 
 
if strcmpi( appDataStruct.link(1).displayLinkStyle, 'edges') 
    %Magnifier Hot points 
    magnifierHotPoints = [... 
        magnifierPosition(1) - 3                            magnifierPosition(2);... 
        magnifierPosition(1) + magnifierPosition(3)     magnifierPosition(2);...     
        magnifierPosition(1) + magnifierPosition(3)     magnifierPosition(2)+magnifierPosition(4);...     
        magnifierPosition(1) - 3                           magnifierPosition(2)+magnifierPosition(4)... 
        ]; 
 
    %Secondary axes Hot points 
    secondaryAxesHotPoints = [... 
        secondaryAxesPosition(1)                                secondaryAxesPosition(2);... 
        secondaryAxesPosition(1) + secondaryAxesPosition(3)     secondaryAxesPosition(2);... 
        secondaryAxesPosition(1) + secondaryAxesPosition(3)     secondaryAxesPosition(2)+secondaryAxesPosition(4);... 
        secondaryAxesPosition(1)                                secondaryAxesPosition(2)+secondaryAxesPosition(4)... 
        ]; 
     
     
    for i=1:4 
        X(1) = magnifierHotPoints(i,1); 
        Y(1) = magnifierHotPoints(i,2); 
        X(2) = secondaryAxesHotPoints(i,1); 
        Y(2) = secondaryAxesHotPoints(i,2); 
 
        %If intersection with secondary Axes bottom edge 
%         intersectionPoint = intersectionPointInPixels(... 
%                     [X(1) Y(1) X(2) Y(2)], ... 
%                     [   secondaryAxesPosition(1)... 
%                         secondaryAxesPosition(2)... 
%                         secondaryAxesPosition(1)+secondaryAxesPosition(3)... 
%                         secondaryAxesPosition(2) ]... 
%                         ); 
%         if ~isempty(intersectionPoint)                       
%             D2_1 = (X(1)-X(2))^2 + (Y(1)-Y(2))^2; 
%             D2_2 = (X(1)-intersectionPoint(1))^2 + (Y(1)-intersectionPoint(2))^2; 
%             if D2_2<D2_1 
%                 %link to intersecting point 
%                 X(2) = intersectionPoint(1); 
%                 Y(2) = intersectionPoint(2); 
%             end 
%         end   
%  
%         %If intersection with secondary Axes top edge 
%         intersectionPoint = intersectionPointInPixels(... 
%                     [X(1) Y(1) X(2) Y(2)], ... 
%                     [   secondaryAxesPosition(1)... 
%                         secondaryAxesPosition(2)+secondaryAxesPosition(4)... 
%                         secondaryAxesPosition(1)+secondaryAxesPosition(3)... 
%                         secondaryAxesPosition(2)+secondaryAxesPosition(4) ]... 
%                         ); 
%         if ~isempty(intersectionPoint) 
%             D2_1 = (X(1)-X(2))^2 + (Y(1)-Y(2))^2; 
%             D2_2 = (X(1)-intersectionPoint(1))^2 + (Y(1)-intersectionPoint(2))^2; 
%             if D2_2<D2_1 
%                 %link to intersecting point 
%                 X(2) = intersectionPoint(1); 
%                 Y(2) = intersectionPoint(2); 
%             end 
%         end  
%  
%         %If intersection with secondary Axes left edge 
%         intersectionPoint = intersectionPointInPixels(... 
%                     [X(1) Y(1) X(2) Y(2)], ... 
%                     [   secondaryAxesPosition(1)... 
%                         secondaryAxesPosition(2)... 
%                         secondaryAxesPosition(1)... 
%                         secondaryAxesPosition(2)+secondaryAxesPosition(4) ]... 
%                         ); 
%         if ~isempty(intersectionPoint) 
%             D2_1 = (X(1)-X(2))^2 + (Y(1)-Y(2))^2; 
%             D2_2 = (X(1)-intersectionPoint(1))^2 + (Y(1)-intersectionPoint(2))^2; 
%             if D2_2<D2_1 
%                 %link to intersecting point 
%                 X(2) = intersectionPoint(1); 
%                 Y(2) = intersectionPoint(2); 
%             end 
%         end  
%  
%         %If intersection with secondary Axes right edge 
%         intersectionPoint = intersectionPointInPixels(... 
%                     [X(1) Y(1) X(2) Y(2)], ... 
%                     [   secondaryAxesPosition(1)+secondaryAxesPosition(3)... 
%                         secondaryAxesPosition(2)... 
%                         secondaryAxesPosition(1)+secondaryAxesPosition(3)... 
%                         secondaryAxesPosition(2)+secondaryAxesPosition(4) ]... 
%                         ); 
%         if ~isempty(intersectionPoint) 
%             D2_1 = (X(1)-X(2))^2 + (Y(1)-Y(2))^2; 
%             D2_2 = (X(1)-intersectionPoint(1))^2 + (Y(1)-intersectionPoint(2))^2; 
%             if D2_2<D2_1 
%                 %link to intersecting point 
%                 X(2) = intersectionPoint(1); 
%                 Y(2) = intersectionPoint(2); 
%             end 
%         end  
 
        %Plot/update line 
        if length(appDataStruct.link)<i || ~isfield(appDataStruct.link(i), 'handle') ||... 
           isempty(appDataStruct.link(i).handle) || ~ishandle(appDataStruct.link(i).handle) 
         
            appDataStruct.link(i).handle = annotation( 'line', X/figurePosition(3), Y/figurePosition(4),... 
                                                   'LineWidth', appDataStruct.global.edgeWidth,... 
                                                   'LineStyle', ':',... 
                                                    'Color', appDataStruct.global.edgeColor ); 
         
        else 
            set(appDataStruct.link(i).handle, 'X', X/figurePosition(3), 'Y', Y/figurePosition(4));  
        end 
             
    end 
end 
 
if strcmpi( appDataStruct.link(1).displayLinkStyle, 'elbow') 
 
    %Magnifier Hot points 
    magnifierHotPoints = [... 
        magnifierPosition(1) + magnifierPosition(3)/2   magnifierPosition(2);... 
        magnifierPosition(1) + magnifierPosition(3)     magnifierPosition(2)+magnifierPosition(4)/2;... 
        magnifierPosition(1) + magnifierPosition(3)/2   magnifierPosition(2)+magnifierPosition(4);... 
        magnifierPosition(1)                            magnifierPosition(2)+magnifierPosition(4)/2;... 
        ]; 
 
    %Secondary axes Hot points 
    secondaryAxesHotPoints = [... 
        secondaryAxesPosition(1) + secondaryAxesPosition(3)/2   secondaryAxesPosition(2) - tightInset(2) - 2;... 
        secondaryAxesPosition(1) + secondaryAxesPosition(3)     secondaryAxesPosition(2)+secondaryAxesPosition(4)/2;... 
        secondaryAxesPosition(1) + secondaryAxesPosition(3)/2   secondaryAxesPosition(2)+secondaryAxesPosition(4);...     
        secondaryAxesPosition(1) - tightInset(1) - 2            secondaryAxesPosition(2)+secondaryAxesPosition(4)/2;...     
        ]; 
     
     
    %Allowed connections 
%     iMagnifierHotPoints(1) = 1;  
%     iSecondaryAxesHotPoints(1) = 4; 
%     iMagnifierHotPoints(2) = 1;  
%     iSecondaryAxesHotPoints(2) = 2; 
%     iMagnifierHotPoints(3) = 2;  
%     iSecondaryAxesHotPoints(3) = 3; 
%     iMagnifierHotPoints(4) = 2;  
%     iSecondaryAxesHotPoints(4) = 1; 
%     iMagnifierHotPoints(5) = 3;  
%     iSecondaryAxesHotPoints(5) = 4; 
%     iMagnifierHotPoints(6) = 3;  
%     iSecondaryAxesHotPoints(6) = 2; 
%     iMagnifierHotPoints(7) = 4;  
%     iSecondaryAxesHotPoints(7) = 1; 
%     iMagnifierHotPoints(8) = 4;  
%     iSecondaryAxesHotPoints(8) = 3; 
%     iMagnifierHotPoints(9) = 1;  
%     iSecondaryAxesHotPoints(9) = 3; 
%     iMagnifierHotPoints(10) = 2;  
%     iSecondaryAxesHotPoints(10) = 4; 
%     iMagnifierHotPoints(11) = 3;  
%     iSecondaryAxesHotPoints(11) = 1; 
%     iMagnifierHotPoints(12) = 4;  
%     iSecondaryAxesHotPoints(12) = 2; 
     
    %Minimize distance between hot spots 
    L1 = size(magnifierHotPoints, 1); 
    L2 = size(secondaryAxesHotPoints, 1); 
    [iMagnifierHotPoints iSecondaryAxesHotPoints] = meshgrid(1:L1, 1:L2); 
    D2  =   ( magnifierHotPoints(iMagnifierHotPoints(:),1) - secondaryAxesHotPoints(iSecondaryAxesHotPoints(:),1) ).^2 + ... 
            ( magnifierHotPoints(iMagnifierHotPoints(:),2) - secondaryAxesHotPoints(iSecondaryAxesHotPoints(:),2) ).^2; 
 
    [C,I] = sort( D2, 'ascend' ); 
 
    X(1) = magnifierHotPoints(iMagnifierHotPoints(I(1)),1); 
    Y(1) = magnifierHotPoints(iMagnifierHotPoints(I(1)),2); 
    X(2) = secondaryAxesHotPoints(iSecondaryAxesHotPoints(I(1)),1); 
    Y(2) = secondaryAxesHotPoints(iSecondaryAxesHotPoints(I(1)),2); 
     
    %Plot/update line 
    if ~isfield(appDataStruct, 'link') ||... 
       (isfield(appDataStruct, 'link') && ~isfield(appDataStruct.link, 'handle')) ||... 
       (isfield(appDataStruct, 'link') && isfield(appDataStruct.link, 'handle') && ~ishandle(appDataStruct.link.handle)) ||... 
       (isfield(appDataStruct, 'link') && isfield(appDataStruct.link, 'handle') && isempty(appDataStruct.link.handle)) 
 
        appDataStruct.link.handle = annotation( 'line', X/figurePosition(3), Y/figurePosition(4),... 
                                                'LineWidth', appDataStruct.global.edgeWidth,... 
                                                'Color', appDataStruct.global.edgeColor ); 
 
    else 
 
        set(appDataStruct.link.handle, 'X', X/figurePosition(3), 'Y', Y/figurePosition(4));  
 
    end 
end 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NAME: intersectionPointInPixels 
%  
% PURPOSE: Computes constrained intersection of two lines in pixels, on the 2D space 
%  
% INPUT ARGUMENTS: 
%                   line1 [double 1x4]: [Xstart Ystart Xend Yend] 
%                   line2 [double 1x4]: [Xstart Ystart Xend Yend] 
% 
% OUTPUT ARGUMENTS:  
%                   intersectionPont [double 1x2]: [X Y] intersection 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function intersectionPoint = intersectionPointInPixels( line1, line2) 
                     
%Cartessian caracterization of line 1                     
X1(1) = line1(1); 
Y1(1) = line1(2);  
X1(2) = line1(3); 
Y1(2) = line1(4);  
 
a1 = (Y1(2) - Y1(1)) / (X1(2) - X1(1)); 
b1 = Y1(1) - X1(1)*a1; 
 
%Cartessian caracterization of line 2                     
X2(1) = line2(1); 
Y2(1) = line2(2);  
X2(2) = line2(3); 
Y2(2) = line2(4);  
 
a2 = (Y2(2) - Y2(1)) / (X2(2) - X2(1)); 
b2 = Y2(1) - X2(1)*a2; 
 
%Intersection 
if isfinite(a1) && isfinite(a2) 
    intersectionPoint(1) = (b2-b1) / (a1-a2); 
    intersectionPoint(2) = intersectionPoint(1)*a1 + b1; 
end 
%Pathologic case 1 (line2 x=constant) 
if isfinite(a1) && ~isfinite(a2) 
    intersectionPoint(1) = X2(1); 
    intersectionPoint(2) = intersectionPoint(1)*a1 + b1; 
end 
%Pathologic case 2 (line1 x=constant) 
if ~isfinite(a1) && isfinite(a2) 
    intersectionPoint(1) = X1(1); 
    intersectionPoint(2) = intersectionPoint(1)*a2 + b2; 
end 
 
if intersectionPoint(1)<min([X1(1) X2(1)]) ||... 
   intersectionPoint(1)>max([X1(2) X2(2)]) ||... 
   intersectionPoint(2)<min([Y1(1) Y2(1)]) ||... 
   intersectionPoint(2)>max([Y1(2) Y2(2)])  
         
    intersectionPoint = []; 
end 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NAME: computeSecondaryAxesDefaultPosition 
%  
% PURPOSE: obtain the default position and size of the secondary axis, relative to the 
% lower left corner of the figure, in pixels. Includes legends and axes 
% numbering 
%  
% INPUT ARGUMENTS: 
%                   none 
% OUTPUT ARGUMENTS:  
%                   position [double 1x4]:  
%                               X of lower left corner of the axis frame 
%                               Y of lower left corner of the axis frame 
%                               Width of the axis frame 
%                               Height of the axis frame 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function defaultPosition = computeSecondaryAxesDefaultPosition() 
 
global appDataStruct; 
 
% If image, defualt aspect ratio of magnifier and secondary axes to [1 1] 
childHandle = get(appDataStruct.mainAxes.handle, 'Children'); 
BOX = appDataStruct.secondaryAxes.initialPosition;
X_LEN = get(childHandle,'XData');
Y_LEN = get(childHandle,'YData');
my_ration = BOX./[X_LEN(2),X_LEN(2),X_LEN(2),X_LEN(2)];
%Get position and size of main Axis (left & bottom relative to figure frame) 
mainAxesPosition = getMainAxesPositionInPixels();
% my_rationmy_ration*my_ration
%Set initial position and size for secondary axis
secondaryAxisPosition_W = mainAxesPosition(3)*my_ration(3);
secondaryAxisPosition_H = mainAxesPosition(4)*my_ration(4);
% secondaryAxisPosition_X = mainAxesPosition(1)+mainAxesPosition(3)-secondaryAxisPosition_W-10;
% secondaryAxisPosition_Y = mainAxesPosition(2)+mainAxesPosition(4)-secondaryAxisPosition_H-10;
secondaryAxisPosition_X = mainAxesPosition(1)*(1+my_ration(1));
secondaryAxisPosition_Y = mainAxesPosition(2)*(1+my_ration(2));


defaultPosition = [... 
                        secondaryAxisPosition_X... 
                        secondaryAxisPosition_Y... 
                        secondaryAxisPosition_W... 
                        secondaryAxisPosition_H... 
                        ]; 
                         
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NAME: computeMagnifierDefaultPosition 
%  
% PURPOSE: obtain the default position and size of the magnifier, relative to the 
% lower left corner of the figure, in pixels. Includes legends and axes 
% numbering 
%  
% INPUT ARGUMENTS: 
%                   none 
% OUTPUT ARGUMENTS:  
%                   position [double 1x4]:  
%                               X of lower left corner of the rectangle 
%                               Y of lower left corner of the rectangle 
%                               Width of the rectangle 
%                               Height of the rectangle 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function defaultPosition = computeMagnifierDefaultPosition() 
 
global appDataStruct; 


% If image, defualt aspect ratio of magnifier and secondary axes to [1 1] 
childHandle = get(appDataStruct.mainAxes.handle, 'Children'); 
 
%Set initial position and size of magnifying rectangle 
mainAxisXLim = get( appDataStruct.mainAxes.handle, 'XLim' ); 
mainAxisYLim = get( appDataStruct.mainAxes.handle, 'YLim' ); 
mainAxesPositionInPixels = getMainAxesPositionInPixels(); 
xMainAxisUnits2PixelsFactor = mainAxesPositionInPixels(3)/determineSpan( mainAxisXLim(1), mainAxisXLim(2) );
yMainAxisUnits2PixelsFactor = mainAxesPositionInPixels(4)/determineSpan( mainAxisYLim(1), mainAxisYLim(2) );


%Get main axis position and dimensions, in 
BOX = appDataStruct.magnifier.initialPosition;

magnifierPosition_W = BOX(3)*xMainAxisUnits2PixelsFactor;
magnifierPosition_H = BOX(4)*yMainAxisUnits2PixelsFactor;
magnifierPosition_X = BOX(1)*xMainAxisUnits2PixelsFactor - magnifierPosition_W;
magnifierPosition_Y = BOX(2)*yMainAxisUnits2PixelsFactor - magnifierPosition_H;

defaultPosition = [... 
                        magnifierPosition_X+mainAxesPositionInPixels(1)... 
                        magnifierPosition_Y+mainAxesPositionInPixels(2)... 
                        magnifierPosition_W... 
                        magnifierPosition_H... 
                        ]; 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% NAME: setDefaultProperties 
%  
% PURPOSE: Set value for default properties 
%  
% INPUT ARGUMENTS: 
%                   none 
% OUTPUT ARGUMENTS:  
%                   none 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                     
function setDefaultProperties() 
 
global appDataStruct 
 
appDataStruct.figure = []; 
appDataStruct.mainAxes = []; 
appDataStruct.secondaryAxes = [];  
appDataStruct.magnifier = []; 
appDataStruct.link = []; 
appDataStruct.global = []; 
 
%Get handle for figure  
appDataStruct.figure.handle = gcf; 
 
%Get handles for main axis  
appDataStruct.mainAxes.handle = []; 
 
%Get handles for secondary axis 
appDataStruct.secondaryAxes.handle = []; 
 
%Default units 
appDataStruct.global.units = 'pixels'; 
 
%Default operation mode 
appDataStruct.global.mode = 'interactive'; 
 
%Default line display mode 
appDataStruct.link.displayLinkStyle = 'straight'; 
 
%Default line width for lines of displayed objects 
appDataStruct.global.edgeWidth = 1; 
 
%Default color for lines of displayed objects 
appDataStruct.global.edgeColor = 'black'; 
 
%Default zoom mode 
appDataStruct.global.zoomMode = 'off'; 
 
%Default color of secondary axes face 
appDataStruct.secondaryAxes.faceColor = 'white'; 
 
%Default secondary axes position 
appDataStruct.secondaryAxes.initialPosition = []; 
 
%Default secondary axes XLim 
appDataStruct.secondaryAxes.initialXLim = []; 
 
%Default secondary axes YLim 
appDataStruct.secondaryAxes.initialYLim = []; 
 
%Default magnifier position 
appDataStruct.magnifier.initialPosition = []; 
 
%Default magnifier shape 
appDataStruct.magnifier.shape = 'rectangle'; 
 
%Default additional zooming factors 
appDataStruct.secondaryAxes.additionalZoomingFactor = [0 0]; 
 

 