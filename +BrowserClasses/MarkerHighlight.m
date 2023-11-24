%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MERIDIAN Detection Viewer - Class "MarkerHighlight"
%   Written by Wilfried Beslin
%   Last update Oct. 21, 2020, using MATLAB R2018b
%
%   Description:
%   Class defining and controlling graphics involved with highlighting
%   a marker while the viewer is in "edit" mode. 
%   "Highlighting" is accomplished by drawing four Lines over a marker: one
%   line for each edge (this approach makes it easy to highlight the edges
%   individually for resizing).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MarkerHighlight < handle
    %% PROPERTIES =========================================================
    properties (SetAccess = private)
        hAxes       % handle of parent Axes object
        hMarker     % handle of marker being highlighted
        signalType  % may be 'Detection', 'UndetectedUpcall', or 'Annotation'
        signalIdx   % index of the highlighted marker relative to the complete list of signals of that type
        screenIdx   % index of the highlighted marker relative to the list of on-screen markers of that type
        lineData    % 4-by-4 table containing line objects and associated info, including edge type (DT or F) and mouse hit information
    end
    properties (Dependent)
        active          % indicates if a marker is highlighted or not
        hasSelection    % indicates if a marker has been selected for highlighting or not
        isValid         % indicates if object has a valid marker handle or not
        isDisplayed     % indicates if highlight graphics exist or not
        editFreq        % determines if FMin and FMax can be adjusted or not
        DTMin           % marker start date/time    
        DTMax           % marker end date/time      
        FMin            % marker min frequency      
        FMax            % marker max frequency      
    end
    properties (Access = private, Constant)
        bufferProp = 0.005; % +/- proportion of the axes dimensions to use as buffer region for edge detection
        edgeParams = table(...
            {':'; '-'},...
            [2.0; 2.5],...
            [1.5; 2.0],...
            'VariableNames',{'LineStyle','LineWidthMultiplier','ColorMultiplier'},...
            'RowNames',{'Normal','Hover'});
    end
    
    %% METHODS - PUBLIC ===================================================
    methods
        % Constructor -----------------------------------------------------
        function obj = MarkerHighlight(ha)
        % Creates a "MarkerHighlight" object.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            obj.hAxes = ha;
            obj.lineData = table(...
                repelem(matlab.graphics.GraphicsPlaceholder,4,1),...
                [true;true;false;false],...
                [false;false;true;true],...
                false(4,1),...
                'RowNames',{'DTMin','DTMax','FMin','FMax'},...
                'VariableNames',{'Object','IsDT','IsF','Hit'});
            obj.reset();
        end
        
        % setHighlight ----------------------------------------------------
        function setHighlight(obj,h,type,iSig,iScreen)
        % Highlights a marker
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            obj.hMarker = h;
            obj.signalType = type;
            obj.signalIdx = iSig;
            obj.screenIdx = iScreen;
            obj.drawLines();
        end
        
        % reset -----------------------------------------------------------
        function reset(obj)
        % Clears highlighting on a marker (if any) and resets properties to
        % placeholder values
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if obj.isDisplayed
                delete(obj.lineData.Object);
            end
            obj.hMarker = matlab.graphics.GraphicsPlaceholder;
            obj.signalType = '';
            obj.signalIdx = NaN;
            obj.screenIdx = NaN;
            obj.lineData.Object(:) =  matlab.graphics.GraphicsPlaceholder;
            obj.lineData.Hit(:) = false;
            %obj.reset();
        end
        
        % refresh ---------------------------------------------------------
        function refresh(obj,varargin)
        % checks if object has a valid marker for highlighting. If so,
        % lines are redrawn; if not, lines are deleted and the status of
        % the object gets updated to reflect this. This function can also 
        % take input arguments that allow updating the marker handle and 
        % screen index.
        %
        % Optional input arguments: 
        %   - hMarker
        %   - screenIdx
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            narginchk(1,3)
        
            % parse input if any
            if nargin > 1
                obj.hMarker = varargin{1};
                obj.screenIdx = varargin{2};
            end
            
            % check marker validity
            if isgraphics(obj.hMarker)
                obj.drawLines();
            else
                type = obj.signalType;
                iSig = obj.signalIdx;
                obj.reset();
                obj.signalType = type;
                obj.signalIdx = iSig;
            end
        end
        
        % checkMouseHits --------------------------------------------------
        function checkMouseHits(obj,x,y)
        % Checks if cursor position is over draggable edges, and updates
        % them accordingly
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            if obj.isDisplayed
                % extract useful data
                ha = obj.hAxes;
                ruler = ha.XAxis;
                DT1 = obj.DTMin;
                DT2 = obj.DTMax;
                F1 = obj.FMin;
                F2 = obj.FMax;
                editF = obj.editFreq;
                pBuff = obj.bufferProp;
                
                % convert x to datetime
                xDT = num2ruler(x,ruler);
                
                % set buffer size
                DTBuffer = diff(ha.XLim)*pBuff;
                FBuffer = diff(ha.YLim)*pBuff;

                % check if cursor has a hit with each highlight line
                withinDT = (xDT >= DT1 - DTBuffer) && (xDT <= DT2 + DTBuffer);
                withinF = (y >= F1 - FBuffer) && (y <= F2 + FBuffer);
                
                obj.lineData.Hit('DTMax') = withinF && (xDT >= DT2 - DTBuffer) && (xDT <= DT2 + DTBuffer);
                obj.lineData.Hit('DTMin') = ~obj.lineData.Hit('DTMax') && withinF && (xDT >= DT1 - DTBuffer) && (xDT <= DT1 + DTBuffer);
                obj.lineData.Hit('FMax') = editF && withinDT && (y >= F2 - FBuffer) && (y <= F2 + FBuffer);
                obj.lineData.Hit('FMin') = editF && ~obj.lineData.Hit('FMax') && withinDT && (y >= F1 - FBuffer) && (y <= F1 + FBuffer);
                %%% Thought: another solution might be to check mouse
                %%% distance from each line, and only check hits for the
                %%% DT/F lines that are closest
                
                obj.updateLineAppearance();
            end
        end
        
        % resize ----------------------------------------------------------
        function resize(obj,x,y)
        % Sets new DT and F limits for the highlighted marker
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            l = obj.lineData;

            % convert x to datetime
            ruler = obj.hAxes.XAxis;
            xDT = num2ruler(x,ruler);

            % check for DT hits
            DTRangeNew = [obj.DTMin,obj.DTMax];
            DTHit = l.Hit(l.IsDT);
            if any(DTHit)
                DTRangeNew(DTHit) = xDT;
                if diff(DTRangeNew) < 0
                    DTRangeNew = sort(DTRangeNew);
                    obj.swapEdgeHit('DT')
                end
            end
            
            % check for F hits
            FRangeNew = [obj.FMin,obj.FMax];
            FHit = l.Hit(l.IsF);
            if any(FHit)
                FRangeNew(FHit) = y;
                if diff(FRangeNew) < 0
                    FRangeNew = sort(FRangeNew);
                    obj.swapEdgeHit('F')
                end
            end

            % change line and patch coordinates
            set(l.Object('DTMin'),{'XData','YData'},{DTRangeNew([1,1]),FRangeNew([1,2])});
            set(l.Object('DTMax'),{'XData','YData'},{DTRangeNew([2,2]),FRangeNew([1,2])});
            set(l.Object('FMin'),{'XData','YData'},{DTRangeNew([1,2]),FRangeNew([1,1])});
            set(l.Object('FMax'),{'XData','YData'},{DTRangeNew([1,2]),FRangeNew([2,2])});
            markerXDataNew = ruler2num(DTRangeNew([1,1,2,2]),ruler);
            markerYDataNew = FRangeNew([1,2,2,1]);
            set(obj.hMarker,{'XData','YData'},{markerXDataNew,markerYDataNew});
        end
    end
    
    %% METHODS - GETTERS (public) =========================================
    methods
        % get.hasSelection ------------------------------------------------
        function s = get.hasSelection(obj)
                s = ~isnan(obj.signalIdx);
        end
        
        % get.isValid -----------------------------------------------------
        function s = get.isValid(obj)
            s = isgraphics(obj.hMarker);
        end
        
        % get.isDisplayed -------------------------------------------------
        function s = get.isDisplayed(obj)
            s = all(isgraphics(obj.lineData.Object));
        end
        
        % get.editFreq ----------------------------------------------------
        function e = get.editFreq(obj)
            e = strcmp(obj.signalType,'Annotations');
        end
        
        % get.DTMin -------------------------------------------------------
        function v = get.DTMin(obj)
            try
                v = obj.lineData.Object('DTMin').XData(1);
            catch
                v = NaN;
            end
        end
        
        % get.DTMax -------------------------------------------------------
        function v = get.DTMax(obj)
            try
                v = obj.lineData.Object('DTMax').XData(1);
            catch
                v = NaN;
            end
        end
        
        % get.FMin --------------------------------------------------------
        function v = get.FMin(obj)
            try
                v = obj.lineData.Object('FMin').YData(1);
            catch
                v = NaN;
            end
        end
        
        % get.FMax --------------------------------------------------------
        function v = get.FMax(obj)
            try
                v = obj.lineData.Object('FMax').YData(1);
            catch
                v = NaN;
            end
        end
    end
    
    %% METHODS - PRIVATE ==================================================
    methods (Access = private)
        % drawLines -------------------------------------------------------
        function drawLines(obj)
        % Draws the highlighting lines
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if obj.isValid
                % extract useful data
                ha = obj.hAxes;
                hp = obj.hMarker;
                ruler = ha.XAxis;
                
                % clear existing lines if there are any
                if obj.isDisplayed
                    delete(obj.lineData.Object);
                end

                % get DT and F limits of marker to be highlighted
                DTRange = num2ruler([min(hp.XData), max(hp.XData)], ruler);
                FRange = [min(hp.YData), max(hp.YData)];

                % plot lines
                [lineStyle, lineWidth, lineCol] = obj.getEdgeSpecs('Normal');
                hLineDTMin = line(ha, DTRange([1,1]), FRange([1,2]),'LineStyle',lineStyle,'LineWidth',lineWidth,'Color',lineCol);
                hLineDTMax = line(ha, DTRange([2,2]), FRange([1,2]),'LineStyle',lineStyle,'LineWidth',lineWidth,'Color',lineCol);
                hLineFMin  = line(ha, DTRange([1,2]), FRange([1,1]),'LineStyle',lineStyle,'LineWidth',lineWidth,'Color',lineCol);
                hLineFMax  = line(ha, DTRange([1,2]), FRange([2,2]),'LineStyle',lineStyle,'LineWidth',lineWidth,'Color',lineCol);

                % assign lines to object
                obj.lineData.Object = [hLineDTMin; hLineDTMax; hLineFMin; hLineFMax];
            else
                warning('Nothing to draw!')
            end
        end
        
        % updateLineAppearance --------------------------------------------
        function updateLineAppearance(obj)
        % updates line appearance based on hit status.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            l = obj.lineData;
            [lineStyleHit,lineWidthHit,colHit] = obj.getEdgeSpecs('Hover');
            [lineStyleNoHit,lineWidthNoHit,colNoHit] = obj.getEdgeSpecs('Normal');
            set(l.Object(l.Hit),{'LineStyle','LineWidth','Color'},{lineStyleHit,lineWidthHit,colHit});
            set(l.Object(~l.Hit),{'LineStyle','LineWidth','Color'},{lineStyleNoHit,lineWidthNoHit,colNoHit});
        end
        
        % swapEdgeHit -----------------------------------------------------
        function swapEdgeHit(obj,edgeType)
        % swaps the edge hit status of the Min vs. Max lines of a
        % specified axis
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            isType = obj.lineData.(['Is',edgeType]);
            obj.lineData.Hit(isType) = ~obj.lineData.Hit(isType);
            obj.updateLineAppearance();
        end
        
        % getEdgeSpecs ----------------------------------------------------
        function [lineStyle,lineWidth,col] = getEdgeSpecs(obj,type)
        % retrieves line specs for Normal or Hover mode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            % extract useful data
            hp = obj.hMarker;
            params = obj.edgeParams(type,:);
            
            % get specs
            lineStyle = params.LineStyle{:};
            lineWidth = hp.LineWidth*params.LineWidthMultiplier;
            col = min([hp.EdgeColor*params.ColorMultiplier; 1 1 1]);
        end
    end
end