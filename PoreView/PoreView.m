classdef PoreView < handle
    %POREVIEW: Analysis suite for streaming signal data
    %
    % ---------------------------------------------------------------------
    %
    % POREVIEW Usage:
    %   Adjusting view
    %       - Scrolling over a plot zooms the x-axis about the cursor position. 
    %       - Scrolling over a y-axis zooms that y axis about the cursor position.
    %       - Middle-click and drag pans the plot.
    %       - Right-clicking brings up the signal context menu.
    %       - Left-click and drag on axes to zoom.
    %       - Press 'a' to autoscale all axes.
    %
    %   Cursors
    %       - Double-click to bring cursors.
    %       - Click and drag cursors to move them.
    %       - Press 'c' to show/hide cursors.
    %
    %   Press 'Escape' to quit at any time.
    %
    % ---------------------------------------------------------------------
    %
    % POREVIEW Methods:
    %   PoreView(fname) - Starts PoreView. Filename can be empty, a
    %       filename, or a directory where your files are stored.
    %   loadFile(fname) - Loads a file, if it can find it. If IV curve is
    %       selected, forwards to plot_iv. Creates default signal panels.
    %
    %   setKeyboardCallback(fun) - Sets the keyboard callback function for
    %       user code. fun should take one argument, a struct, containing
    %       the key information.
    %   waitKey() - Blocks until a key is pressed. Returns character of the
    %       key (not a struct), eg 'k'.
    %
    %   addSignalPanel(sigs) - Add a signal panel, at the bottom, that
    %       displays the signals specified in sigs. Can be [].
    %   removeSignalPanel(panel) - Remove signal panel indexed by panel
    %   setSignalPanel(panel,sigs) - Sets panel to display signals sigs.
    %
    %   getAxes(pan) - Get the axes object handle for a given panel
    %   clearAxes() - Clear all user-drawn objects from all axes
    %
    %   autoscaleY() - Rescale Y-axes on all panels
    %
    %   getCursors() - Return cursor positions, or [] if hidden
    %   setCursors(trange) - Sets cursor positions and makes them visible
    %   toggleCursors() - Toggles visibility of cursors
    %
    %   refresh() - Redraws all plot displays.
    %   getView() - Returns time range visible in window.
    %   setView(trange) - Sets visible time range (clipping appropriately) and redraws.
    %
    % ---------------------------------------------------------------------
    %
    % POREVIEW Properties:
    %   data - Internal SignalData class, or [] if not loaded
    %   fig - Handle to figure object of program
    %   psigs - Struct array with signal panel information. Click for more.
    %
    % ---------------------------------------------------------------------
    %
    %
    %
    
    properties
        data        % Internal SignalData class, or [] if not loaded
        fig         % Handle to figure object
        % psigs - Struct array containing information about signal panels,
        %       as well as some internally defined helper functions.
        %   
        %   psigs(i).sigs - List of signals to display.
        %   psigs(i).axes - Axes object displaying signals
        %   psigs(i).setY(yrange) - Sets y-axis.
        %   psigs(i).getY() - Gets y-axis.
        %   psigs(i).resetY() - Autoscales y intelligently to look nice.
        %   psigs(i).shiftY(scale,offset) - Scales and moves y-axis. Offset
        %       is in units of the current y-range.
        %   
        %   There are some others that should not be important.
        psigs
    end
    properties (Hidden=true)
        DEFS        % UI definitions (widths etc.)
        
        panels      % Main window panel handle struct
        xaxes       % x-axis object
        cursors     % Main window cursor object (blue arrow thingies)
        ctext       % Cursor display text handle

        hcmenu      % Handle to context menu
        sigmenus    % Handles to all the add/remove signal menus
        cursig      % Most recent signal clicked
    end
    
    methods
        function obj = PoreView(fname)
            % pv = POREVIEW() - Create a new default instance of PoreView
            % pv = POREVIEW(filename) - Start by loading specified file
            % pv = POREVIEW(directory) - Start in a directory (for open file dialog)
            
            % some UI defs
            obj.DEFS = [];
            obj.DEFS.LEFTWID        = 60;
            obj.DEFS.BOTHEIGHT      = 60;
            obj.DEFS.BUTWID         = 20;
            obj.DEFS.BUTLEFT        = 3;
            obj.DEFS.BUTBOT         = 3;
            obj.DEFS.LABELWID       = 200;
            obj.DEFS.CURSCOLOR      = [0.5 0.3 0.0];
            
            % start making GUI objects
            obj.fig = figure('Name','PoreView','MenuBar','none',...
                'NumberTitle','off','DockControls','off');
            
            % set its position
            oldunits = get(obj.fig,'Units');
            set(obj.fig,'Units','normalized');
            set(obj.fig,'Position',[0.1,0.1,0.8,0.8]);
            set(obj.fig,'Units',oldunits);
            
            % Build a new color map
            CO = [       0    0.4470    0.7410;
                    0.8500    0.3250    0.0980;
                    0.9290    0.6940    0.1250;
                    0.4940    0.1840    0.5560;
                    0.4660    0.6740    0.1880;
                    0.3010    0.7450    0.9330;
                    0.6350    0.0780    0.1840];
                
            % Set the color order 
            set(obj.fig, 'DefaultAxesColorOrder', CO);
            
            if nargin == 0
                fname = '';
            end
            
            % how we load a new file
            function openFileFcn(~,~)
                % check if we've loaded a file yet
                if ~isempty(obj.data)
                    % if we did, let's use that filename
                    fn = obj.data.filename;
                else
                    % were we passed a directory, not file?
                    fn = fname;
                end
                % get a filename from dialog box
                [FileName,PathName] = uigetfile('*.abf;*.cbf;*.fast5','PoreView',fn);
                % and load (or attempt to)
                obj.loadFile([PathName FileName]);
            end
            
            % make the menu bar
            f = uimenu('Label','File');
            uimenu(f,'Label','Open','Callback',@openFileFcn);
            uimenu(f,'Label','Quit','Callback',@(~,~) close(obj.fig));
            c = uimenu('Label','Cursors');
            uimenu(c,'Label','Bring','Callback',@(~,~) obj.bringCursors());
            uimenu(c,'Label','Show/Hide','Callback',@(~,~) obj.toggleCursors());
            hm = uimenu('Label','Help');
            uimenu(hm,'Label','PoreView','Callback',@(~,~) doc('PoreView.m'));
            uimenu(hm,'Label','SignalData','Callback',@(~,~) doc('SignalData.m'));
            uimenu(hm,'Label','About','Callback',@(~,~) msgbox({'PoreView v1.0 - written by Tamas Szalay, April 2014.' '' ...
                'This program and its author are not affiliated with Molecular Devices.' ''},'About PoreView'));
            
            % and the context menus
            obj.hcmenu = uicontextmenu();
            obj.sigmenus = [];
            
            sa = uimenu(obj.hcmenu,'Label','Add Signal');
            st = uimenu(obj.hcmenu,'Label','Set Signal');
            snew = uimenu(obj.hcmenu,'Label','Add Panel','Separator','on',...
                'Callback',@(~,~) obj.addSignalPanel(obj.psigs(obj.cursig).sigs));
            srem = uimenu(obj.hcmenu,'Label','Remove Panel',...
                'Callback',@(~,~) obj.removeSignalPanel(obj.cursig));
            
            % populate signals submenu function
            function addsig(sig)
                % append sig to the list of this signal panel's signals
                obj.psigs(obj.cursig).sigs = [obj.psigs(obj.cursig).sigs sig];
                % and redraw
                obj.refresh();
            end
            function popSigMenuFcn(~,~)
                % delete the old ones
                delete(obj.sigmenus);
                if isempty(obj.data)
                    return
                end
                % get signal list, if we have one
                slist = obj.data.getSignalList();
                
                obj.sigmenus = [];
                % create the menus
                for i=2:length(slist)
                    obj.sigmenus(end+1) = uimenu(sa,'Label',slist{i});
                    set(obj.sigmenus(end),'Callback',@(~,~) addsig(i));
                end
                
                % and to set
                for i=2:length(slist)
                    obj.sigmenus(end+1) = uimenu(st,'Label',slist{i});
                    set(obj.sigmenus(end),'Callback',@(~,~) obj.setSignalPanel(obj.cursig,i));
                end
            end
            
            set(obj.hcmenu,'Callback',@popSigMenuFcn);
            
            obj.cursig = 0;
            
            % the main layout components
            obj.panels = [];
            obj.panels.Middle = uipanel('Parent',obj.fig,'Position',[0 0.5 1 0.5],'Units','Pixels');
            obj.panels.Bottom = uipanel('Parent',obj.fig,'Position',[0 0.5 1 0.5],'Units','Pixels');
            
            set(obj.panels.Bottom,'BorderType','none');
            set(obj.panels.Middle,'BorderType','none');
            
            % handles the resizing of the main panels
            function mainResizeFcn(~,~)
                sz = getPixelPos(obj.fig);

                set(obj.panels.Bottom,'Position',[1,0,sz(3)+2,obj.DEFS.BOTHEIGHT]);
                % set this guy outside the edges of the figure by one pixel
                % horizontally, to hide the border on the sides
                set(obj.panels.Middle,'Position',[0,obj.DEFS.BOTHEIGHT,sz(3)+2,sz(4)-obj.DEFS.BOTHEIGHT]);
            end
            set(obj.fig,'ResizeFcn',@mainResizeFcn);
            mainResizeFcn
            
            % ========== X AXIS CODE ===========
            % make an x-axis for display porpoises only. this will need to
            % get resized correctly later, sadly :-/
            hxaxes = axes('Parent',obj.panels.Bottom,'TickDir','out',...
                'Position',[0 1 1 0.01],'YTickLabel','',...
                'YLimMode','manual','XLimMode','manual');
            
            obj.xaxes = hxaxes;
            
            % create dummy cursor lines in the x-axis that the real
            % cursor objects can copy
            % and, conveniently, use them to display little arrows
            obj.cursors = [line() line()];
            set(obj.cursors,'Parent',obj.xaxes,'XData',[0 0],...
                'YData',[-3 -3],'Visible','off','Tag','PVCURS',...
                'Color',obj.DEFS.CURSCOLOR,'Marker','^',...
                'MarkerFaceColor',obj.DEFS.CURSCOLOR,'MarkerSize',5,...
                'Clipping','off','LineStyle','none');
            
            
            % again, screw scroll bars
            function shiftX(zoom,offset)
                xlim = get(hxaxes,'XLim');
                dx = xlim(2)-xlim(1);
                xm = mean(xlim);
                xlim = xm+zoom*(xlim-xm);
                xlim = xlim+dx*offset;
                obj.setView(xlim);
            end
            
            % now make the buttons
            nbut = 5;
            buts = zeros(nbut,1);
            
            buts(1) = uicontrol('Parent', obj.panels.Bottom, 'String','<html>-</html>',...
                'callback', @(~,~) shiftX(2,0));
            buts(2) = uicontrol('Parent', obj.panels.Bottom, 'String','<html>&larr;</html>',...
                'callback', @(~,~) shiftX(1,-0.25));
            buts(3) = uicontrol('Parent', obj.panels.Bottom, 'String','<html>R</html>',...
                'callback', @(~,~) obj.setView());
            buts(4) = uicontrol('Parent', obj.panels.Bottom, 'String','<html>&rarr;</html>',...
                'callback', @(~,~) shiftX(1,0.25));
            buts(5) = uicontrol('Parent', obj.panels.Bottom, 'String','<html>+</html>',...
                'callback', @(~,~) shiftX(0.5,0));
            
            % and create a text-box in the bottom-right to hold the dt
            % (time between cursors)
            obj.ctext = uicontrol('Parent',obj.panels.Bottom,...
                'Style','text','Units','Pixels','Position',[0 0 1 1],...
                'Tag','Cursor dt: %0.3f %s   ','FontSize',10,...
                'HorizontalAlignment','right','Visible','off');
            uistack(obj.ctext,'bottom');
            l = linkprop([obj.ctext obj.cursors(1)],'Visible');
            set(obj.ctext,'UserData',l);
            
            % how to move buttons when thingy gets resized
            function resizeFcn(~,~)
                % get height of panel in pixels
                sz = getPixelPos(obj.panels.Bottom);
                % figure out where the middle is
                mid = sz(3)/2;
                for i=1:nbut
                    % position the buttons
                    set(buts(i),'Position',[mid+(i-nbut/2-1)*obj.DEFS.BUTWID,obj.DEFS.BUTBOT,obj.DEFS.BUTWID,obj.DEFS.BUTWID]);
                end
                % also need to resize x-axis labels
                set(obj.xaxes,'Units','Pixels');
                set(obj.xaxes,'Position',[sz(1)+obj.DEFS.LEFTWID,sz(4)+2,sz(3)-obj.DEFS.LEFTWID,1]);
                % and tick length
                s = 6/max(sz(3:4));
                set(obj.xaxes,'TickLength',s*[1 1]);
                
                % now position the stupid label thing
                set(obj.ctext,'Position',[sz(3)-obj.DEFS.LABELWID 1 obj.DEFS.LABELWID 20]);
            end
            % set the resize function
            set(obj.panels.Bottom, 'ResizeFcn', @resizeFcn);
            % and call it to set default positions
            resizeFcn
            
            
            
            % ========== FINISHING TOUCHES ===========
            
            % create drag/drop functions, etc
            obj.setMouseCallbacks();
            % and a dummy keyboard callback
            obj.setKeyboardCallback(@(e) []);
            
            % create a blank signal panel and blank data
            obj.data = [];
            obj.addSignalPanel([]);
            
            % load data if we were called with a filename
            % this also creates default signal panels
            if (nargin > 0)
                obj.loadFile(fname);
            end
        end
        
    end
    
    methods (Hidden = true)
        function setMouseCallbacks(obj)
            % obj.SETMOUSECALLBACKS()
            %   Creates mouse callback interface, by defining a ton of fns
            
            % point-rectangle hit test utility function
            function b = isIn(pos,p)
                b = false;
                if (pos(1) > p(1) && pos(1) < (p(1)+p(3)) &&...
                            pos(2) > p(2) && pos(2) < (p(2)+p(4)))
                    b = true;
                end
            end
            % function that steps through all relevant objects to figure
            % out which one is at a given position
            function [hnd,ind,pt,s] = getHandleAt(pos)
                % pos should be in pixels
                if nargin < 1
                    pos = get(obj.fig,'CurrentPoint');
                end
                % first, let's check signals
                nsig = length(obj.psigs);
                for i=1:nsig
                    ind = i;
                    % did we click on a plot?
                    if isIn(pos,getpixelposition(obj.psigs(i).axes,true))
                        hnd = obj.psigs(i).axes;
                        pt = get(hnd,'CurrentPoint');
                        pt = pt(1,1:2);
                        s = 'a';
                        obj.cursig = ind;
                        return;
                    end
                    % are we over the y axis part?
                    if (isIn(pos,getpixelposition(obj.psigs(i).panel,true)) &&...
                            pos(1) > 0.6*obj.DEFS.LEFTWID)
                        % return the handle to the main axes anyway, s'more
                        % useful that way
                        hnd = obj.psigs(i).axes;
                        pt = get(hnd,'CurrentPoint');
                        pt = pt(1,1:2);
                        s = 'y';
                        obj.cursig = ind;
                        return;
                    end
                end
                % ok, not signals, so let's check x-axis
                if (isIn(pos,getpixelposition(obj.panels.Bottom)) &&...
                        pos(2) > 0.6*obj.DEFS.BOTHEIGHT)
                    hnd = obj.xaxes;
                    ind = -1;
                    pt = get(hnd,'CurrentPoint');
                    pt = pt(1,1:2);
                    s = 'x';
                    obj.cursig = 0;
                    return;
                end
                hnd = -1;
                ind = -1;
                pt = [];
                s = '';
                obj.cursig = 0;
            end
            % handles scrolling zoom in-out
            function scrollCallback(~,e)
                % figure out what we are over, if anything
                [hnd,ind,pt,s] = getHandleAt();
                
                % get the scaling factor to zoom by
                scl = 1.4 .^ (e.VerticalScrollCount);
                
                % quit if found nothing
                if (hnd==-1)
                    return
                end
                
                if (s == 'y')
                    % we're scrolling inside a y-axis, scroll y-lims
                    pty = pt(2);
                    ylim = obj.psigs(ind).getY();
                    ylim = sort(pty + scl*(ylim-pty));
                    obj.psigs(ind).setY(ylim);
                elseif (s == 'a' || s == 'x')
                    % we're scrolling in a plot, scroll the time axis
                    ptx = pt(1);
                    xlim = get(obj.xaxes,'XLim');
                    xlim = sort(ptx + scl*(xlim-ptx));
                    obj.setView(xlim);
                end
            end
            
            function mouseMoveX(r)
                % get start and current point
                pt0 = get(obj.xaxes,'UserData');
                pt1 = get(obj.xaxes,'CurrentPoint');
                % convert to x-range
                xr = sort([pt0(1) pt1(1,1)]);
                % update the line
                if ishandle(r)
                    set(r,'YData',[0,0],'XData',xr,'LineWidth',8);
                end
            end
            function mouseMoveY(r,ind)
                % get start and current point
                pt0 = get(obj.psigs(ind).yaxes,'UserData');
                pt1 = get(obj.psigs(ind).yaxes,'CurrentPoint');
                % y-range, sorted
                yr = sort([pt0(2) pt1(1,2)]);
                % set x limits of y-line
                xl = get(obj.psigs(ind).yaxes,'XLim');
                % and update the line, if it exists
                if ishandle(r)
                    set(r,'XData',[xl(1),xl(1)],'YData',yr,'LineWidth',5);
                end
            end
            function mouseUpX(r)
                % kill the line
                delete(r);
                % get the x-range one last time
                pt0 = get(obj.xaxes,'UserData');
                pt1 = get(obj.xaxes,'CurrentPoint');
                xrange = sort([pt0(1) pt1(1,1)]);
                % make sure we zoomed at all
                if xrange(1)==xrange(2)
                    return
                end
                % then update the view
                obj.setView(xrange);
                % and clear callbacks
                set(obj.fig,'WindowButtonUpFcn','');
                set(obj.fig,'WindowButtonMotionFcn','');
            end
            function mouseUpY(r,ind)
                % kill the line
                delete(r);
                % y-range as usual
                pt0 = get(obj.psigs(ind).yaxes,'UserData');
                pt1 = get(obj.psigs(ind).yaxes,'CurrentPoint');
                yrange = sort([pt0(2) pt1(1,2)]);
                if yrange(1)==yrange(2)
                    return
                end
                % zoom in
                obj.psigs(ind).setY(yrange);
                % and clear callbacks
                set(obj.fig,'WindowButtonUpFcn','');
                set(obj.fig,'WindowButtonMotionFcn','');
            end
            function mouseMovePan(ind)
                % this is the doozy
                pt0 = get(obj.psigs(ind).yaxes,'UserData');
                pt1 = get(obj.psigs(ind).axes,'CurrentPoint');
                pt1 = pt1(1,1:2);
                % so we want to get pt0 and pt1 equal
                % to get that guy under the mouse
                dpt = pt1-pt0;
                % and x axis
                xl = get(obj.xaxes,'XLim');
                obj.setView(xl - dpt(1));
                % now update y-axis
                yl = obj.psigs(ind).getY();
                obj.psigs(ind).setY(yl-dpt(2));
            end
            function mouseUpPan()
                % clear callbacks
                set(obj.fig,'WindowButtonUpFcn','');
                set(obj.fig,'WindowButtonMotionFcn','');
            end
            % and create their mouse down callbacks
            function mouseMoveCursor(ind,cursind)
                % get dragged point
                pt = get(obj.psigs(ind).axes,'CurrentPoint');
                % now update cursor positions
                curs = obj.getCursors();
                curs(cursind) = pt(1);
                obj.setCursors(curs);
            end
            function mouseUpCursor()
                % clear callbacks
                set(obj.fig,'WindowButtonUpFcn','');
                set(obj.fig,'WindowButtonMotionFcn','');
            end
            % the main mouse down dispatch function
            function mouseDown(~,~)
                % figure out what we are over, if anything
                [hnd,ind,pt,s] = getHandleAt();
                
                sel = get(obj.fig,'SelectionType');
                if (hnd == -1)
                    return
                end
                
                % see if we are dragging on x or y axes
                if (s == 'x' && strcmp(sel,'normal'))
                    % make a line
                    r = line();
                    set(r,'Parent',obj.xaxes);
                    % store the drag start point
                    set(obj.xaxes,'UserData',pt);
                    % and set appropriate callbacks, with line as data
                    set(obj.fig,'WindowButtonUpFcn',@(~,~) mouseUpX(r));
                    set(obj.fig,'WindowButtonMotionFcn',@(~,~) mouseMoveX(r));
                elseif (s == 'y' && strcmp(sel,'normal'))
                    % make a line
                    r = line();
                    set(r,'Parent',obj.psigs(ind).yaxes);
                    % store the drag start point
                    set(obj.psigs(ind).yaxes,'UserData',pt);
                    % and set appropriate callbacks, with line as data
                    set(obj.fig,'WindowButtonUpFcn',@(~,~) mouseUpY(r, ind));
                    set(obj.fig,'WindowButtonMotionFcn',@(~,~) mouseMoveY(r, ind));
                elseif (s == 'a' && strcmp(sel,'extend'))
                    % store the drag start point
                    set(obj.psigs(ind).yaxes,'UserData',pt);
                    % and set the callbacks
                    set(obj.fig,'WindowButtonUpFcn',@(~,~) mouseUpPan());
                    set(obj.fig,'WindowButtonMotionFcn',@(~,~) mouseMovePan(ind));
                elseif (s == 'a' && strcmp(sel,'open'))
                    obj.bringCursors();
                elseif (s == 'a' && strcmp(sel,'normal'))
                    % click-drag a cursor?
                    cursind = find(get(obj.fig,'CurrentObject') == obj.psigs(ind).cursors,1);
                    if ~isempty(cursind)
                         % set the callbacks
                        set(obj.fig,'WindowButtonUpFcn',@(~,~) mouseUpCursor());
                        set(obj.fig,'WindowButtonMotionFcn',@(~,~) mouseMoveCursor(ind,cursind));
                    end
                end
            end
            
            % set global figure callbacks
            set(obj.fig,'WindowScrollWheelFcn',@scrollCallback);
            set(obj.fig,'WindowButtonDownFcn',@mouseDown);
        end
    end     
        
        
        
        
        
        
    methods (Hidden = false)

        function loadFile(obj, fname)
            % obj.LOADFILE(filename)
            %   Loads the specified file, and initializes the signal panels
            
            % attempt to load data
            d = SignalData(fname);
            
            if d.ndata < 0
                if (d.ndata == -2)
                    disp('IV curve found, attempting load...');
                    plot_iv(fname);
                end
                return
            end
            
            % get last bit of filename
            [~,fn,ext] = fileparts(fname);
            % and set title bar
            set(obj.fig,'Name',['PoreView - ' fn ext]);
            
            % set internal data
            obj.data = d;
            
            % first, delete all panels
            while ~isempty(obj.psigs)
                obj.removeSignalPanel();
            end
            % then, create right number of panels
            for i=1:obj.data.nsigs
                obj.addSignalPanel(i+1);
            end
            
            % and reset view
            obj.setView();
            
            % and cursors
            obj.setCursors(obj.getView());
            obj.toggleCursors();
        end
        
        function setKeyboardCallback(obj, fun)
            % obj.SETKEYBOARDCALLBACK(fun)
            %   PoreView will call fun whenever a key is pressed (except
            %   for default keys). The function should take a single
            %   argument, which is a Matlab keyboard event struct with
            %   fields e.Characer, e.Key, e.Modifier.
            
            % define our own internal keyboard function
            function keyboardFcn(~,e)
                if e.Character == 'a'
                    % autoscale axes
                    for i=1:length(obj.psigs)
                        obj.psigs(i).resetY();
                    end
                    return
                end
                if e.Character == 'c'
                    obj.toggleCursors();
                    return
                end
                if strcmp(e.Key,'escape')
                    close(obj.fig);
                    return
                end
                
                % just in case we're waiting for a blocking keyboard event
                % we should make sure it happens
                uiresume(obj.fig);
                
                % and then pass it on to the user-defined function
                % (but only if a real key was pressed, not just shift etc)
                if ~any(strcmp(e.Key,{'shift','control','alt'}))
                    fun(e);
                end
            end
            set(obj.fig,'WindowKeyPressFcn',@keyboardFcn);
        end
        
        function k = waitKey(obj)
            % k = obj.WAITKEY()
            %   Blocks until a (non-builtin) key is pressed. Returns the
            %   key character only (eg. 'k').
            
            % wait for uiresume event
            uiwait(obj.fig);
            
            % key was pressed, fetch it
            k = get(obj.fig,'CurrentCharacter');
        end
        
        function addSignalPanel(obj, sigs)
            % obj.addSignalPanel(sigs)
            %   Add a new signal panel at the bottom with the specified
            %   signals (can be []).
            
            % make a new one at the end
            i = length(obj.psigs)+1;
            if (i==1)
                obj.psigs = obj.makeSignalPanel();
            else
                obj.psigs(i) = obj.makeSignalPanel();
            end
            % save which signals we want to draw
            obj.psigs(i).sigs = sigs;
            % and redo their sizes
            obj.sizeSignalPanels()
            % and refresh the view, without changing x-limits
            if i==1
                % or just reset, if it's the first panel
                % (which it really shouldn't be?)
                obj.setView();
            else
                obj.refresh();
            end
            % autoset the y-limits
            obj.psigs(i).resetY();
        end
        
        function removeSignalPanel(obj, pan)
            % obj.removeSignalPanel()
            % obj.removeSignalPanel(pan)
            %   Remove the specified signal panel, or the bottommost one.
            %   Cannot remove all of them.
            
            % just remove the topmost one if none specified
            if (nargin < 2)
                pan = 1;
            end
            % already deleted all of them
            if isempty(obj.psigs) || (nargin == 2 && length(obj.psigs)==1)
                return
            end
            % otherwise, delete the panel object
            delete(obj.psigs(pan).panel);
            % and remove it from the array
            obj.psigs = obj.psigs((1:length(obj.psigs))~=pan);
            % and resize!
            obj.sizeSignalPanels();
        end
        function setSignalPanel(obj, pan, sigs)
            % obj.setSignalPanel(pan,sigs)
            %   Set the specified signal panel to have signal sigs.
            
            obj.psigs(pan).sigs = sigs;
            obj.refresh();
        end
        
    end
    
    methods (Hidden=true)
        
        function sizeSignalPanels(obj)
            % obj.sizeSignalPanels()
            %   Set the sizes of the panels on the screen to be evenly
            %   split vertically. Internal function.
            np = length(obj.psigs);
            for i=1:np
                set(obj.psigs(i).panel,'Position',[0 (np-i)/np 1 1/np]);
            end
        end
        
        function sig = makeSignalPanel(obj)
            % sig = obj.makeSignalPanel()
            %   This makes a single panel. Positioning/tiling is taken
            %   care of externally after creation. Returns a struct
            %   containing handles and view functions.
        
            % create and ultimately return a struct containing panel info
            sig = [];
            
            % which data to view? (none by default)
            sig.sigs = [];
            
            % make a panel to hold the entire thing
            sig.panel = uipanel('Parent',obj.panels.Middle,'Position',[0 0 1 1],'Units','Normalized');
            % give it a stylish border
            set(sig.panel,'BorderType','line','BorderWidth',1);
            set(sig.panel,'HighlightColor','black');
            
            % first make a fake axes object that just displays the label
            sig.yaxes = axes('Parent',sig.panel,'Position',[0 0 1 1],...
                'TickDir','out','Box','off','XLimMode','manual');
            % and then the main axes object for showing data
            sig.axes = axes('Parent',sig.panel,'Position',[0 0 1 1],...
                'XTickLabel','','YTickLabel','','GridLineStyle','-',...
                'XColor', 0.85*[1 1 1],'YColor', 0.85*[1 1 1]);
            % let's bind the context menu here as well
            set(sig.axes,'uicontextmenu',obj.hcmenu);
            % equivalent to 'hold on'
            set(sig.axes,'NextPlot','add','XLimMode','manual');
            % and gridify it
            set(sig.axes,'XGrid','on','YGrid','on','Tag','PVAXES');

            % magic function to make Y-axes consistent, saving me some
            % bookkeeping headaches and stuff
            l = linkprop([sig.yaxes sig.axes],{'YLim'});
            set(sig.axes,'UserData',l);
            
            % screw scroll bars, you never need to scroll the current trace
            % so we'll just do buttons instead
            % the callback scales the y limits of things
            function ylim = getY()
                ylim = get(sig.axes,'YLim');
            end
            function setY(ylim)
                ylim = sort(ylim);
                if (ylim(2)-ylim(1)) < 1e-10
                    return
                end
                set(sig.axes,'YLimMode','manual');
                set(sig.axes,'YLim',ylim);
            end
            function resetY()
                % let Matlab scale our axes, but don't scale to user-added
                % objects or to data cursors (unless there is no data!)
                % so first we make them all invisible
                cvis = get(obj.cursors(1),'Visible');
                set(obj.cursors,'Visible','off');
                hs = findobj(sig.axes,'-not','Tag','PVPLOT');
                if (~isempty(obj.data))
                    % if no data, scale to user-plotted stuff
                    set(hs,'Visible','off');
                end
                set(sig.axes,'YLimMode','auto');
                ylim = getY();
                % and then make them visible again
                set(hs,'Visible','on');
                % and restore data cursor visibility to previous value
                set(obj.cursors,'Visible',cvis);
                setY(ylim);
                % zoom out a tiny bit to make it visually pleasing
                shiftY(1.25,0);
            end
            function shiftY(zoom,offset)
                ylim = getY();
                dy = ylim(2)-ylim(1);
                ym = mean(ylim);
                ylim = ym+zoom*(ylim-ym);
                ylim = ylim+dy*offset;
                
                setY(ylim);
            end
            % save the y-limit functions to our little struct
            sig.shiftY = @shiftY;
            sig.setY = @setY;
            sig.resetY = @resetY;
            sig.getY = @getY;
            
            
            % make the cursors, and start them off invisible
            sig.cursors = [line() line()];
            set(sig.cursors,'Parent',sig.axes,'XData',[5 5],'YData',1e3*[-1 1],...
                'Color',obj.DEFS.CURSCOLOR,'Visible','off','Tag','PVCURS');
            
            % link their properties to the dummy x-axis cursor
            % this way we never have to think about it
            for j=1:2
                l = linkprop([obj.cursors(j) sig.cursors(j)],{'XData','Visible'});
                set(sig.cursors(j),'UserData',l);
                set(sig.cursors(j),'Visible','off');
            end
            
            % now make the buttons
            nbut = 5;
            % store the handles in a little array
            buts = zeros(nbut,1);
            % buts isn't getting put in sig, because we don't need to
            % touch it from the outside...
            
            buts(1) = uicontrol('Parent', sig.panel, 'String','<html>-</html>',...
                'callback', @(~,~) sig.shiftY(2,0));
            buts(2) = uicontrol('Parent', sig.panel, 'String','<html>&darr;</html>',...
                'callback', @(~,~) sig.shiftY(1,-0.25));
            buts(3) = uicontrol('Parent', sig.panel, 'String','<html>A</html>',...
                'callback', @(~,~) sig.resetY());
            buts(4) = uicontrol('Parent', sig.panel, 'String','<html>&uarr;</html>',...
                'callback', @(~,~) sig.shiftY(1,0.25));
            buts(5) = uicontrol('Parent', sig.panel, 'String','<html>+</html>',...
                'callback', @(~,~) sig.shiftY(0.5,0));
            
            % how to move buttons when thingy gets resized
            function resizeFcn(~,~)
                % get height of panel in pixels
                sz = getPixelPos(sig.panel);
                % figure out where the middle is
                mid = sz(4)/2;
                for i=1:nbut
                    % position the buttons
                    set(buts(i),'Position',...
                        [obj.DEFS.BUTLEFT,mid+(i-nbut/2-1)*obj.DEFS.BUTWID,obj.DEFS.BUTWID,obj.DEFS.BUTWID]);
                end
                % now position the axes objects too
                set(sig.axes,'Units','Pixels');
                set(sig.yaxes,'Units','Pixels');
                sz(1) = sz(1) + obj.DEFS.LEFTWID;
                sz(3) = sz(3) - obj.DEFS.LEFTWID;
                % set bottom to 1
                sz(2) = 1;
                set(sig.axes,'Position',sz);
                sz(1) = sz(1) + 1;
                set(sig.yaxes,'Position',sz);
                % now do some ticklength calcs
                s = 5/max(sz(3:4));
                set(sig.yaxes,'TickLength',s*[1 1]);
            end
            % set the resize function
            set(sig.panel, 'ResizeFcn', @resizeFcn);
            % and call it to set default positions
            resizeFcn
        end
        
    end
    
    methods

        function h = getAxes(obj, ind)
            % h = obj.getAxes(panel)
            %   Returns a handle to the axes object on a particular signal
            %   panel, so user can draw on it.
            h = obj.psigs(ind).axes;
        end
        
        function clearAxes(obj)
            % obj.clearAxes()
            %   Clears all user-added objects from all axes.
            
            for i=1:length(obj.psigs)
                % find objects with empty tags, and kill them
                hs = findobj(obj.psigs(i).axes,'Tag','');
                delete(hs);
            end
        end
        
        function autoscaleY(obj)
            % obj.autoscaleY()
            %   Autoscales y-axes in all panels, ignoring user-drawn
            %   objects.
            
            for i=1:length(obj.psigs)
                obj.psigs(i).resetY();
            end
        end
        
        function bringCursors(obj)
            % obj.bringCursors()
            %   Brings cursors to the view window and makes them visible.
            
            xlim = obj.getView();
            xm = mean(xlim);
            xlim = xm+0.5*(xlim-xm);
            obj.setCursors(xlim);
        end
        
        function toggleCursors(obj)
            % obj.toggleCursors()
            %   Toggle cursor visibility
            
            cvis = get(obj.cursors(1),'Visible');
            if strcmp(cvis,'on')
                set(obj.cursors,'Visible','off');
            else
                set(obj.cursors,'Visible','on');
            end
        end
        
        function curs = getCursors(obj)
            % curs = obj.getCursors()
            %   Gets cursor positions, or [] if they are not visible.
            
            % return empty if they're invisible
            if strcmp(get(obj.cursors(1),'Visible'),'off')
                curs = [];
                return
            end
            
            % otherwise return their positions
            curs = zeros(1,2);
            for i=1:2
                curs(i) = mean(get(obj.cursors(i),'XData'));
            end
            % return in proper time-order i guess
            curs = sort(curs);
        end 
        
        function setCursors(obj, curs)
            % obj.setCursors(curs)
            %   Set cursor positions and make them visible.

            curs = sort(curs);
            for i=1:2
                set(obj.cursors(i),'XData',[curs(i) curs(i)],'Visible','on');
            end
            % now update the dt label
            dt = curs(2) - curs(1);
            s = 's';
            if (dt < 1e-3)
                dt = dt * 1e6;
                s = 'us';
            elseif (dt < 1)
                dt = dt * 1e3;
                s = 'ms';
            end
            % the format string is hidden in ctext's tag
            set(obj.ctext,'String',sprintf(get(obj.ctext,'Tag'),dt,s));
        end 
        
        function refresh(obj)
            % obj.refresh()
            %   Just refresh the view and force a redraw.
            
            obj.setView(obj.getView());
        end
        
        function xlim=getView(obj)
            % xlim = obj.getView()
            %   Returns x-limits of current screen.
            
            xlim = get(obj.xaxes,'XLim');
        end
        
        function setView(obj,rng)
            % obj.setView(xlim)
            %   Sets the x-limits (with bounding, of course). Also
            %   reloads/redraws everything.
            
            % if we don't have anything loaded, just limit to plotted
            % objects et cetera
            if isempty(obj.data)
                if (nargin < 2)
                    % reset view
                    set([obj.psigs.axes],'XLimMode','auto');
                    % and find the max of them
                    mm = [inf -inf];
                    for i=1:numel(obj.psigs)
                        xlim = get(obj.psigs(i).axes,'XLim');
                        mm(1) = min(mm(1),xlim(1));
                        mm(2) = max(mm(2),xlim(2));
                    end
                    set(obj.xaxes,'XLimMode','manual');
                    set([obj.psigs.axes],'XLimMode','manual');
                    set(obj.xaxes,'XLim',mm);
                    set([obj.psigs.axes],'XLim',mm);
                else
                    % set view
                    set(obj.xaxes,'XLim',rng);
                    set([obj.psigs.axes],'XLim',rng);
                end
                return
            end
            
            % also, check if we got a range at all or not
            if (nargin < 2)
                range = [obj.data.tstart obj.data.tend];
            else
                range = rng;
            end
            
            % are we too zoomed-out?
            dr = range(2)-range(1);
            if (dr > obj.data.tend)
                range = [obj.data.tstart obj.data.tend];
                dr = obj.data.tend;
            end
            % or too zoomed-in?
            if (dr == 0)
                return
            end
            
            % are we too far to the left?
            if (range(1) < 0)
                % then shift back
                range = range - range(1);
            end
            % or too far to the right?
            if (range(2) > obj.data.tend)
                range = [-dr 0] + obj.data.tend;
            end
            
            % now we can set all the xlims properly
            set(obj.xaxes,'XLim',range);
            set([obj.psigs.axes],'XLim',range);
            
            % get the data, and whether we're looking at reduced
            [d, isred] = obj.data.getViewData(range);
            % set the axes color maps to be lighter if using reduced
            CO = get(obj.fig, 'DefaultAxesColorOrder');
            % normal CO is dark
            if isred
                % make it lighter
                CO = 1 - 0.80*(1-CO);
            end
            
            % and replot everything
            for i=1:length(obj.psigs)
                % now, don't clear everything (using cla)
                % instead, just delete the lines we drew previous
                delete(findobj(obj.psigs(i).axes,'Tag','PVPLOT'));
                
                % plot the selected signals, if any
                if isempty(obj.psigs(i).sigs)
                    continue
                end
                % set the axes color order
                set(obj.psigs(i).axes,'ColorOrder',CO);
                
                % get the plot handles
                hps = plot(obj.psigs(i).axes,d(:,1),d(:,obj.psigs(i).sigs));
                % and tag them to clear them next time, also make them
                % non-clickable
                set(hps,'Tag','PVPLOT','HitTest','off');
                % and move the plotted lines to the bottom of axes
                % this line slows it down a bit, but oh well...
                uistack(flipud(hps),'bottom');

                if (nargin < 2)
                    % if we did setView(), reset Y
                    obj.psigs(i).resetY();
                end
            end
        end
    end
    
end

% some useful helper functions
function sz = getPixelPos(hnd)
    old_units = get(hnd,'Units');
    set(hnd,'Units','Pixels');
    sz = get(hnd,'Position');
    set(hnd,'Units',old_units);
end