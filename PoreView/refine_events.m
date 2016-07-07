function DoubleEvents = refine_events(pv, DoubleEvents)

    pv.addSignalPanel([]);
    pv.addSignalPanel([]);

    % define the things we need for drawing and such
    selectedEvent = -1;
    
    function setSel(j)
        selectedEvent = j;
        updateEvent(j);
        redrawEvents();
    end
    function setEventView()
        if selectedEvent > 0 && selectedEvent <= numel(DoubleEvents)
            ev = DoubleEvents(selectedEvent);
            dt = abs(mean(ev.t1) - mean(ev.t0));
            if (dt < 0.002)
                dt = 0.002;
            end
            pv.setView([mean(ev.t1)-dt, mean(ev.t0)+dt]);            
            pv.autoscaleY();
            
            pv.setCursors(ev.t0);
        end
    end
    function redrawEvents()
        pv.clearAxes();
        
        for j=1:numel(DoubleEvents)
            ev = DoubleEvents(j);
            fs = 'g';
            if ev.t1 < 0
                fs = 'r';
            end
            if j==selectedEvent
                fs = [fs 'o-'];
            end
            
            pp = [0 0];
            pp(1) = plot(pv.getAxes(1),mean(ev.t0)*[1 1],[ev.m0,ev.m0+0.4],fs);
            pp(2) = plot(pv.getAxes(2),mean(ev.t1)*[1 1],[ev.m1,ev.m1-1],fs);
            set(pp,'ButtonDownFcn',@(~,~) setSel(j));
    
            % do we draw full event thingy? draw beginning and end
            if j==selectedEvent
                pv.psigs(3).setY([ev.peak0 ev.m0] +  0.1*[-1 1]);
                pv.psigs(4).setY([ev.peak1 ev.m1] +  0.1*[1 -1]);
                
                dt = 3e-3*[-1 1];
                
                dd = pv.data.getByTime(ev.t0 + dt);
                fake0 = gen_event(dd(:,1),ev.t0,ev.m0,ev.peak0-ev.m0);
                
                t0_calc = abs(ev.ecd0 / (ev.peak0-ev.m0));
                t0_calc = mean(ev.t0) + [-0.5 0.5]*t0_calc;
                fake0_calc = gen_event(dd(:,1),t0_calc,ev.m0,ev.peak0-ev.m0);

                plot(pv.getAxes(3),dd(:,1),fake0,'r');
                plot(pv.getAxes(3),dd(:,1),fake0_calc,'g');
                plot(pv.getAxes(3),dd(:,1),filt_event(dd(:,2),fake0));
                
                
                dd = pv.data.getByTime(ev.t1 + dt);
                fake1 = gen_event(dd(:,1),ev.t1,ev.m1,ev.peak1-ev.m1);
                
                t1_calc = abs(ev.ecd1 / (ev.peak1-ev.m1));
                t1_calc = mean(ev.t1) + [-0.5 0.5]*t1_calc;
                fake1_calc = gen_event(dd(:,1),t1_calc,ev.m1,ev.peak1-ev.m1);
                
                plot(pv.getAxes(4),dd(:,1),fake1,'r');
                plot(pv.getAxes(4),dd(:,1),fake1_calc,'g');
                plot(pv.getAxes(4),dd(:,1),filt_event(dd(:,3),fake1));
            end
        end
    end

    function updateEvent(ind)
        
        if (ind < 1 || ind > numel(DoubleEvents))
            return
        end
        ev = DoubleEvents(ind);
        
        dt = 3e-3*[-1 1];
                
        dd = pv.data.getByTime(ev.t0 + dt);
        fake0 = gen_event(dd(:,1),ev.t0,ev.m0,ev.peak0-ev.m0);
        filt0 = filt_event(dd(:,2),fake0) - ev.m0;
        inds = dd(:,1) >= ev.t0(1) & dd(:,1) <= ev.t0(2);
        ev.ecd0 = sum(filt0(inds))*pv.data.si;
        
        dd = pv.data.getByTime(ev.t1 + dt);
        fake1 = gen_event(dd(:,1),ev.t1,ev.m0,ev.peak1-ev.m1);
        filt1 = filt_event(dd(:,3),fake1) - ev.m1;
        inds = dd(:,1) >= ev.t1(1) & dd(:,1) <= ev.t1(2);
        ev.ecd1 = sum(filt1(inds))*pv.data.si;
        
        DoubleEvents(ind) = ev;
        
    end

    pv.setKeyboardCallback(@(x) ([]));
    
    % now start and loop through all the events, manually adjusting
    setSel(1);
    setEventView();
    
    while 1
        redrawEvents();
        
        try
            k = pv.waitKey();
        catch
            break
        end
        
        if k=='q'
            break
        elseif k=='d'
            % delete selected event
            if selectedEvent>0 && selectedEvent <= numel(DoubleEvents)
                DoubleEvents = DoubleEvents([1:selectedEvent-1, selectedEvent+1:numel(DoubleEvents)]);
                % leave selectedEvent var as is, for next event, and move on
                setEventView();
            end            
        elseif k=='0'
            % set bounds on first event
            if selectedEvent>0 && selectedEvent <= numel(DoubleEvents)
                DoubleEvents(selectedEvent).t0 = pv.getCursors();
                updateEvent(selectedEvent);
            end
        elseif k=='1'
            % set bounds on first event
            if selectedEvent>0 && selectedEvent <= numel(DoubleEvents)
                DoubleEvents(selectedEvent).t1 = pv.getCursors();
                updateEvent(selectedEvent);
            end
        elseif k=='m'
            % reset the mean/median of the event (which event?)
            ts = pv.getCursors();
            dd = pv.data.getByTime(ts);
            if selectedEvent>0 && selectedEvent <= numel(DoubleEvents)
                ev = DoubleEvents(selectedEvent);
                if ts(1) > ev.t0(2)
                    ev.m0 = mean(dd(:,6));
                elseif ts(2) < ev.t1(1)
                    ev.m1 = mean(dd(:,7));
                end
                DoubleEvents(selectedEvent) = ev;
                updateEvent(selectedEvent);
            end
        elseif k=='n'
            if selectedEvent < numel(DoubleEvents)
                selectedEvent = selectedEvent + 1;
                fprintf(1,'%d/%d\n',selectedEvent,numel(DoubleEvents));
            end
            setEventView();
        elseif k=='p'
            if selectedEvent > 1
                selectedEvent = selectedEvent - 1;
            end
            setEventView();
        elseif k=='s'
            selectedEvent = -1;
        elseif k=='v'
            t = mean(pv.getCursors());
            if selectedEvent>0 && selectedEvent <= numel(DoubleEvents)
                if abs(t-mean(DoubleEvents(selectedEvent).t0)) < ...
                        abs(t-mean(DoubleEvents(selectedEvent).t1))
                    pv.setView(DoubleEvents(selectedEvent).t1 + 2e-3*[-1 1]);
                    pv.setCursors(DoubleEvents(selectedEvent).t1);
                else
                    pv.setView(DoubleEvents(selectedEvent).t0 + 2e-3*[-1 1]);
                    pv.setCursors(DoubleEvents(selectedEvent).t0);
                end
            end
        elseif k=='o'
            plot_signals(pv);
        end
    end
end