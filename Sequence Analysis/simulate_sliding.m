function out = simulate_sliding(seq, voltage, temp, varargin)
% SIMULATE_SLIDING(sequence, voltage, temperature, optional params)
% generates simulated sliding data for a nanopore experiment with a 
% helicase sliding along DNA
% optional parameters:
% 'duplex' can be set to true or false.  true indicates the helicase is up
% against a dsDNA duplex.  flase indicates sliding along ssDNA.
% Stephen Fleming, 2/20/16

    % get estimate of current levels for sequence from Oxford's predictions
    levs = get_model_levels_oxford(seq, 2.3*voltage*linspace(0.15,0.35,20), 2.3*voltage, voltage);
    
    % generate a sequence of states based on transition probabilities
    
    % transition_prob = [forward prob, back prob]
    % for now, this is just linearly proportional to voltage...
    transition_prob = [8*(voltage-100)/(180-100)+1, 1]; % equal at 100mV, almost all forward at 180mV
    if transition_prob(1)<=0
        transition_prob = transition_prob(1)*-1 + transition_prob + 1;
    end
    transition_prob = transition_prob / sum(transition_prob);
    
    state = 1;
    t = 0;
    sequence = state;
    while state <= numel(levs)
        
        % at each time step, there is a Botzmann factor probability of
        % moving at all
        do_step = try_take_step(temp, rand(1));
        
        if do_step
            if rand(1) < transition_prob(1)
                state = state+1;
            else
                state = max(1,state-1);
            end
        end
        
        % add current state to sequence of states, and advance time
        sequence(end+1) = state;
        t = t+1;
        if mod(t,10)==0
            fprintf('.')
        end
        
    end
    
    fprintf('\n')
    
    out.t = 0:1:t;
    out.sequence = sequence;
    out.current_levels = levs;
    
    if nargout == 0
        figure()
        
        subplot(2,1,1)
        plot(0:1:t,sequence,'o-')
        xlim([0, max(t)])
        ylim([0, numel(levs)+1])
        set(gca,'fontsize',18)
        ylabel('state')
        xlabel('time step')
        title(['Voltage = ' num2str(voltage) ' mV, Temp = ' num2str(temp) ' C'])
        
        subplot(2,1,2)
        line([0:1:t-1;1:1:t]*1e-1,repmat(levs(sequence(1:end-1)),[2,1]),'LineWidth',3)
        xlim([0, max(t)*1e-1])
        ylim([min(levs)-10, max(levs)+10])
        set(gca,'fontsize',18)
        ylabel('current (pA)')
        xlabel('time (s)')
    end
    
    function bool = try_take_step(temp, seed)
        
        energy_barrier = 37;
        energy = temp;
        
        cooperativity = 8;
        factor = 100; % if timestep is 100ms and we want an attempt every 1ms
        prob = min(1, factor * (expcdf(energy, energy_barrier) .^ cooperativity));
        
        if seed <= prob
            bool = true;
        else
            bool = false;
        end
    end

end