function pA = oxford_simulator(seq,model,plotFlag)

% Takes in a sequence in the form of 'ATCGTTCAAAGC...' and plots what
% Oxford's simulated squiggles would look like.

% generate "states"
states = get_states(nt2int(seq), 5); % k=5 for k-mer

% get the current from Oxford's models
load('models.mat')
pA = model_data{model}.level_mean(states);

% plot it
if plotFlag==1
    figure(2)
    clf(2)
    plot(1:numel(pA),pA,'o--');
    ylim([min(model_data{model}.level_mean) max(model_data{model}.level_mean)])
    title('Oxford Nanopore predicted "squiggle" from molecule','FontSize',20)
    xlabel('Base Number','FontSize',20)
    ylabel('Current (pA)','FontSize',20)
    set(gca,'FontSize',22)
end

end