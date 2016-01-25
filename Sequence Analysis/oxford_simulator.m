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
    plot(1:numel(pA),pA,'o-');
    ylim([min(model_data{model}.level_mean) max(model_data{model}.level_mean)])
    title(seq,'FontSize',16,'FontName','Courier')
    xlabel('Base Number')
    xlim([-1 numel(seq)-2])
    ylabel('Current (pA)')
    set(gca,'FontSize',16)
    set(gca,'LooseInset',[0 0 0 0]) % the all-important elimination of whitespace!
    set(gca,'OuterPosition',[0.01 0.01 0.98 0.98]) % fit everything in there
    set(gcf,'Position',[-1034 431 1021 278])
end

end