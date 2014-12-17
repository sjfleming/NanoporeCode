function pA = oxford_simulator(seq,plotFlag)

% Takes in a sequence in the form of 'ATCGTTCAAAGC...' and plots what
% Oxford's simulated squiggles would look like.

% generate "states"
states = get_states(nt2int(seq), 5); % k=5 for k-mer

% get the current from Oxford's models
load('models.mat')
pA = model_data{1}.level_mean(states);
pA2 = model_data{3}.level_mean(states);
pA3 = model_data{5}.level_mean(states);
pA4 = model_data{7}.level_mean(states);

% plot it
if plotFlag==1
    figure(2)
    clf(2)
    axes('FontSize',20)
    plot(1:numel(states),pA,'o');
%     hold on
%     plot(1:numel(states),pA2,'or');
%     plot(1:numel(states),pA3,'og');
%     plot(1:numel(states),pA4,'ok');
    ylim([20 70])
    title('Oxford Nanopore predicted "squiggle" from molecule','FontSize',20)
    xlabel('Base Number','FontSize',20)
    ylabel('Current (pA)','FontSize',20)
end

end