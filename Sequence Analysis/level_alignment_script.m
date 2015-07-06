%% Script to align level data to a known reference
% Stephen Fleming
% 4/8/15

%% Align level data to reference

% load data file
clear all
file = '/Users/Stephen/Documents/Stephen/Research/Analysis/Biopore/20150311/0023_discreteData_1843.mat';
currentCutoffs = [40 70]; % cutoffs in pA for "real" levels
load(file)

% get reference sequence
f = fastaread('/Users/Stephen/GitHub/NanoporeCode/POISSON-code/References/puc19.fasta'); % pUC19
refseq = f.Sequence;
refseq = [refseq(397:end) refseq(1:396)]; % the EcoRI cut, linearized
clear f;

% f = fastaread('/Users/Stephen/GitHub/NanoporeCode/POISSON-code/References/Lambda_NEB.fasta'); % lambda
% refseq = f.Sequence;
% refseq = refseq(22426:26104); % the EcoRI, BglII double digest 3.6kb fragment
% clear f;

SK23 = 'GGTTGTTTCTGTTGGTGCTGATATTGCGGCGTCTGCTTGGGTGTTTAACCT'; % SK23 leader sequence
template = [SK23 refseq]; % add on the SK23 leader
complement = [SK23 seqrcomplement(refseq)]; % reverse completment of sequence

% set alignment parameters [sd = 2, si = -0.1]
sd = 1.75;
si = -0.3;

% plot data squiggle
figure(1)
clf(1)
plot(discreteData.levels,'o-')
hold on
line([1 numel(discreteData.levels)],[currentCutoffs(1) currentCutoffs(1)],'Color','r','LineStyle','--')
line([1 numel(discreteData.levels)],[currentCutoffs(2) currentCutoffs(2)],'Color','r','LineStyle','--')
ylabel('Level current (pA)')
xlabel('Level number')
set(gca,'FontSize',22)
figure(5)
clf(5)

for j = 1:4
    
    switch j
        case 1
            thetitle = 'template, model #1';
            seq = template;
            m = 1;
        case 2
            thetitle = 'template, model #2';
            seq = template;
            m = 2;
        case 3
            thetitle = 'complement, model #1';
            seq = complement;
            m = 1;
        case 4
            thetitle = 'complement, model #2';
            seq = complement;
            m = 2;
    end
    
    display([thetitle ' '])
    
    % create Oxford model for reference sequence
    refCurrentRaw = oxford_simulator(seq,m,0);
    refCurrent = normlevels(refCurrentRaw,discreteData.levels(and(discreteData.levels>currentCutoffs(1),discreteData.levels<currentCutoffs(2))));
    
%     figure(1)
%     clf
%     [f,x] = ecdf(discreteData.levels);
%     plot(x,f)
%     hold on
%     [f2,x2] = ecdf(refCurrent);
%     plot(x2,f2)
%     legend('Data','Model 1')
%     xlabel('Level Current (pA)')
%     ylabel('Cumulative Distribution Function')
%     set(gca,'FontSize',22)
    
    % example data
    %discreteData.levels = refCurrent([1000:1100 1102:1110 1113:1200]) + randn(198,1)/0.2;
    
    % align data to reference
    measuredLevels = discreteData.levels(and(discreteData.levels>currentCutoffs(1),discreteData.levels<currentCutoffs(2)));
    dpath = align_local(refCurrent, measuredLevels, sd, si);
    dpath = dpath - 1;
    newref = nan*dpath;
    newref(:,1) = refCurrent(dpath(:,1));
    newref(:,2) = measuredLevels(dpath(:,2));
    
    % plot the alignment
    figure(5)
    subplot(4,1,j)
    h = gcf;
    plot(dpath(:,1),newref(:,1),'o--','Color',[0.7,0.7,0.7],'LineWidth',1) % reference
    hold on
    plot(dpath(:,1),newref(:,2),'o-','LineWidth',2) % don't re-plot reference
    h_legend = legend('Model prediction','Data');
    set(h_legend,'FontSize',20);
    legend('boxoff')
    xlim([dpath(1,1)-1 dpath(end,1)+1])
    xlabel('Level number in reference','FontSize',20)
    ylabel('Current (pA)','FontSize',20)
    title(['E5/pUC19, ' thetitle],'FontSize',20)
    set(gca,'FontSize',20)
    
    display(' ')
    
end

annotation('textbox', [0.5 0.92 0 0], 'String', [name ' sd=' num2str(sd) ' si=' num2str(si)], 'FontSize', 20);
