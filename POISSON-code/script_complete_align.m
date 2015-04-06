% Align all the molecules completely

% %% initial
% 
% clear all
% pd = PoreData('/Users/Stephen/Documents/Stephen/Research/Data/MinION/Lambda-73/');
% bm = BioMap('lambda.bam');
% % reference sequence
% f = fastaread('./References/Lambda_NEB.fasta');
% refseq = f.Sequence;
% clear f;
% params = PoreParams.Load('CS_params.conf');
% 
% %% align (takes a long time)
% 
% inds = bm.getIndex(1,numel(refseq));
% inds = unique(inds);
% inds = inds(bm.Flag(inds) == 0);
% pdinds = pd.getIndex(bm.Header(inds));
% events = pd.getEvents(pdinds);
% %% (can iterate this)
% events = seedaligns(refseq,events,params);
% %save('/Users/Stephen/Documents/Stephen/Research/Analysis/MinION/20150323/aligned_events.mat','events','refseq')

%% load previously aligned data

clear all
load('/Users/Stephen/Documents/Stephen/Research/Analysis/MinION/20150323/aligned_events.mat')

%% coverage

centers = 1:numel(refseq);
total = [];
% for i = 1:numel(events)
%     total = [total; events(i).ref_align(events(i).ref_align>0)];
% end
coverage = arrayfun(@(x) hist(x.ref_align(x.ref_align>0),centers), events, 'UniformOutput', false);
coverage = cell2mat(coverage');
figure(1)
bar(centers,sum(coverage,1))
xlabel('Position along lambda')
ylabel('Number of molecules read')
title('Lambda coverage')
set(gca,'FontSize',24)

%% systematic errors of the model

refstates = seqtostates(refseq);

figure(2)
colormap parula
positions = (1:20) + 2161;
levels = 0*positions;
leverr = 0*positions;
levelsScaled = nan(numel(positions),numel(events));
for i = 1:numel(positions)
    
    evlevels = nan(numel(events),1);
    for j = 1:numel(events)
        if events(j).model.complement
            continue
        end
        lev = mean(events(j).mean(events(j).ref_align==positions(i)));
        
        scaling = max(events(j).model.level_mean) - min(events(j).model.level_mean);
        offset = min(events(j).model.level_mean);
        
        if ~isempty(lev)
            levelsScaled(i,j) = (lev-offset)/scaling; % levels themselves, scaled
            evlevels(j) = lev - events(j).model.level_mean(refstates(positions(i))); % discrepancy
        end
    end
    levels(i) = nanmean(evlevels);
    leverr(i) = nanstd(evlevels);
   
end
errorbar(positions,levels,leverr,'o')
ylabel('Measurement - model current (pA)')
xlabel('Position along lambda')
xlim([positions(1)-5 positions(end)+5])
title('Discrepancy between model and measurement')
set(gca,'FontSize',24)

% the levels themselves

figure(3)
errorbar(positions,nanmean(levelsScaled,2)',nanstd(levelsScaled,1,2)','ko')
hold on
plot(positions,levelsScaled,'o-')
ylabel('Scaled, measured current')
xlabel('Position along lambda')
xlim([positions(1)-5 positions(end)+5])
title('Levels')
set(gca,'FontSize',24)

%% specific k-mers

% segments of interest
clear segmentData n s
delta = 0; % start's offset from first base
segment = 'AAAACAAAA';
mer = 5;
for i = 1:numel(segment)-mer+1
    segmentData(i).sequence = segment(i:i+mer-1);
end
%segmentData(1).sequence = 'GGGAA';
% segmentData(1).sequence = 'CGGGCT';
% segmentData(2).sequence = 'AGGGCT';
% segmentData(3).sequence = 'GGGGCT';
% segmentData(4).sequence = 'TGGGCT';
% segmentData(5).sequence = 'AAGGAA';
% segmentData(6).sequence = 'ACGGAA';
% segmentData(7).sequence = 'CAGGAA';
%segmentData(8).sequence = 'GGGGG';

% find them in the reference sequence and the data
for i = 1:numel(segmentData)
    segmentData(i).start = strfind(refseq,segmentData(i).sequence)+delta;
    n(i) = numel(segmentData(i).start);
    s{i} = segmentData(i).sequence;
    segmentData(i).mean = [];
    
    positions = segmentData(i).start;
    levels = 0*positions;
    leverr = 0*positions;
    levelsScaled = nan(numel(positions),numel(events));
    for k = 1:numel(positions)
        
        evlevels = nan(numel(events),1);
        for j = 1:numel(events)
            if events(j).model.complement
                continue
            end
            lev = mean(events(j).mean(events(j).ref_align==positions(k)));
            
            scaling = max(events(j).model.level_mean) - min(events(j).model.level_mean);
            offset = min(events(j).model.level_mean);
            
            if ~isempty(lev)
                levelsScaled(k,j) = (lev-offset)/scaling; % levels themselves, scaled
                evlevels(j,k) = lev - events(j).model.level_mean(refstates(positions(k))); % discrepancy
            end
        end
        
    end
    %segmentData(i).mean = [segmentData(i).mean; reshape(levelsScaled,[],1)];
    segmentData(i).mean = levelsScaled;
    segmentData(i).discrepancies = evlevels';
    
end
% errorbar(positions,levels,leverr,'o')
% ylabel('Measurement - model current (pA)')
% xlabel('Position along lambda')
% xlim([positions(1)-5 positions(end)+5])
% title('Discrepancy between model and measurement')
% set(gca,'FontSize',24)

% plot the number of each
figure(4)
clf
bar(n)
set(gca,'XTick',1:numel(segmentData),'XTickLabel',s,'FontSize',18)
title('Number of of times each 5-mer occurs in lambda reference')
xlabel([num2str(mer) '-mer sequence, 5'' first'])
ylabel('Number of occurrences in reference')

% plot the normalized current levels
h = figure(5);
clf
colormap default
c = get(groot,'DefaultAxesColorOrder');
x3 = (1:0.001:numel(segmentData));
scaledCurrents = arrayfun(@(x) nanmean(reshape(x.mean,[],1)), segmentData)*scaling+offset;
scaledCurrentsE = arrayfun(@(x) nanstd(reshape(x.mean,[],1)*scaling+offset), segmentData);
errorbar(1:numel(segmentData),scaledCurrents,scaledCurrentsE,'o')
hold on
line([0 numel(segmentData)+1],[scaledCurrents(1) scaledCurrents(1)],'LineStyle','--','Color','k')
yy = spline(-1:numel(segmentData)+2,[scaledCurrents(1) scaledCurrents(1) scaledCurrents scaledCurrents(end) scaledCurrents(end)],x3);
plot(x3,yy,':','Color',c(1,:))
set(gca,'XTick',1:numel(segmentData),'XTickLabel',s,'FontSize',22)
title('MinION data from lambda DNA')
xlabel([num2str(mer) '-mer sequence, 5'' first'])
ylabel('Scaled current (pA)')
ylim([0 1]*scaling+offset)
set(gca,'LooseInset',[0 0 0 0]) % the all-important elimination of whitespace!
set(gca,'OuterPosition',[0 0 0.99 1]) % fit everything in there
%set(h,'Position',[100 500 1000 400]) % size the figure

% histogram
clear fullHist
div = 0.01;
x = (-0.1:div:1.2) * scaling + offset;
x2 = (-0.1:div/10:1.2) * scaling + offset;
h = figure(6);
clf
for i = 1:numel(segmentData)
    fullHist(i,:) = hist(reshape(segmentData(i).mean,[],1) * scaling + offset,x);
    fullHist(i,:) = fullHist(i,:)./(sum(fullHist(i,:))*div*scaling);
    %histogram(segmentData(i).mean,x,'FaceColor',c(i,:),'EdgeColor',c(i,:),'DisplayStyle','stairs','Normalization','pdf');
    %stairs(x,fullHist(i,:),'Color',c(mod(i,size(c,1)),:),'LineWidth',3)
    stairs(x,fullHist(i,:),'LineWidth',3)
    hold on
end
%title('Multi-modal levels: MinION data, lambda DNA')
xlabel('Scaled current (pA)')
ylabel('Probability density')
set(gca,'LooseInset',[0 0 0 0]) % the all-important elimination of whitespace!
set(gca,'OuterPosition',[0 0 0.99 1]) % fit everything in there
%set(h,'Position',[100 500 850 400]) % size the figure
set(gca,'FontSize',22)
xlim([-0.1 1.2]*scaling + offset)
l = legend(s);
set(l,'FontName','Courier New','FontWeight','Bold')
set(gca,'ColorOrderIndex',1);
for i = 1:numel(segmentData)
    plot(x2,normpdf(x2,nanmean(reshape(segmentData(i).mean,[],1)*scaling + offset), ...
        nanstd(reshape(segmentData(i).mean,[],1)*scaling + offset)),'LineWidth',1,'LineStyle','--')
end
