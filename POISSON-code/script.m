%% Stephen Fleming, March 9, 2015
% Script to analyze current level data from a MinION run using POISSON.

%% Open data and create relevant objects
clear all
pd = PoreData('/Users/Stephen/Documents/Stephen/Research/Data/MinION/Lambda-73/');
bm = BioMap('lambda.bam');
% reference sequence
f = fastaread('./References/Lambda_NEB.fasta');
seq = f.Sequence;
clear f;

%% Determine regions of interest

% segments of interest
segmentData(1).sequence = 'TTTTT';
segmentData(2).sequence = 'TTTTA';
segmentData(3).sequence = 'TTTAT';
segmentData(4).sequence = 'TTATT';
segmentData(5).sequence = 'TATTT';
segmentData(6).sequence = 'ATTTT';
segmentData(7).sequence = 'TTTTT';
segmentData(8).sequence = 'TTTTC';
segmentData(9).sequence = 'TTTCT';
segmentData(10).sequence = 'TTCTT';
segmentData(11).sequence = 'TCTTT';
segmentData(12).sequence = 'CTTTT';
segmentData(13).sequence = 'TTTTT';
segmentData(14).sequence = 'TTTTG';
segmentData(15).sequence = 'TTTGT';
segmentData(16).sequence = 'TTGTT';
segmentData(17).sequence = 'TGTTT';
segmentData(18).sequence = 'GTTTT';

% find them in the reference sequence
for i = 1:numel(segmentData)
    segmentData(i).start = strfind(seq,segmentData(i).sequence);
    n(i) = numel(segmentData(i).start);
    s{i} = segmentData(i).sequence;
end
clear i candidates

% plot the number of each
figure(1)
clf
bar(n)
set(gca,'XTick',1:numel(segmentData),'XTickLabel',s,'FontSize',18)
title('Number of of times each 5-mer occurs in lambda reference')
xlabel('5-mer sequence')
ylabel('Number of occurrences in reference')

params = PoreParams.Load('CS_params.conf');

%% Align event levels to nucleotides in reference sequence
a = 1;
for i = 1:numel(segmentData) % go through each k-mer
    for j = 1:numel(segmentData(i).start) % go through each location in lambda (all molecules in each loop)
        
        if a==1
            tic;
        end
        % pull out sequences that align there (plus minus 200)
        i0 = segmentData(i).start(j)-200;
        i1 = segmentData(i).start(j)+200;
        i0 = max(i0,1);
        i1 = min(i1,numel(seq));
        inds = bm.getIndex(i0,i1);
        inds = inds(bm.Flag(inds) == 0); % template only
        events = arrayfun(@(x) pd.getEvent(x,'t'),pd.getIndex(bm.Header(inds)));
        curseq = seq(i0:i1);
        [~,events] = align_likes(curseq,events,params);
        kmerind = segmentData(i).start(j)-i0+1;
        segmentData(i).mean(j) = nanmean(arrayfun(@(x) mean(x.mean(x.ref_align==kmerind)) - x.model.level_mean(end), events)); % subtract lowest (TTTTT) level in each bc of diff scaling
        segmentData(i).individualMoleculeMeans{j} = arrayfun(@(x) mean(x.mean(x.ref_align==kmerind)) - x.model.level_mean(end), events); % subtract lowest (TTTTT) level in each bc of diff scaling
        % times a scaling maybe max - min...
        if a==1
            t = toc;
            display(['This will take ' num2str(round(t*sum(n)/60)) ' minutes'])
        end
        display([num2str(round(a/sum(n)*100,2)) '% complete'])
        a = a+1;
        
    end
end

%% Plot the result

figure(2)
clf
polyT = mean(arrayfun(@(x) x.model.level_mean(end), events)); % add mean TTTTT level back
errorbar(arrayfun(@(x) mean(x.mean), segmentData) + polyT, arrayfun(@(x) std(x.mean), segmentData),'o','MarkerSize',8)
hold on
line([0 numel(segmentData)+1], [mean(segmentData(1).mean) mean(segmentData(1).mean)] + polyT,'LineStyle','--','Color',[0.6 0.6 0.6])
set(gca,'XTick',1:numel(segmentData),'XTickLabel',s,'FontSize',18)
ylim([50 70])
title('Scaled current levels corresponding to 5-mers','FontSize',20)
xlabel('5-mer sequence','FontSize',20)
ylabel('Current (pA)','FontSize',20)

figure(3)
clf
%hold on
polyT = mean(arrayfun(@(x) x.model.level_mean(end), events)); % add mean TTTTT level back
%for i = 1:numel(segmentData)
%    plot(i*ones(1,numel(segmentData(i).mean)),segmentData(i).mean + polyT,'o','Color',[0.9 0.9 0.9])
%end
%hold on
errorbar(arrayfun(@(x) mean(x.mean), segmentData(1:6)) + polyT, arrayfun(@(x) std(x.mean), segmentData(1:6)),'o:','MarkerSize',8,'Color',[0,0.4470,0.7410])
hold on
errorbar(arrayfun(@(x) mean(x.mean), segmentData(7:12)) + polyT, arrayfun(@(x) std(x.mean), segmentData(7:12)),'*:','MarkerSize',8,'Color',[0.8500,0.3250,0.0980])
errorbar(arrayfun(@(x) mean(x.mean), segmentData(13:18)) + polyT, arrayfun(@(x) std(x.mean), segmentData(13:18)),'s:','MarkerSize',8,'Color',[0.4660,0.6740,0.1880])
line([0 7], [mean(segmentData(1).mean) mean(segmentData(1).mean)] + polyT,'LineStyle','--','Color',[0.6,0.6,0.6])
xlim([0 7])
ylim([50 70])
box on
set(gca,'XTick',1:6,'XTickLabel',{'TTTTT','TTTTX','TTTXT','TTXTT','TXTTT','XTTTT'},'FontSize',18)
legend('X = A','X = C','X = G')
title('Scaled current levels corresponding to 5-mers','FontSize',20)
xlabel('5-mer sequence','FontSize',20)
ylabel('Current (pA)','FontSize',20)

%% Make the plot with the single A substitution

load('/Users/Stephen/Documents/Stephen/Research/Analysis/MinION/20150309/lambda73_polyT_singleA.mat')
for i = 1:numel(segmentData)
    s{i} = segmentData(i).sequence;
end
figure(4)
clf
polyT = mean(arrayfun(@(x) x.model.level_mean(end), events)); % add mean TTTTT level back
x = 1:6;
y = arrayfun(@(x) mean(x.mean), segmentData) + polyT;
errorbar(x, y, arrayfun(@(x) std(x.mean), segmentData),'o','MarkerSize',8,'Color',[0 0.4470 0.7410])
hold on
xx = 0:0.1:7;
yy = spline([0 x 7],[polyT y polyT],xx);
plot(xx,yy,':')
line([0 numel(segmentData)+1], [mean(segmentData(1).mean) mean(segmentData(1).mean)] + polyT,'LineStyle','--','Color',[0.6 0.6 0.6])
set(gca,'XTick',1:numel(segmentData),'XTickLabel',s,'FontSize',18)
ylim([50 70])
xlim([0.3 6.7])
title('Scaled current levels corresponding to 5-mers','FontSize',20)
xlabel('5-mer sequence','FontSize',26)
ylabel('Current (pA)','FontSize',26)
%annotation('doublearrow',([1 4]-0.3)/(6.7-0.3),([mean([y(2) y(5)]) mean([y(2) y(5)])]-polyT)/(70-50))
annotation('doublearrow',[0.335 0.705],[0.4 0.4])
annotation('textbox',[0.42 0.47 0 0],'String','4-positions','FontSize',20)
