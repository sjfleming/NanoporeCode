% Align molecules to a template

%% initial

clear all
%pd = PoreData('/Users/Stephen/Documents/Stephen/Research/Data/MinION/Lambda-73/');
bm = BioMap('lambda.bam');
% reference sequence
f = fastaread('./References/puc19.fasta');
refseq = f.Sequence;
refseq = [refseq(397:end) refseq(1:396)]; % the EcoRI cut
clear f;
params = PoreParams.Load('CS_params.conf');

%% align (takes a long time)

inds = bm.getIndex(1,numel(refseq));
inds = unique(inds);
inds = inds(bm.Flag(inds) == 0);
pdinds = pd.getIndex(bm.Header(inds));
events = pd.getEvents(pdinds);
%% (can iterate this)
events = seedaligns(refseq,events,params);
%save('/Users/Stephen/Documents/Stephen/Research/Analysis/MinION/20150323/aligned_events.mat','events','refseq')

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

