%%

clear all
iterations = 100;
mean_length = 100e-3;
filt = 2000;
do_filter = true;
FPS = 1e-5;
levels = 100;
sample_freq = 1e4;

num_levs = zeros(1,iterations);

for i = 1:iterations
    seq = int2nt(randi(4,1,levels+4));
    [pA, pAstdmean, pAstd] = oxford_simulator(seq,1,0);
    timings = exprnd(mean_length,size(pA));
    pts = round(timings*sample_freq);
    current = [];
    for j = 1:numel(pts)
        current = [current; randn(pts(j),1)*pAstd(j) + pA(j)];
    end
    time = linspace(0,sum(timings),numel(current))';
    
    if do_filter
        d = filt_lpb([time, current], 4, filt);
        time = d(:,1);
        current = d(:,2);
    end
    
    lev = karplus_levels([time,current], 1/mean_length/10000, FPS, filt);
    num_levs(i) = numel(lev);
    disp(i)
    
end

%%

figure(1)
stairs(1:200,hist(num_levs,1:200))
hold on

%%

figure(2)
clf
for i = 1:numel(lev)
    inds = [find(time>=lev{i}.start_time,1,'first'), find(time>=lev{i}.end_time,1,'first')];
    plot(time(inds(1):inds(2)),current(inds(1):inds(2)))
    hold on
end
% plot(time,current)
% hold on
% line(cell2mat(cellfun(@(x) [x.start_time, x.end_time], lev, 'uniformoutput', false))',(cellfun(@(x) x.current_mean, lev)*ones(1,2))','linewidth',2,'color','r')

