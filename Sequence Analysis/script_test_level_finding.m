%%

clear all
iterations = 100;
mean_length = 100e-3;
filt = 500;
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
    
    lev = karplus_levels([time,current], 1/mean_length, FPS, filt);
    num_levs(i) = numel(lev);
    disp(i)
    
end
