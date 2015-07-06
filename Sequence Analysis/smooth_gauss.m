function v = smooth_gauss(v, sig)
% Smooths v with a Gaussian kernel

    krn = normpdf(linspace(-3*sig,3*sig,ceil(1+6*sig)),0,sig);
    krn = krn / sum(krn);
    v = conv(v,krn,'same');

end