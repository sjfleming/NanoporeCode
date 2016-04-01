function L = lorentzianpdf(x, mean, fwhm)
% LORENTZIANPDF(x, mean, fwhm) returns the Lorentzian function pdf
% http://mathworld.wolfram.com/LorentzianFunction.html
    
    L = (1/pi) * (fwhm/2) ./ ((x-mean).^2 + (fwhm/2).^2);
    
end