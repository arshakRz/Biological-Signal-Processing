function a = myAR(sig, p)
    autocorr_sig = xcorr(sig, p, 'biased');
    R = toeplitz(autocorr_sig(p+1:end-1));
    r = autocorr_sig(p+2:end);
    a = -R\r;
end