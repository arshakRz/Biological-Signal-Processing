function [a, b] = myARMA(sig, p, q)
    a = myPAR(sig, p); 
    ar_filtered_signal = filter(-1 * a(2:end), 1, sig);
    b =myMA(ar_filtered_signal, q);
end