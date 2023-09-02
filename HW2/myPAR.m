function a = myPAR(sig, p)
    N = length(sig);
    R = zeros(p+1, 1);
    for k = 0:p
        R(k+1) = sum(sig(1:N-k) .* sig(k+1:N)) / N;
    end    
    R_mat = toeplitz(R(1:end-1));
    a = -R_mat \ R(2:end);
end