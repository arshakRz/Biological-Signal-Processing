%{
function b = myMA(sig, q)
    N = length(sig);
    b = zeros(q, 1);
    gamma = zeros(q+1, 1);
    gamma(1) = sum(sig.^2) / N;
    for k = 1:q
        sum_term = 0;
        for j = 1:k-1
            sum_term = sum_term + b(j) * sig(k-j);
        end
        error_term = sig(k+1:N) - sum_term;
        b(k) = sum(error_term .* sig(k+1:N)) / sum(error_term.^2);
        b(1:k-1) = b(1:k-1) - b(k) * flipud(b(1:k-1));
        gamma(k+1) = (1 - b(k)^2) * gamma(k);
    end
end
%}
function ma_coeffs = myMA(signal, ma_order)
    N = length(signal);
    ma_coeffs = zeros(ma_order, 1);
    E = zeros(ma_order+1, 1);
    E(1) = signal' * signal / N;
    
    for k = 1:ma_order
        numerator = signal(k+1:N)' * signal(1:N-k);
        denominator = signal(1:N-k)' * signal(1:N-k);
        ma_coeffs(k) = numerator / denominator;
        
        E(k+1) = (1 - ma_coeffs(k)^2) * E(k);
        
        for j = 1:k-1
            ma_coeffs(j) = ma_coeffs(j) - ma_coeffs(k) * ma_coeffs(k-j);
        end
    end
end