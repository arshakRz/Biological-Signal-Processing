function [a, b] = fit_param(sig, model_type, p, q)
    if model_type == "AR"
        a = myAR(sig, p);
        b = 1;
    elseif model_type == "MA"
        b = myMA(sig, q);
        a = 1;
    else 
        [a, b] = myARMA(sig, p, q);
    end
end