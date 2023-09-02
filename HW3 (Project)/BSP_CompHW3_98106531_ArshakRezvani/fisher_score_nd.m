function [score] = fisher_score_nd(features, labels)
    v = features;
    v1 = v(labels==0,:);
    v2 = v(labels==1,:);
    s1 = zeros(size(v,2));
    s2 = zeros(size(v,2));
    for i = 1:length(v1)
        s1 = s1 + (v1(i,:)-mean(v1,1))' * (v1(i,:)-mean(v1,1));
    end
    for i = 1:length(v2)
        s2 = s2 + (v2(i,:)-mean(v2,1))' * (v2(i,:)-mean(v2,1));
    end
    sw = (s1/length(v1)+s2/length(v2));
    sb = (mean(v1,1)-mean(v,1))'*(mean(v1,1)-mean(v,1)) ...
        +(mean(v2,1)-mean(v,1))'*(mean(v2,1)-mean(v,1));
    %score = trace(sb)/trace(sw);
    score = det(sb)/det(sw);
    %score = trace(sb .* sw.^-1);
end