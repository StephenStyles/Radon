function [countlist] = decay_counts(n,sample_time, n_samples, half_lives, count)
    if size(count) ~= size(half_lives)
        count = ones(size(half_lives))
    end
    countlist = zeros(n_samples, 1);
    for i = 1:n
        templist = exprnd(half_lives/log(2));
        templist = floor(cumsum(templist)/sample_time);
        for j=1:length(templist)
            t=templist(j)+1;
            if t <= n_samples
                countlist(t) = countlist(t) + count(j);
            end
        end
    end
end