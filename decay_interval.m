function proportions = decay_interval(t1,t2,half_lives)
    rates = log(2)./half_lives;
    proportions = zeros(size(half_lives));
    for j = 1:length(half_lives)
        for r = 1:j
            tmp = 1;
            for q=1:j
                if q ~= r
                    tmp = tmp * rates(q) / (rates(q)-rates(r));
                end
            end
            tmp = tmp * (exp(-rates(r)*t1) - exp(-rates(r)*t2));
            proportions(j) = proportions(j) + tmp;
        end
    end
end

