function [outState] = decay_state(t, half_lives)
    state = zeros(size(half_lives));
    rates = log(2)./half_lives;
    for j=1:length(state)
        for r=1:j
            tmp = 1;
            for q=1:j
                if q ~= r
                   tmp = tmp * rates(q)/(rates(q)-rates(r));
                end
            end
            tmp = tmp * rates(r);
            tmp = tmp * exp(-rates(r)*t);
            state(j) = state(j) + tmp;
        end
        state(j) = state(j) / rates(j);
    end
    outState = state;
end



