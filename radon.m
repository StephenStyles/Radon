clear;
close all;

t_s=15; % in seconds
total_time = 60 * 60 * 1/2; % in seconds
radon_level = 8000; % in pCi/L

t = 0:t_s:total_time;
ns1 = radon_level * 3.7e-2 * 0.3 * t_s; 

% For Rn-222
half_lives222=[3.825*24*60*60 3.05*60 26.8*60 19.9*60 164.3e-6];
alphas222    =[             1       1       0       0        1]; 
counts222 = zeros(length(t), length(half_lives222));
for i = 1:length(t)
    counts222(i,:) = decay_interval(t(i), t(i)+t_s, half_lives222).*alphas222;
end

sum_counts222 = sum(counts222,2);
exp_counts222 = decay_counts(ns1/sum_counts222(1),t_s,length(t),half_lives222,alphas222);
exp_counts222 = exp_counts222./ns1;
plot(t/60, sum_counts222./sum_counts222(1))
hold on;
plot(t/60, exp_counts222)

% For Rn-220
half_lives220 = [54.5 0.158 10.64*60*60 60.55*60];
alphas220     = [   1     1           0        1];
counts220 = zeros(length(t), length(half_lives220));
for i=1:length(t)
    counts220(i,:) = decay_interval(t(i), t(i)+t_s, half_lives220).*alphas220; 
end
sum_counts220 = sum(counts220,2);
exp_counts220 = decay_counts(ns1/sum_counts220(1),t_s,length(t),half_lives220,alphas220);
exp_counts220 = exp_counts220./ns1;
plot(t/60, sum_counts220./sum_counts220(1))
plot(t/60, exp_counts220)

exp_countsmix = decay_counts(0.5*ns1/sum_counts222(1),t_s,length(t),half_lives222,alphas222) + decay_counts(0.5*ns1/sum_counts220(1),t_s,length(t),half_lives220,alphas220);
exp_countsmix = exp_countsmix./ns1;

plot(t/60, (sum_counts220./sum_counts220(1) + sum_counts222./sum_counts222(1))/2)
plot(t/60, exp_countsmix)



