clear;
close all;

t = 0:15:3000;

% For Rn-222
half_lives=[3.825*24*60*60 3.05*60 26.8*60 19.9*60 164.3e-6];
counts    =[             1       1       0       0        1]; 
counts222 = zeros(length(t), length(half_lives));
for i = 1:length(t)
    counts222(i,:) = decay_interval(t(i), t(i)+15, half_lives).*counts;
end

sum_counts222 = sum(counts222,2);
plot(t, sum_counts222./sum_counts222(1))
hold on;

% For Rn-220
half_lives = [54.5 0.158 10.64*60*60 60.55*60];
counts     = [   1     1           0        1];
counts220 = zeros(length(t), length(half_lives));
for i=1:length(t)
    counts220(i,:) = decay_interval(t(i), t(i)+15, half_lives).*counts; 
end
sum_counts220 = sum(counts220,2);
plot(t, sum_counts220./sum_counts220(1))

plot(t, (sum_counts220./sum_counts220(1) + sum_counts222./sum_counts222(1))/2)




