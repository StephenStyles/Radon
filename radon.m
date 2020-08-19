clear;
close all;

t = 0:15:3000;

% For Rn-222
half_lives=[3.825*24*60*60 3.05*60 ];%26.8*60+19.7*60+1.5e-4];
counts = zeros(length(t), length(half_lives));
for i = 1:length(t)
    counts(i,:) = decay_interval(t(i), t(i)+15, half_lives);
    fprintf("t: %d-%d : %d counts\n", t(i), t(i)+15, sum(counts(i,:)));
end

sum_counts = sum(counts,2);
plot(t, sum_counts./sum_counts(1))
hold on;

% For Rn-220
half_lives = [54.5 0.158];
counts = zeros(length(t), length(half_lives));
for i=1:length(t)
    counts(i,:) = decay_interval(t(i), t(i)+15, half_lives);
end
sum_counts = sum(counts,2);
plot(t, sum_counts./sum_counts(1))


