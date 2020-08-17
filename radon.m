clear;
close all;

t = 0:15:165;
% For Rn-220
half_lives=[54.5 0.158];
counts = zeros(length(t), length(half_lives));
for i = 1:length(t)
    counts(i,:) = decay_interval(t(i), t(i)+15, half_lives);
    fprintf("t: %d-%d : %d counts", t(i), t(i)+15, sum(counts(i,:)));
end