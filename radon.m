clear;
close all;

t_s=15; % in seconds
total_time = 60 * 15; % in seconds
radon_level = 8000; % in pCi/L
num_runs = 10;

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
%exp_counts222 = decay_counts(ns1/sum_counts222(1),t_s,length(t),half_lives222,alphas222);
%exp_counts222 = exp_counts222./ns1;
plot(t/60, sum_counts222./sum_counts222(1),'LineWidth',2,'DisplayName','Rn-222 Expected')
hold on;
%plot(t/60, exp_counts222,'LineWidth',2,'DisplayName','Rn-222 Simulated')

% For Rn-220
half_lives220 = [54.5 0.158 10.64*60*60 60.55*60];
alphas220     = [   1     1           0        1];
counts220 = zeros(length(t), length(half_lives220));
for i=1:length(t)
    counts220(i,:) = decay_interval(t(i), t(i)+t_s, half_lives220).*alphas220; 
end
sum_counts220 = sum(counts220,2);
%exp_counts220 = decay_counts(ns1/sum_counts220(1),t_s,length(t),half_lives220,alphas220);
%exp_counts220 = exp_counts220./ns1;
plot(t/60, sum_counts220./sum_counts220(1),'LineWidth',2,'DisplayName','Rn-220 Expected')
%plot(t/60, exp_counts220,'LineWidth',2,'DisplayName','Rn-220 Simulated')

%exp_countsmix = decay_counts(0.5*ns1/sum_counts222(1),t_s,length(t),half_lives222,alphas222) + decay_counts(0.5*ns1/sum_counts220(1),t_s,length(t),half_lives220,alphas220);
%exp_countsmix = exp_countsmix./ns1;
for p =0.1:0.1:0.9
    plot(t/60, p*sum_counts220./sum_counts220(1) + (1-p)*sum_counts222./sum_counts222(1),'LineWidth',2,'HandleVisibility','off')
end
%plot(t/60, exp_countsmix,'LineWidth',2,'DisplayName','Equal Activity Simulated')
legend({},'FontSize',20)
xlabel("Time (minutes)",'FontSize',20)
ylabel("Alpha decay counts per sample period",'FontSize',20)
title("Normalized expected count",'FontSize',20)
%title(sprintf("Activity: %d pCi/L",radon_level),'FontSize',20)
% exp_countslpf = movmean(exp_countsmix,11);
% plot(t/60, exp_countslpf)
% 
% actual = [0.5*ns1/sum_counts222(1);0.5*ns1/sum_counts220(1)];
% input = [sum_counts222 sum_counts220];
% lr_mat = (input'*input)\input';
% errors_unf = zeros(2,num_runs);
% errors_lpf = zeros(2,num_runs);
% for i=1:num_runs
%     exp_countsmix = decay_counts(0.5*ns1/sum_counts222(1),t_s,length(t),half_lives222,alphas222) + decay_counts(0.5*ns1/sum_counts220(1),t_s,length(t),half_lives220,alphas220);
%     exp_countslpf = movmean(exp_countsmix, 11);
%     guess_unf = lr_mat*exp_countsmix;
%     guess_lpf = lr_mat*exp_countslpf;
%     errors_unf(:,i) = guess_unf-actual;
%     errors_lpf(:,i) = guess_lpf-actual;
% end
% 
% fprintf("Unfiltered %% error in Radon-222 measurement: mean = %.2f%%, std = %.2f%%\n",mean(100*errors_unf(1,:)/actual(1)),std(100*errors_unf(1,:)/actual(1)))
% fprintf("Low pass filtered %% error in Radon-222 measurement: mean = %.2f%%, std = %.2f%%\n",mean(100*errors_lpf(1,:)/actual(1)),std(100*errors_lpf(1,:)/actual(1)))
% fprintf("Unfiltered %% error in Radon-220 measurement: mean = %.2f%%, std = %.2f%%\n",mean(100*errors_unf(2,:)/actual(2)),std(100*errors_unf(2,:)/actual(2)))
% fprintf("Low pass filtered %% error in Radon-220 measurement: mean = %.2f%%, std = %.2f%%\n",mean(100*errors_lpf(2,:)/actual(2)),std(100*errors_lpf(2,:)/actual(2)))
% 
% 
% 
