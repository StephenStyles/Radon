close all
% You need to import the data here, sorry!
% You also need to set up X and Y correctly.
surf(X, Y, rn220std./150);
xlabel("Sampling Duration","FontSize",20) % step size variable
ylabel("Sampling Period","FontSize",20)
zlabel("(SD of estimated amount) / (actual amount)","FontSize",20)
title("SD in estimate using 150 Rn-220 particles","FontSize",20)
ax=gca;
ax.FontSize=20;
[caz,cel] = view;
view(-caz,cel)

figure()
surf(X,Y,rn222std./750000)
xlabel("Sampling Duration (minutes)","FontSize",20)
ylabel("Sampling Period (seconds)","FontSize",20)
zlabel("(SD of estimated amount) / (actual amount)","FontSize",20)
title("SD in estimate using 750,000 Rn-222 particles","FontSize",20)
ax=gca;
ax.FontSize=20;
[caz,cel] = view;
view(-caz,cel)