%%% Coastline Optimization 16.346
%%% Inclination data plot
%%% Andrew Adams & Tesla Wells


%% collect data
incangles = [0, 10, 20, 28.5, 30, 40, 50, 51.6, 60, 70, 80, 90];
percents = [2.70, 2.8866, 2.2921, 2.2253, 2.0946, 2.458, 2.7769, 3.0445, 3.1895, 4.39, 5.28, 4.3495];

plot(incangles, percents,'k-*')
set(gca,'fontsize',18)
axis([0 90 0 7])
xlabel('Inclination Angle')
ylabel('Coastline Time (%)')
title('Coverage Time (%) by inclination')
    