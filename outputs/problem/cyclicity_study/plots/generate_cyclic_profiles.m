%=====================
# PLOTTING FUNCTION
%=====================
% 
function generate_cyclic_profiles ( n, data_shmax, data_shvol, data_bghsz , t , z_max, plt_title, fig_name, visibility ) 

nr=1; #num of rows in plot
nc=3; #num of cols in plot

figure(n,'visible',visibility);
set(gcf, 'PaperSize', [15 5], 'PaperPosition', [0 0 15 5]);

subplot(nr, nc, 1);
plot (t, data_shmax, 'k','linewidth',4);
% add label to the x-axis
xlabel('time [Ma BP]', "fontsize", 25);
% add label to the y-axis
ylabel('max. sh [-]', "fontsize", 25);
ylim ([0 1]);
% add title for the plot
#title("...","fontsize", 25);
set(gca, 'FontSize', 15)
set(gca, "xdir", "reverse")

subplot(nr, nc, 2);
plot (t, data_shvol, 'k','linewidth',4);
% add label to the x-axis
xlabel('time [Ma BP]', "fontsize", 25);
% add label to the y-axis
ylabel('sh sp. vol. [m^3/m^2]', "fontsize", 25);
ymax = max(data_shvol) + 0.5*(max(data_shvol)-min(data_shvol));
ylim([0 max(1,ymax)]);
% add title for the plot
title(plt_title,"fontsize", 25);
set(gca, 'FontSize', 15)
set(gca, "xdir", "reverse")

subplot(nr, nc, 3);
plot (t, data_bghsz, 'k','linewidth',4);
% add label to the x-axis
xlabel('time [Ma BP]', "fontsize", 25);
% add label to the y-axis
ylabel('bGHSZ [mbsf]', "fontsize", 25);
ylim([0 z_max]);
% add title for the plot
#title("...","fontsize", 25);
set(gca, 'FontSize', 15)
set(gca, "xdir", "reverse")
set(gca, "ydir", "reverse")

saveas (n, fig_name);
  
endfunction