%=====================
# PLOTTING FUNCTION
%=====================
% 
function generate_timeseries_plot ( n, data , z , t , plt_title, z_label, t_label, c_label, fig_name ) 

figure(n);
set(gcf, 'PaperSize', [15 10], 'PaperPosition', [0 0 15 10]);
imagesc(t, z, data);

colormap (jet(32));
c = colorbar;% get current colormap
set(c,'fontsize',25);
title(c, c_label, 'FontSize', 25);
xlabel( t_label );
ylabel( z_label );
title( plt_title, "fontsize", 25);
set(gca, 'FontSize', 25)
set(gca, "xdir", "reverse")

saveas (n, fig_name);
  
endfunction
