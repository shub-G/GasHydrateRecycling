## FOR MOVIE:
## ffmpeg -r 5 -i conc_%04d.png -b:v 1000k test.mp4

graphics_toolkit gnuplot

%=====================
# USER INPUTS
%=====================
K="1e-16";
kr="1e-18";
vs=0.05;
data_name   =strcat("../test1_K=",K,"_kr=",kr,"_vs=",num2str(vs,'%.4f'));
output_name =strcat("movie/test1_K=",K,"_kr=",kr,"_vs=",num2str(vs,'%.4f'));
t_END       = 27.0e6;   #annum #total calculation time 
t_recording = 0.0e6;    #annum #time when fine output is recorded
t_POI       = 0.0e6;    #annum #time when period of interest begins
t_offset    = 27.0e6;   #annum #offsets time axis s.t. t=0 corresponsts to t_offset
tOP_coarse  = 500000.0; #annum #coarse output-interval (before t_recording)
tOP_fine    = 2500.0;   #annum #fine output-interval (after t_recording)
Z           = 600.0;    #cm    #sediment column depth
Z_bottom    = 600.0;    #cm    #sediment column depth for plotting

%=====================
% READ DATA
%=====================
d_sh=strcat(char(data_name), "_sh_profile", ".txt");
data_sh=dlmread(d_sh, ' ', 0, 0);
                
[val,cz]=min(abs(data_sh(1,:)-Z_bottom));
z = data_sh(1,cz:end);
z_new = z;
[valt,rt]=min(abs(data_sh(:,1)-t_POI)); #row number where fine output starts
t = data_sh(rt:end,1);
t_new = t - t_offset ;
t_new *= 1.e-6;

d_sg=strcat(char(data_name), "_sg_profile", ".txt");
data_sg=dlmread(d_sg, ' ', 0, 0);

d_sw=strcat(char(data_name), "_sw_profile", ".txt");
data_sw=dlmread(d_sw, ' ', 0, 0);

d_T=strcat(char(data_name), "_T_profile", ".txt");
data_T=dlmread(d_T, ' ', 0, 0);

d_pg=strcat(char(data_name), "_pg_profile", ".txt");
data_pg=dlmread(d_pg, ' ', 0, 0);

d_pw=strcat(char(data_name), "_pw_profile", ".txt");
data_pw=dlmread(d_pw, ' ', 0, 0);

d_pc=strcat(char(data_name), "_pc_profile", ".txt");
data_pc=dlmread(d_pc, ' ', 0, 0);

d_pe=strcat(char(data_name), "_pe_profile", ".txt");
data_pe=dlmread(d_pe, ' ', 0, 0);

d_K=strcat(char(data_name), "_K_profile", ".txt");
data_K=dlmread(d_K, ' ', 0, 0);

d_sal=strcat(char(data_name), "_S_profile", ".txt");
data_sal=dlmread(d_sal, ' ', 0, 0);

d_zf=strcat(char(data_name), "_zf_profile", ".txt");
data_zf=dlmread(d_zf, ' ', 0, 0);

d_CeqHYD=strcat(char(data_name), "_Ceq_HYD_profile", ".txt");
data_CeqHYD=dlmread(d_CeqHYD, ' ', 0, 0);

d_CeqCH4=strcat(char(data_name), "_Ceq_CH4_profile", ".txt");
data_CeqCH4=dlmread(d_CeqCH4, ' ', 0, 0);

d_SO4_2n=strcat(char(data_name), "_SO4^2-_profile", ".txt");
data_SO4_2n=dlmread(d_SO4_2n, ' ', 0, 0);

d_NH4_1p=strcat(char(data_name), "_NH4^1+_profile", ".txt");
data_NH4_1p=dlmread(d_NH4_1p, ' ', 0, 0);

d_CH4=strcat(char(data_name), "_CH4_profile", ".txt");
data_CH4=dlmread(d_CH4, ' ', 0, 0);

d_CO2=strcat(char(data_name), "_CO2_profile", ".txt");
data_CO2=dlmread(d_CO2, ' ', 0, 0);

d_Cl_1n=strcat(char(data_name), "_Cl^1-_profile", ".txt");
data_Cl_1n=dlmread(d_Cl_1n, ' ', 0, 0);

d_OM=strcat(char(data_name), "_organic_matter_profile", ".txt");
data_OM=dlmread(d_OM, ' ', 0, 0);

d_CO3_2n=strcat(char(data_name), "_CO3^2-_profile", ".txt");
data_CO3_2n=dlmread(d_CO3_2n, ' ', 0, 0);

d_HCO3_1n=strcat(char(data_name), "_HCO3^1-_profile", ".txt");
data_HCO3_1n=dlmread(d_HCO3_1n, ' ', 0, 0);

%=====================
% POSTPROCESSED QOI
%=====================
potential1 = (data_pe - data_pg);

potential2 = (data_CeqCH4-data_CH4);

for n=1:1:size(t_new)
  %FIND BOTTOM OF STABILITY ZONE 
  BSR=0.;
  i=1;
  while (potential1(rt+n-1,cz+i)>0.1)
    BSR=z_new(i);
    i++;
  endwhile
  zBSR(n,1:columns(z_new)) = BSR;
  
  %FIND TOP OF EQUILIBRIUM ZONE
  TSZ=0.;
  j=columns(z_new);
  while (potential2(rt+n-1,cz+j-1)>0.)
    TSZ=z_new(j);
    j=j-1;
  endwhile
  zTSZ(n,1:columns(z_new)) = TSZ;
endfor


%=====================
k=1;
kmax=rows(t_new)

while (k<=kmax)

    k

    clf;
    figure('visible','off');
    set(gcf, 'PaperSize', [25 10], 'PaperPosition', [0 0 25 10]);
    nr_plots=1;   #no. of plot-rows
    nc_plots=5;   #plots per row
    %
    %PLOT: EQB AND GAS PRESSURES
    %=====================
    subplot (nr_plots, nc_plots, 1)
    Pmax=80.;
    xP = [0 : Pmax/(columns(z_new)-1) : Pmax]';
    y1=zBSR(k,:);                   #create lower curve
    y2=zTSZ(k,:);                   #create upper curve
    XP=[xP',fliplr(xP')];           #create continuous x value array for plotting
    Y=[y1,fliplr(y2)];            #create y values for out and then back
    fill(XP,Y,'y','facealpha',0.30,'LineStyle','none');                  #plot filled area
    hold on;
    plot(xP, zBSR(k,:), '-.k','linewidth',2);
    hold on;
    plot(xP, zTSZ(k,:), 'y','linewidth',1,'HandleVisibility','off');
    hold on;
    plot (data_pg(rt-1+k,cz:end), z_new, 'r','linewidth',4);
    hold on;
    plot (data_pe(rt-1+k,cz:end), z_new, 'b','linewidth',4);
    % add label to the x-axis
    xlabel('pressure [MPa]');
    % add label to the y-axis
    ylabel('depth below seafloor [m]');
    % add title for the plot
    #title(strcat("Time = ",int2str(t_new)," kyrs"),"fontsize", 30);
    % add legends
    h = legend('eqb.','BSR', 'pg', 'pe');
    legend (h, "location", "northeast");
    set (h, "fontsize", 12);
    set (h, 'box', 'off');
    % first 2 parameter sets the x-axis and next 2 will set the y-axis
    axis([0 Pmax 0 Z_bottom])
    set(gca,'YDir','reverse')
    set(gca, 'FontSize', 20)
    %
    % PLOT: SATURATIONS
    %=============
    subplot (nr_plots, nc_plots, 2)
    Smax=1.;
    xS = [0 : Smax/(columns(z_new)-1) : Smax]';
    XS=[xS',fliplr(xS')];       
    fill(XS,Y,'y','facealpha',0.30,'LineStyle','none');
    hold on;
    plot(xS, zBSR(k,:), '-.k','linewidth',2);
    hold on;
    plot(xS, zTSZ(k,:), 'y','linewidth',1,'HandleVisibility','off');
    hold on;
    plot (data_sg(rt-1+k,cz:end), z_new, 'r','linewidth',4);
    hold on;
    plot(data_sh(rt-1+k,cz:end), z_new, 'b','linewidth',4);
    % add label to the x-axis
    xlabel('saturation [-]');
    % add label to the y-axis
    #ylabel('depth below seafloor [m]');
    % add title for the plot
    #title(strcat("time=",int2str(t)," kyrs"));
    % add legends
    h = legend('eqb.','BSR','sg', 'sh');
    legend (h, "location", "northeast");
    set (h, "fontsize", 12);
    set (h, 'box', 'off');
    % first 2 parameter sets the x-axis and next 2 will set the y-axis
    axis([0 Smax 0 Z_bottom])
    set(gca,'YDir','reverse')
    set(gca, 'FontSize', 20)
    %
    % PLOT: CH4 CONC.
    %============
    subplot (nr_plots, nc_plots, 3)
    Cmax=600.;
    xC = [0 : Cmax/(columns(z_new)-1) : Cmax]';
    XC=[xC',fliplr(xC')];       
    fill(XC,Y,'y','facealpha',0.30,'LineStyle','none');
    hold on;
    plot(xC, zBSR(k,:), '-.k','linewidth',2);
    hold on;
    plot(xC, zTSZ(k,:), 'y','linewidth',1,'HandleVisibility','off');
    hold on;
    plot (data_CH4(rt-1+k,cz:end), z_new, 'r','linewidth',4);
    hold on;
    plot(data_CeqCH4(rt-1+k,cz:end), z_new, '-.r','linewidth',2);
    hold on;
    plot(data_CeqHYD(rt-1+k,cz:end), z_new, '-.b','linewidth',2);
    % add label to the x-axis
    xlabel('concentrations [mmol/Lpw]');
    % add label to the y-axis
    #ylabel('depth below seafloor [m]');
    % add title for the plot
    title(strcat("time =  ",num2str((t_new(end,1)-t_new(k,1)),'%.4f')," Ma BP"));
    % add legends
    h = legend('eqb.','BSR','CH4', 'CH4eq','HYDeq');
    legend (h, "location", "northeast");
    set (h, "fontsize", 12);
    set (h, 'box', 'off');
    % first 2 parameter sets the x-axis and next 2 will set the y-axis
    axis([0 Cmax 0 Z_bottom])
    set(gca,'YDir','reverse')
    set(gca, 'FontSize', 20)
    set(gca,'XTick',0:200:Cmax)
    %
    % PLOT: SALINITY
    %===========
    subplot (nr_plots, nc_plots, 4)
    SALmax=200.;
    xSAL = [0 : SALmax/(columns(z_new)-1) : SALmax]';
    XSAL=[xSAL',fliplr(xSAL')];       
    fill(XSAL,Y,'y','facealpha',0.30,'LineStyle','none');
    hold on;
    plot(xSAL, zBSR(k,:), '-.k','linewidth',2);
    hold on;
    plot(xSAL, zTSZ(k,:), 'y','linewidth',1,'HandleVisibility','off');
    hold on;
    plot (data_sal(rt-1+k,cz:end), z_new, 'k','linewidth',4);
    % add label to the x-axis
    xlabel('salinity [g/kg]',"fontsize", 16);
    % add label to the y-axis
    #ylabel('depth below seafloor [m]',"fontsize", 16);
    % add title for the plot
    #title(strcat("time=",int2str(t)," kyrs"));
    % add legends
    h = legend('eqb.','BSR','salinity');
    legend (h, "location", "northeast");
    set (h, "fontsize", 12);
    set (h, 'box', 'off');
    % first 2 parameter sets the x-axis and next 2 will set the y-axis
    axis([0 SALmax 0 Z_bottom])
    set(gca,'YDir','reverse')
    set(gca, 'FontSize', 20)
    %
    % PLOT: TEMPERATURE
    %==============
    subplot (nr_plots, nc_plots, 5)
    Tmax=310.0;
    Tmin=270.0;
    xT = [Tmin : (Tmax-Tmin)/(columns(z_new)-1) : Tmax]';
    XT=[xT',fliplr(xT')];       
    fill(XT,Y,'y','facealpha',0.30,'LineStyle','none');
    hold on;
    plot(xT, zBSR(k,:), '-.k','linewidth',2);
    hold on;
    plot(xT, zTSZ(k,:), 'y','linewidth',1,'HandleVisibility','off');
    hold on;
    plot (data_T(rt-1+k,cz:end), z_new,'k','linewidth',4);
    % add label to the x-axis
    xlabel('temperature [K]',"fontsize", 16);
    % add label to the y-axis
    #ylabel('depth below seafloor [m]',"fontsize", 16);
    % add title for the plot
    #title(strcat("time=",int2str(t)," kyrs"));
    % add legends
    h = legend('eqb.','BSR','T');
    legend (h, "location", "northeast");
    set (h, "fontsize", 12);
    set (h, 'box', 'off');
    % first 2 parameter sets the x-axis and next 2 will set the y-axis
    axis([Tmin Tmax 0 Z_bottom])
    set(gca,'YDir','reverse')
    set(gca, 'FontSize', 20)
    set(gca,'XTick',0:20:Tmax)

    %=====================
    % SAVE PLOTS
    %=====================
    name=strcat((output_name),"_profiles_", num2str(k,"%04d"), ".png");
    print (name);

    k=k+1;

endwhile

%=====================
% END
%=====================