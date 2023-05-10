clear;
%=====================
# USER INPUTS
%=====================
K="1e-16";
kr="1e-18";
vs=0.05;
data_name   =strcat("../test_K=",K,"_kr=",kr,"_vs=",num2str(vs,'%.4f'));
output_name =strcat("test_K=",K,"_kr=",kr,"_vs=",num2str(vs,'%.4f'));
t_END       = 27.0e6;   #annum #total calculation time 
t_recording = 0.0e6;    #annum #time when fine output is recorded
t_POI       = 26.0e6;    #annum #time when period of interest begins
t_offset    = 26.0e6;   #annum #offsets time axis s.t. t=0 corresponsts to t_offset
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
potential = (data_pg - data_pe + 0.1 );
GHSZ=potential./abs(potential);

[dz_T dt_T] = gradient(data_T);
[dz_pg dt_pg] = gradient(data_pg);

potential1 = (data_pe - data_pg);

potential2 = (data_CeqCH4-data_CH4);

for n=1:1:size(t_new)
  %FIND BOTTOM OF STABILITY ZONE 
  BSR(n)=0.;
  i=1;
  while (potential1(rt+n-1,cz+i)>0.1)
    BSR(n)=z_new(i);
    i++;
  endwhile
  zBSR(n,1:columns(z_new)) = BSR(n);
  
  %FIND TOP OF EQUILIBRIUM ZONE
  TSZ(n)=0.;
  j=columns(z_new);
  while (potential2(rt+n-1,cz+j-1)>0.)
    TSZ(n)=z_new(j);
    j=j-1;
  endwhile
  zTSZ(n,1:columns(z_new)) = TSZ(n);
endfor

%=====================
% AUTOCORRELATION
%=====================
% 
for n=1:1:size(t_new)
  %sh_max
  tmp = max(data_sh(rt+n-1,cz:end));
  sh_max(n) = round(100 * tmp) / 100;
  %total hydrate volume
  tmp = sum(data_sh(rt-1+n,cz:end));
  sh_vol(n) = round(100 * tmp) / 100;
  %base of GHSZ
  bGHSZ(n)=0.0;
  for k=2:1:columns(z_new)-1
    if( abs(GHSZ(rt+n-1,cz+k)-GHSZ(rt+n-1,cz+k-1)) >0.0 )
      bGHSZ(n) = 0.5*(z_new(1,k)+z_new(1,k-1));
      k=columns(z_new);
    elseif( isnan(GHSZ(rt+n-1,cz+k)) )
      bGHSZ(n) = z_new(1,k);
      k=columns(z_new);
    endif
  endfor
endfor

if( (max(sh_vol)-min(sh_vol)) > 0.001 )
  [P,f,alpha] = fastlomb(sh_vol,t_new);
  [maxpower,maxpowerindex] = max(P);
  rightfreq = f(maxpowerindex);
  period = 1/rightfreq;
else 
  period = NaN;
endif


printf('QOI: \t (i) peak sh[-]:[%.3f %.3f] DIFF:%.3f \t (ii) sh volume[m^2/m^3]:[%.2f %.2f] DIFF:%.2f \t (iii) bGHSZ[mbsf]:[%.2f %.2f] DIFF:%.2f \t (iv) Time period: %.2f years \n', 
        min(sh_max), max(sh_max), (max(sh_max)-min(sh_max)), 
        min(sh_vol), max(sh_vol), (max(sh_vol)-min(sh_vol)), 
        min(bGHSZ), max(bGHSZ), (max(bGHSZ)-min(bGHSZ)), 
        period 
);

%=====================
% PLOTTING
%=====================
% 
clf;
%

generate_cyclic_profiles( 1, 
                          sh_max, sh_vol, bGHSZ,
                          1.e-6*(t_END-t_offset)-t_new,
                          Z_bottom,
                          strcat("K_0=",K,", k_r=",kr,", v_s=",num2str(vs,'%.2f')),
                          strcat(char(output_name),"_cyclic_profile", ".png"),
                          "on"
                         );

%
generate_timeseries_plot( 2, 
                          data_sh(rt:end,cz:end)', 
                          z_new, 
                          1.e-6*(t_END-t_offset)-t_new, 
                          "Hydrate Saturation",
                          "depth [mbsf]", 
                          "time [Ma BP]", 
                          "[-]",
                          strcat(char(output_name),"_sh_profile", ".png")
                         );
                         
generate_timeseries_plot( 3,
                          data_sg(rt:end,cz:end)', 
                          z_new, 
                          1.e-6*(t_END-t_offset)-t_new, 
                          "Free Gas Saturation",
                          "depth [mbsf]", 
                          "time [Ma BP]", 
                          "[-]",
                          strcat(char(output_name),"_sg_profile", ".png")
                         );
                         
                         
generate_timeseries_plot( 4,
                          GHSZ(rt:end,cz:end)', 
                          z_new, 
                          1.e-6*(t_END-t_offset)-t_new, 
                          "Gas Hydrate Stability Zone",
                          "depth [mbsf]", 
                          "time [Ma BP]", 
                          "[-]",
                          strcat(char(output_name),"_GHSZ_profile", ".png")
                         );
                         
generate_timeseries_plot( 5,
                          data_pg(rt:end,cz:end)', 
                          z_new, 
                          1.e-6*(t_END-t_offset)-t_new, 
                          "Gas Phase Pressure",
                          "depth [mbsf]", 
                          "time [Ma BP]", 
                          "[MPa]",
                          strcat(char(output_name),"_pg_profile", ".png")
                         );
                         
generate_timeseries_plot( 6,
                          dt_pg(rt:end,cz:end)', 
                          z_new, 
                          1.e-6*(t_END-t_offset)-t_new, 
                          "d_t Pg",
                          "depth [mbsf]", 
                          "time [Ma BP]", 
                          "[MPa/a]",
                          strcat(char(output_name),"_dtpg_profile", ".png")
                         );
                         
generate_timeseries_plot( 7,
                          data_T(rt:end,cz:end)', 
                          z_new, 
                          1.e-6*(t_END-t_offset)-t_new, 
                          "Temperature",
                          "depth [mbsf]", 
                          "time [Ma BP]", 
                          "[^oC]",
                          strcat(char(output_name),"_T_profile", ".png")
                         );

generate_timeseries_plot( 8,
                          dt_T(rt:end,cz:end)', 
                          z_new, 
                          1.e-6*(t_END-t_offset)-t_new, 
                          "d_t T",
                          "depth [mbsf]", 
                          "time [Ma BP]", 
                          "[^oC/a]",
                          strcat(char(output_name),"_dtT_profile", ".png")
                         );
                        

%=====================
# END
%=====================