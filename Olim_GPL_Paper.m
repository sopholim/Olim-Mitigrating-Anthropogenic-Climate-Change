%% Mitigating Anthropogenic Climate Change with Aqueous Green Energy
% Sophia Olim
% May 16 2023
% Define Variables
% 100 years using diagnosed emissions from IPCC RCP 8.5

% Read ncfiles
ncfile_rcp = 'tsi.rcp.100yr.nc';
ncfile_3TW_bm = 'tsi.02001.01.01.nc';

%Time
time1 = ncread('tsi.rcp.100yr.nc', 'time'); %300 yrs
% gws - global warming scenario with diagnosed emissions from IPCC
% RCP 8.5, no OTEC or DACCS
rcp_ACO2 = ncread(ncfile_rcp, 'A_co2');
rcp_Asat = ncread(ncfile_rcp, 'A_sat');
rcp_Fco2emit = ncread(ncfile_rcp, 'F_co2emit');

% 3 TW power
bm_3TW_ACO2 = ncread(ncfile_3TW_bm, 'A_co2');
bm_3TW_Asat = ncread(ncfile_3TW_bm, 'A_sat');
% global flux co2 emissions
bm_3TW_Fco2emit = ncread(ncfile_3TW_bm, 'F_co2emit');
% Ocean 
bm_3TW_Otempsur = ncread(ncfile_3TW_bm, 'O_tempsur');
bm_3TW_Ophyt = ncread(ncfile_3TW_bm, 'O_phyt');
bm_3TW_Oalksur = ncread(ncfile_3TW_bm, 'O_alksur');
bm_3TW_Oo2sur= ncread(ncfile_3TW_bm, 'O_o2sur');
% OTEC power
bm_3TW_OTECpwrtot = ncread(ncfile_3TW_bm, 'OT_pwrtot');
% OTEC emt red
bm_3TW_OTECemtred = ncread(ncfile_3TW_bm, 'OT_emtred');

% difference files
% OTEC - gws (change - control)
% ACO2
diff_bm_3TW_ACO2 = bm_3TW_ACO2 - rcp_ACO2;
diff_bm_3TW_Asat = bm_3TW_Asat - rcp_Asat; 
% F co2 emit
diff_bm_3TW_Fco2emit = bm_3TW_Fco2emit - rcp_Fco2emit;

%% Asat
figure(1);clf;
set(gca,'fontsize',14);
newcolors = [229/255 56/255 59/255; 157/255 2/255 8/255];
colororder(newcolors);
plot(time1, bm_3TW_Asat, ':', time1, rcp_Asat, '-', time1, rcp_Asat, 'ko', 'MarkerIndices', 31, LineWidth=6);
ax = gca;
ax.FontSize = 24;
%set(gca,'ycolor',[229/255 56/255 59/255]);
ylabel('Surface Temperature [^oC]');
xlabel('Time [years]');
legend('With OTEC & DACCS', 'No OTEC', 'OTEC begins in 2030', Location='northwest', fontsize = 30);
title('Global Average Surface Temperature', FontSize=24);
%% F co2 emit
figure(2);
newcolors = [4/255 150/255 255/255; 0/255 78/255 137/255];
colororder(newcolors);
plot(time1, bm_3TW_Fco2emit, ':', time1, rcp_Fco2emit, '-', time1, rcp_Fco2emit, 'ko','MarkerIndices',31, LineWidth=6);
ax = gca;
ax.FontSize = 24;
ylabel('Total CO_2 Emissions [Gt/yr]', fontsize = 24)
xlabel('Time [years]');
legend('With OTEC & DACCS', 'No OTEC', 'OTEC begins in 2030', Location='northwest', fontsize = 30);
%set(gca,'ycolor',[4/255 150/255 255/255]);
title('Global Total CO_2 Emissions', FontSize=24);
%% A co2
figure(3);
newcolors = [123/255 44/255 191/255; 60/255 9/255 108/255];
colororder(newcolors);
plot(time1, bm_3TW_ACO2, ':', time1, rcp_ACO2, '-', time1, rcp_ACO2, 'ko','MarkerIndices',31, LineWidth=6);
ax = gca;
ax.FontSize = 24;
ylabel('Atmospheric CO_2 [ppm]', fontsize = 24)
xlabel('Time [years]');
legend('With OTEC & DACCS', 'No OTEC', 'OTEC begins in 2030', Location='northwest', fontsize = 30);
%set(gca,'ycolor',[4/255 150/255 255/255]);
title('Global Atmospheric CO_2', FontSize=24);
%% Gt of co2 sequestered
%OTEC 3TW CO_2 flux - no intervention CO_2 flux
diff_Fco2emit = bm_3TW_Fco2emit - rcp_Fco2emit;
x = time1;
y = abs(diff_Fco2emit);
figure(4);
ax = gca;
ax.FontSize = 24;
plot(x,y,'b', LineWidth=6);
ylabel('[Gt]', fontsize = 24);
title('Amount of total CO_2', fontsize = 24);
xlabel('[years]', fontsize = 24);
% A = trapz takes sums the area under the curve
A = trapz(x,y);
legend('19.6 Gt of CO_2 reduced', Location='northwest', fontsize = 24);
%% temp diff
%OTEC 3TW A sat - no intervention A sat
diff_Asat = bm_3TW_Asat - rcp_Asat;
m = time1;
n = abs(diff_bm_3TW_Asat);
figure(5);
ax = gca;
ax.FontSize = 24;
plot(m,n,'r', LineWidth=6);
ylabel('[^oC]', fontsize = 24);
title('Temperature Diff', fontsize = 24);
xlabel('[years]', fontsize = 24);
% A = trapz takes sums the area under the curve
A = trapz(m,n);
legend('1.18 ^oC Temp Diff', Location='northwest', fontsize = 24);

%% A c02 diff
%OTEC 3TW A co2 - no intervention A co2
diff_Aco2 = bm_3TW_ACO2 - rcp_ACO2;
a = time1;
b = abs(diff_bm_3TW_ACO2);
figure(6);
ax = gca;
ax.FontSize = 24;
plot(a,b,'m', LineWidth=6);
ylabel('[ppm]', fontsize = 24);
title('Atmos CO_2 Diff', fontsize = 24);
xlabel('[years]', fontsize = 24);
% A = trapz takes sums the area under the curve
A = trapz(m,n);
legend('277 ppm Diff', Location='northwest', fontsize = 24);

