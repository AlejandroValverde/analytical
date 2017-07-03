% Variation of parameters
clear all

%% Parameters
geom.L = 200; %mm
geom.B = 80; %mm
geom.H = 40; %mm
geom.t1 = 1;
geom.t2 = 1;

geom.nPointsPerSection = 100; %For analytical model

mat.E1 = 69000; %N/mm2, aluminium
mat.G1 = 26000; %N/mm2, aluminium: 26 GPa
mat.E2 = mat.E1/1; %N/mm2, steel: 200 GPa
mat.G2 = mat.E2 / ( 2*(0.3269 + 1) ); %N/mm2, steel: 79.3 GPa

loadCase.Q_z_total = 4000; %N

loadCase.posForceAdim = 0.5;

plotSettings.plotAnalytical = false;
plotSettings.savePlot = false;
plotSettings.MarkerSize = 30; %Marker size for scattered points, specified as a positive value in points.
plotSettings.LineWidth = 1.5; %Line width, specified as a positive value in points.
plotSettings.axGridAlpha = 0.2; %Grid-line transparency, specified as a value in the range [0,1].
plotSettings.axFontSize = 14; %Font size for axis labels, specified as a scalar numeric value.
plotSettings.axLineWidth = 1.5; %Width of axes outline, tick marks, and grid lines, specified as a scalar value in point units.
plotSettings.TitleFontSizeMultiplier = 1.5; %Scale factor for title font size, specified as a numeric value greater than 0.
%The axes applies this scale factor to the value of the FontSize property to determine the font size for the title.

%% Organize folder
[dirWork] = FsClass.organizeFolders();

%% Variation of B/H - Torsional stiffness
study.BoverH = linspace(0.5, 4, 10);
study.E1overE2 = [1, 10^0.5, 10^1, 10^1.5, 10^2, 10^2.5, 10^3, 10^3.5];
mat_init = mat;
geom_init = geom;

study.result_GIt_BoverH = zeros(length(study.E1overE2), length(study.BoverH));
study.result_ySCoverB_BoverH = zeros(length(study.E1overE2), length(study.BoverH));
study.Phi_y = zeros(length(study.E1overE2), length(study.BoverH));

%Update variable
for i_study= 1:length(study.E1overE2)
for j_study= 1:length(study.BoverH)
mat.E2 = mat.E1 / study.E1overE2(i_study);
mat.G2 = mat.E2 / (2*(0.3269 + 1) ); %N/mm2,
geom.B = geom.H * study.BoverH(j_study);
mainBeam %Execute analytical model script
study.result_GIt_BoverH(i_study, j_study) = oper.torStiff;
study.result_ySCoverB_BoverH(i_study, j_study) = oper.y_sc_closed / geom.B;
study.Phi_y(i_study, j_study) = oper.Phi_y * mat.E2;
end
end

figure('Units', 'normalized', 'Position', [0.15 0.1 0.7 0.75])
set(gcf, 'Name', 'Torsional stiffness variation for different cross sectional aspect ratio B/H and stiffness ratio E_1/E_2)')
ax = gca;
ylabel('G I_t [N mm^2]')
xlabel(['B/H'])
hold on

y_plots = zeros(1, length(study.E1overE2));
legendStr = cell(1, length(study.E1overE2));
for i= 1:length(study.E1overE2)
	y_plots(i) = plot(ax, study.BoverH, study.result_GIt_BoverH(i, :), 'LineWidth', plotSettings.LineWidth);
	legendStr{i} = ['E1/E2=10^{' num2str(round(log10(study.E1overE2(i)), 2)) '}'];
end


legend(ax, y_plots, legendStr, 'location','Best')

FsClass.SetAxisProp(ax, plotSettings);

figure('Units', 'normalized', 'Position', [0.15 0.1 0.7 0.75])
set(gcf, 'Name', 'Shear center for closed section location for different cross sectional aspect ratio B/H and stiffness ratio E_1/E_2)')
ax = gca;
ylabel('y_{SC} / B')
xlabel(['B/H'])
hold on

y_plots = zeros(1, length(study.E1overE2));
legendStr = cell(1, length(study.E1overE2));
for i= 1:length(study.E1overE2)
	y_plots(i) = plot(ax, study.BoverH, study.result_ySCoverB_BoverH(i, :), 'LineWidth', plotSettings.LineWidth);
	legendStr{i} = ['E1/E2=10^{' num2str(round(log10(study.E1overE2(i)), 2)) '}'];
end


legend(ax, y_plots, legendStr, 'location','Best')

FsClass.SetAxisProp(ax, plotSettings);

figure('Units', 'normalized', 'Position', [0.15 0.1 0.7 0.75])
set(gcf, 'Name', 'Shear center for closed section location for different cross sectional aspect ratio B/H and stiffness ratio E_1/E_2)')
ax = gca;
ylabel('\phi_y [N mm^2]')
xlabel(['B/H'])
hold on

y_plots = zeros(1, length(study.E1overE2));
legendStr = cell(1, length(study.E1overE2));
for i= 1:length(study.E1overE2)
	y_plots(i) = plot(ax, study.BoverH, study.Phi_y(i, :), 'LineWidth', plotSettings.LineWidth);
	legendStr{i} = ['E1/E2=10^{' num2str(round(log10(study.E1overE2(i)), 2)) '}'];
end


legend(ax, y_plots, legendStr, 'location','Best')

FsClass.SetAxisProp(ax, plotSettings);