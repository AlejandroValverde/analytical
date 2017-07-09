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

loadCase.Q_z_total = 2000; %N

loadCase.posForceAdim = 0.5;

plotSettings.plotAnalytical = false;
plotSettings.savePlot = true;

plotSettings.lineStyle = {'-', '--', ':', '-.'};
plotSettings.lineColor = {'k', 'b', 'r', 'y'};
plotSettings.MarkerSize = 30; %Marker size for scattered points, specified as a positive value in points.
plotSettings.LineWidth = 1.5; %Line width, specified as a positive value in points.
plotSettings.axGridAlpha = 0.2; %Grid-line transparency, specified as a value in the range [0,1].
plotSettings.axFontSize = 14; %Font size for axis labels, specified as a scalar numeric value.
plotSettings.axLineWidth = 1.5; %Width of axes outline, tick marks, and grid lines, specified as a scalar value in point units.
plotSettings.TitleFontSizeMultiplier = 1.5; %Scale factor for title font size, specified as a numeric value greater than 0.
%The axes applies this scale factor to the value of the FontSize property to determine the font size for the title.

%% Organize folder
[dirWork] = FsClass.organizeFolders();

%Save initial configuration
mat_init = mat;
geom_init = geom;

%%%%%%%%%%%%%%%%%%%%%%%
%% Variation of E1/E2
% study.E1overE2 = logspace(0,5,100);
study.n_E1overE2 = 300;
study.step_E1overE2 = 1.025;
study.E1overE2 = 1; %Initial value
study.twistTip = zeros(1, study.n_E1overE2);
study.storage_E1overE2 = zeros(1, study.n_E1overE2);
operCell = cell(1, study.n_E1overE2);

%Update variable
mat = mat_init;
geom = geom_init;
for i_study= 1:study.n_E1overE2
	study.storage_E1overE2(i_study) = study.E1overE2;
	mat.E2 = mat.E1 / study.E1overE2;
	mat.G2 = mat.E2 / (2*(0.3269 + 1) ); %N/mm2,
	mainBeam %Execute analytical model script
	study.twistTip(i_study) = twist_concentratedLoad(end) .* (180/pi);
	operCell{i_study} = oper;
	study.E1overE2 = study.E1overE2 * study.step_E1overE2;
end

figure('Units', 'normalized', 'Position', [0.15 0.1 0.7 0.75])
set(gcf, 'Name', 'Twist at tip as a function of stiffness ratio E_1/E_2')
ax = gca;
loglog(ax, study.storage_E1overE2, study.twistTip, 'LineWidth', plotSettings.LineWidth);
% set(ax, 'YScale', 'log')
% set(ax, 'XScale', 'log')
% set(ax, 'YLabel', '\phi_{tip} [deg]')
% set(ax, 'XLabel', 'E_1/E_2')
ylabel('\phi_{tip} [deg]')
xlabel(['E_1/E_2'])

FsClass.SetAxisProp(ax, plotSettings);

%Save figure
if plotSettings.savePlot
    saveas(gcf, [dirWork.figures 'twist-E1overE2.png'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variation of B/H
clearvars study mat geom
study.BoverH = linspace(0.5, 4, 100);
study.E1overE2 = [1, 10^0.5, 10^1, 10^1.5, 10^2, 10^2.5, 10^3];

study.result_GIt_BoverH = zeros(length(study.E1overE2), length(study.BoverH));
study.result_ySC_BoverH = zeros(length(study.E1overE2), length(study.BoverH));
study.Phi_y = zeros(length(study.E1overE2), length(study.BoverH));

%Update variable
mat = mat_init;
geom = geom_init;
for i_study= 1:length(study.E1overE2)
for j_study= 1:length(study.BoverH)
mat.E2 = mat.E1 / study.E1overE2(i_study);
mat.G2 = mat.E2 / (2*(0.3269 + 1) ); %N/mm2,
geom.B = (study.BoverH(j_study) * 120) / (study.BoverH(j_study) + 1);
geom.H = 120 - geom.B;
mainBeam %Execute analytical model script
study.result_GIt_BoverH(i_study, j_study) = oper.torStiff;
study.result_ySC_BoverH(i_study, j_study) = -oper.y_sc_closed / geom.B;
study.Phi_y(i_study, j_study) = oper.Phi_y;
end
end

%%%%%%%%%%%%%%%%%%%%%
%% Torsional stiffness

figure('Units', 'normalized', 'Position', [0.15 0.1 0.7 0.75])
set(gcf, 'Name', 'Torsional stiffness variation for different cross sectional aspect ratio B/H and stiffness ratio E_1/E_2')
ax = gca;
ylabel('G I_t [N mm^2]')
xlabel(['B/H'])
hold on

y_plots = zeros(1, length(study.E1overE2));
legendStr = cell(1, length(study.E1overE2));
i_plot = 1;
j_plot = 1;
for i= 1:length(study.E1overE2)
	y_plots(i) = plot(ax, study.BoverH, study.result_GIt_BoverH(i, :), ...
		'Color', plotSettings.lineColor{i_plot}, 'LineStyle',plotSettings.lineStyle{j_plot}, ...
		'LineWidth', plotSettings.LineWidth);
	legendStr{i} = ['E1/E2=10^{' num2str(round(log10(study.E1overE2(i)), 2)) '}'];
    if i_plot == 4
        i_plot = 1;
        j_plot = j_plot + 1;
    else
        i_plot = i_plot + 1;
    end
end


legend(ax, y_plots, legendStr, 'location','Best')

FsClass.SetAxisProp(ax, plotSettings);

%Save figure
if plotSettings.savePlot
    saveas(gcf, [dirWork.figures 'GIt-E1overE2-BoverH.png'])
end

%%%%%%%%%%%%%%%%%%%%%
%% Shear center

figure('Units', 'normalized', 'Position', [0.15 0.1 0.7 0.75])
set(gcf, 'Name', 'Shear center for closed section location for different cross sectional aspect ratio B/H and stiffness ratio E_1/E_2')
ax = gca;
ylabel('y_{SC} / B')
xlabel(['B/H'])
hold on

y_plots = zeros(1, length(study.E1overE2));
legendStr = cell(1, length(study.E1overE2));
i_plot = 1;
j_plot = 1;
for i= 1:length(study.E1overE2)
	y_plots(i) = plot(ax, study.BoverH, study.result_ySC_BoverH(i, :), ...
		'Color', plotSettings.lineColor{i_plot}, 'LineStyle',plotSettings.lineStyle{j_plot}, ...
		'LineWidth', plotSettings.LineWidth);
	legendStr{i} = ['E1/E2=10^{' num2str(round(log10(study.E1overE2(i)), 2)) '}'];
	if i_plot == 4
        i_plot = 1;
        j_plot = j_plot + 1;
    else
        i_plot = i_plot + 1;
    end
end


legend(ax, y_plots, legendStr, 'location','Best')

FsClass.SetAxisProp(ax, plotSettings);

%Save figure
if plotSettings.savePlot
    saveas(gcf, [dirWork.figures 'SC-E1overE2-BoverH.png'])
end

%%%%%%%%%%%%%%%%%%%%%
%% Flexural stiffness

figure('Units', 'normalized', 'Position', [0.15 0.1 0.7 0.75])
set(gcf, 'Name', 'Flexural stiffness for different cross sectional aspect ratio B/H and stiffness ratio E_1/E_2')
ax = gca;
ylabel('\phi_y [N mm^2]')
xlabel(['B/H'])
hold on

y_plots = zeros(1, length(study.E1overE2));
legendStr = cell(1, length(study.E1overE2));
i_plot = 1;
j_plot = 1;
for i= 1:length(study.E1overE2)
	y_plots(i) = plot(ax, study.BoverH, study.Phi_y(i, :), ...
		'Color', plotSettings.lineColor{i_plot}, 'LineStyle',plotSettings.lineStyle{j_plot}, ...
		'LineWidth', plotSettings.LineWidth); 
	legendStr{i} = ['E1/E2=10^{' num2str(round(log10(study.E1overE2(i)), 2)) '}'];
	if i_plot == 4
	    i_plot = 1;
	    j_plot = j_plot + 1;
	else
	    i_plot = i_plot + 1;
	end
end


legend(ax, y_plots, legendStr, 'location','Best')

FsClass.SetAxisProp(ax, plotSettings);

%Save figure
if plotSettings.savePlot
    saveas(gcf, [dirWork.figures 'EIy-E1overE2-BoverH.png'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variation of t_2/t_1
clearvars study mat geom
study.t2overt1 = linspace(0.1, 8, 100);
study.E1overE2 = [1, 10^0.5, 10^1, 10^1.5, 10^2, 10^2.5, 10^3];

study.result_GIt_t2overt1 = zeros(length(study.E1overE2), length(study.t2overt1));
study.result_ySC_t2overt1 = zeros(length(study.E1overE2), length(study.t2overt1));
study.Phi_y = zeros(length(study.E1overE2), length(study.t2overt1));

%Update variable
mat = mat_init;
geom = geom_init;
for i_study= 1:length(study.E1overE2)
for j_study= 1:length(study.t2overt1)
mat.E2 = mat.E1 / study.E1overE2(i_study);
mat.G2 = mat.E2 / (2*(0.3269 + 1) ); %N/mm2,
geom.t2 = (study.t2overt1(j_study) * 1) / (study.t2overt1(j_study) + 1);
geom.t1 = 1 - geom.t2;
mainBeam %Execute analytical model script
study.result_GIt_t2overt1(i_study, j_study) = oper.torStiff;
study.result_ySC_t2overt1(i_study, j_study) = -oper.y_sc_closed / geom.B;
study.Phi_y(i_study, j_study) = oper.Phi_y;
end
end

%%%%%%%%%%%%%%%%%%%%%
%% Torsional stiffness

figure('Units', 'normalized', 'Position', [0.15 0.1 0.7 0.75])
set(gcf, 'Name', 'Torsional stiffness variation for different thickness ratio t_2/t_1 and stiffness ratio E_1/E_2')
ax = gca;
ylabel('G I_t [N mm^2]')
xlabel(['t_2/t_1'])
hold on

y_plots = zeros(1, length(study.E1overE2));
legendStr = cell(1, length(study.E1overE2));
i_plot = 1;
j_plot = 1;
for i= 1:length(study.E1overE2)
	y_plots(i) = plot(ax, study.t2overt1, study.result_GIt_t2overt1(i, :), ...
		'Color', plotSettings.lineColor{i_plot}, 'LineStyle',plotSettings.lineStyle{j_plot}, ...
		'LineWidth', plotSettings.LineWidth);
	legendStr{i} = ['E1/E2=10^{' num2str(round(log10(study.E1overE2(i)), 2)) '}'];
	if i_plot == 4
	    i_plot = 1;
	    j_plot = j_plot + 1;
	else
	    i_plot = i_plot + 1;
	end
end


legend(ax, y_plots, legendStr, 'location','Best')

FsClass.SetAxisProp(ax, plotSettings);

%Save figure
if plotSettings.savePlot
    saveas(gcf, [dirWork.figures 'GIt-E1overE2-t2overt1.png'])
end

%%%%%%%%%%%%%%%%%%%%%
%% Shear center

figure('Units', 'normalized', 'Position', [0.15 0.1 0.7 0.75])
set(gcf, 'Name', 'Shear center for closed section location for different thickness ratio t_2/t_1 and stiffness ratio E_1/E_2')
ax = gca;
ylabel('y_{SC} / B')
xlabel(['t_2/t_1'])
hold on

y_plots = zeros(1, length(study.E1overE2));
legendStr = cell(1, length(study.E1overE2));
i_plot = 1;
j_plot = 1;
for i= 1:length(study.E1overE2)
	y_plots(i) = plot(ax, study.t2overt1, study.result_ySC_t2overt1(i, :), ...
		'Color', plotSettings.lineColor{i_plot}, 'LineStyle',plotSettings.lineStyle{j_plot}, ...
		'LineWidth', plotSettings.LineWidth);
	legendStr{i} = ['E1/E2=10^{' num2str(round(log10(study.E1overE2(i)), 2)) '}'];
	if i_plot == 4
        i_plot = 1;
        j_plot = j_plot + 1;
    else
        i_plot = i_plot + 1;
    end
end


legend(ax, y_plots, legendStr, 'location','Best')

FsClass.SetAxisProp(ax, plotSettings);

%Save figure
if plotSettings.savePlot
    saveas(gcf, [dirWork.figures 'SC-E1overE2-t2overt1.png'])
end

%%%%%%%%%%%%%%%%%%%%%
%% Flexural stifness

figure('Units', 'normalized', 'Position', [0.15 0.1 0.7 0.75])
set(gcf, 'Name', 'Flexural stiffness for different thickness ratio t_2/t_1 and stiffness ratio E_1/E_2')
ax = gca;
ylabel('\phi_y [N mm^2]')
xlabel(['t_2/t_1'])
hold on

y_plots = zeros(1, length(study.E1overE2));
legendStr = cell(1, length(study.E1overE2));
i_plot = 1;
j_plot = 1;
for i= 1:length(study.E1overE2)
	y_plots(i) = plot(ax, study.t2overt1, study.Phi_y(i, :), ...
		'Color', plotSettings.lineColor{i_plot}, 'LineStyle',plotSettings.lineStyle{j_plot}, ...
		'LineWidth', plotSettings.LineWidth);
	legendStr{i} = ['E1/E2=10^{' num2str(round(log10(study.E1overE2(i)), 2)) '}'];
	if i_plot == 4
        i_plot = 1;
        j_plot = j_plot + 1;
    else
        i_plot = i_plot + 1;
    end
end


legend(ax, y_plots, legendStr, 'location','Best')

FsClass.SetAxisProp(ax, plotSettings);

%Save figure
if plotSettings.savePlot
    saveas(gcf, [dirWork.figures 'EIy-E1overE2-t2overt1.png'])
end