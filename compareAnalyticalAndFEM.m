%Import data from Abaqus
clear all
% close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
geom.L = 800; %mm
geom.B = 80; %mm
geom.H = 40; %mm
geom.t1 = 1;
geom.t2 = 1;

geom.nPointsPerSection = 100; %For analytical model


geom.nInnerRibs = 1; %For the abaqus model
mat.E1 = 69000; %N/mm2, aluminium
mat.G1 = 26000; %N/mm2, aluminium: 26 GPa
mat.E2 = mat.E1/100; %N/mm2, steel: 200 GPa
mat.G2 = mat.E2 / ( 2*(0.3269 + 1) ); %N/mm2, steel: 79.3 GPa

% Real materials
% mat.E1 = 200000.0; %69000; %N/mm2, aluminium
% mat.G1 = 79300; %N/mm2, aluminium: 26 GPa
% mat.E2 = 69000.0; % mat.E1/10; %N/mm2, steel: 200 GPa
% mat.G2 = 26000; %mat.E2 / ( 2*(0.3269 + 1) ); %N/mm2, steel: 79.3 GPa

loadCase.Q_z_total = 2000; %N

loadCase.posForceAdim = 0.0;

%% Plotting settings
plotSettings.plotAnalytical = false;
plotSettings.plotDistributedLoad = false;
plotSettings.plotBending = false;
plotSettings.plotTwistAlongZ = false;
plotSettings.plotParametricStudy = true;
plotSettings.shearCenterPos = false;
options.executeAbaqus = false;

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

%% Plot analytical

if plotSettings.plotAnalytical

	mainBeam

end

%% Shear center position shift with G2 and t2

study.G1overG2 = linspace(1, 50, 10);
study.t1overt2 = linspace(1, 50, 10);

operCell_y_sc = cell(1, length(study.t1overt2));

if plotSettings.shearCenterPos
	for t =1:length(study.t1overt2)

		operCell_y_sc{t} = zeros(1, length(study.G1overG2));

		for G = 1:length(study.G1overG2)

			%Update variable
			mat.G2 = mat.G1 / study.G1overG2(G);
			geom.t2 = geom.t1 / study.t1overt2(t);
			mainBeam %Execute analytical model script

			operCell_y_sc{t}(G) = oper.y_sc_closed;

		end
	end

	figure('Units', 'normalized', 'Position', [0.15 0.1 0.7 0.75])
	set(gcf, 'Name', 'Twist due to distributed load')
	ax = gca;
	ylabel('G_1/G_2')
	xlabel(['y_{sc} [mm], ' 'B=' num2str(geom.B) 'mm'])
	title(['Closed section shear center y_{sc} shifting for G_1=' num2str(mat.G1) 'N/mm^2 and t_1=' num2str(geom.t1) 'mm'])
	hold on

	y_plots = zeros(1, length(study.t1overt2));
	legendStr = cell(1, length(study.t1overt2));
	for t =1:length(study.t1overt2)
		y_plots(t) = plot(ax, operCell_y_sc{t}, study.G1overG2, '.', 'MarkerSize', plotSettings.MarkerSize);
		legendStr{t} = ['t1/t2=' num2str(study.t1overt2(t))];
	end

	legend(ax, y_plots, legendStr, 'location','Best')

	FsClass.SetAxisProp(ax, plotSettings);

end

%% Twist along z
Iter = 0;
if plotSettings.plotTwistAlongZ

	cd(dirWork.main)
	mainBeam %Execute analytical model script

	%Update values shear center location and other parameters
	fprintf(['Writing Abaqus input to file...' '\n']) %Show progress
	FsClass.writeToFile('inputAbaqusAnalytical.txt', dirWork, {'E1overE2', 'y_sc', 'numberOfRibs'}, {mat.E1/mat.E2, oper.y_sc_closed, geom.nInnerRibs}, Iter)
	
	%Execute abaqus
	if options.executeAbaqus
		fprintf(['Running Abaqus...' '\n']) %Show progress
		cd(dirWork.abaqus)
		[status,cmout] = dos(strcat('abaqus cae noGUI=basicBeamScript.py'));
		assert(status == 0, ['Abaqus execution failed: ' cmout])
	end

	cd(dirWork.main)
	[dirWork] = FsClass.organizeFolders(num2str(Iter));
	[xPos, dataU2, dataUR3] = FsClass.loadAbaqus(dirWork, num2str(Iter));

	if plotSettings.plotDistributedLoad
		figure('Units', 'normalized', 'Position', [0.15 0.1 0.7 0.75])
		set(gcf, 'Name', 'Twist due to distributed load')
		ax_distributed = gca;
		xlabel('z/L')
		ylabel('\theta_{tip} [deg]')
		hold on
		FsClass.SetAxisProp(ax_distributed, plotSettings);
	end

	if plotSettings.plotBending
		figure('Units', 'normalized', 'Position', [0.15 0.1 0.7 0.75])
		set(gcf, 'Name', 'Beam bending')
		ax_bending = gca;
		xlabel('z/L')
		ylabel('w_{b} [mm]')
		hold on
		FsClass.SetAxisProp(ax_bending, plotSettings);
	end

	figure('Units', 'normalized', 'Position', [0.15 0.1 0.7 0.75])
	set(gcf, 'Name', ['Twist due to moment applied on the rib for E1/E2:' num2str(mat.E1/mat.E2) ' and ' num2str(geom.nInnerRibs) ' inner ribs'])
	ax_concentrated = gca;
	xlabel('z/L')
	ylabel('\theta_{tip} [deg]')
	hold on
	FsClass.SetAxisProp(ax_concentrated, plotSettings);

	for point =1:length(xPos)

	    twist = ((dataU2{point}(end) - dataU2{point}(1)) / geom.B) * (180/pi);
	    bending = dataU2{point}(floor(length(dataU2{point})/2));%mean( dataU2{point} );

	    if point == 1
	        y2 = plot(ax_concentrated, xPos(point) / geom.L, twist,'.r');
	        if plotSettings.plotDistributedLoad
	        	y1 = plot(ax_distributed, xPos(point) / geom.L, twist,'.b');
	        end
	        if plotSettings.plotBending
	        	y1_bend = plot(ax_bending, xPos(point) / geom.L, bending,'.r');
	        end
	    else
	        plot(ax_concentrated, xPos(point) / geom.L, twist,'.r')
	        if plotSettings.plotDistributedLoad
	        	plot(ax_distributed, xPos(point) / geom.L, twist,'.b')
	        end
	        if plotSettings.plotBending
	        	plot(ax_bending, xPos(point) / geom.L, bending,'.r');
	        end
	    end

	end

	meanTwist = (dataUR3.ur3_up + dataUR3.ur3_dn + dataUR3.ur3_ri + dataUR3.ur3_lf) ./4;

	if plotSettings.plotDistributedLoad
		y3 = plot(ax_distributed, dataUR3.zAdim, dataUR3.ur3_up .* (180/pi),'^k');
		y7 = plot(ax_distributed, dataUR3.zAdim, dataUR3.ur3_dn .* (180/pi),'vk');
		y9 = plot(ax_distributed, dataUR3.zAdim, dataUR3.ur3_ri .* (180/pi),'sk');
		y11 = plot(ax_distributed, dataUR3.zAdim, dataUR3.ur3_lf .* (180/pi),'dk');
		y13 = plot(ax_distributed, dataUR3.zAdim, meanTwist .* (180/pi),'ok');
	end

	y4 = plot(ax_concentrated, dataUR3.zAdim, dataUR3.ur3_up .* (180/pi),'^k');
	y8 = plot(ax_concentrated, dataUR3.zAdim, dataUR3.ur3_dn .* (180/pi),'vk');
	y10 = plot(ax_concentrated, dataUR3.zAdim, dataUR3.ur3_ri .* (180/pi),'sk');
	y12 = plot(ax_concentrated, dataUR3.zAdim, dataUR3.ur3_lf .* (180/pi),'dk');
	y14 = plot(ax_concentrated, dataUR3.zAdim, meanTwist .* (180/pi),'ok');

	if plotSettings.plotDistributedLoad
		y5 = plot(ax_distributed, dataUR3.zAdim, twist_fun(dataUR3.zAdim .* geom.L) .* (180/pi), 'b');
	end
	y6 = plot(ax_concentrated, xAdimSec, twist_concentratedLoad .* (180/pi), 'b');

	legend(ax_concentrated, [y2 y4 y8 y10 y12 y14 y6], 'FEM-U2','FEM-UR3_{up}', 'FEM-UR3_{down}', 'FEM-UR3_{right}', 'FEM-UR3_{left}', 'FEM-UR3_{mean}', 'Analytical', 'location','Best')
	if plotSettings.plotDistributedLoad
		legend(ax_distributed, [y1 y3 y7 y9 y11 y13 y5], 'FEM-U2','FEM-UR3_{up}', 'FEM-UR3_{down}', 'FEM-UR3_{right}', 'FEM-UR3_{left}', 'FEM-UR3_{mean}', 'Analytical', 'location','Best')
	end

	if plotSettings.plotBending
		y2_bend = plot(ax_bending, xAdimSec, -oper.w_b, 'b'); %Need to change sign because different axis

		legend(ax_bending, [y1_bend y2_bend], 'FEM-U2', 'Analytical', 'location','Best')
	end

end

cd(dirWork.main)

%% Parametric study

if plotSettings.plotParametricStudy

%Study for range of E2
% Parameter for abaqus: E1overE2
% Values
study.E1overE2 = linspace(1, 30, 5);
% study.E1overE2 = logspace(0,5,5);
% study.E1overE2 = [1, 20, 10];

twistTipAbaqus_from_U2 = zeros(1, length(study.E1overE2));
twistTipAbaqus_from_UR3 = zeros(1, length(study.E1overE2));
twistTipAbaqus_from_UR3_mean = zeros(1, length(study.E1overE2));
twistTipAnalytical = zeros(1, length(study.E1overE2));
operCell = cell(1, length(study.E1overE2));

for i_study = 1:length(study.E1overE2)

    %Update variable
    mat.E2 = mat.E1 / study.E1overE2(i_study);
    mat.G2 = mat.E2 / (2*(0.3269 + 1) ); %N/mm2,
	mainBeam %Execute analytical model script
	operCell{i_study} = oper;
	twistTipAnalytical(i_study) = twist_concentratedLoad(end) .* (180/pi); %.* dataUR3.zAdim(end)
	
	%Update input Abaqus input file
	fprintf(['Writing Abaqus input to file...' '\n']) %Show progress
	FsClass.writeToFile('inputAbaqusAnalytical.txt', dirWork, {'E1overE2', 'y_sc', 'numberOfRibs'}, {mat.E1/mat.E2, oper.y_sc_closed, geom.nInnerRibs}, i_study)

	%Execute Abaqus
	if options.executeAbaqus
		fprintf(['Running Abaqus...' '\n']) %Show progress
		cd(dirWork.abaqus)
		[status,cmout] = dos(strcat('abaqus cae noGUI=basicBeamScript.py'));
		assert(status == 0, ['Abaqus execution failed: ' cmout])
	end

	%Update folder with results
	cd(dirWork.main)
	[dirWork] = FsClass.organizeFolders(num2str(i_study));

	%Load code from Abaqus
	fprintf(['Loading results from Abaqus...' '\n']) %Show progress
	[xPos, dataU2, dataUR3] = FsClass.loadAbaqus(dirWork, num2str(i_study));

	%%%Twist from U2
	% Find index for last position 
	[~, index] = find(xPos == max(xPos));
	twistTipAbaqus_from_U2(i_study) = ((dataU2{index}(end) - dataU2{index}(1)) / geom.B) * (180/pi);
	
	%%%Twist from UR3
	twistTipAbaqus_from_UR3(i_study) = dataUR3.ur3_up(end) .* (180/pi); %Upper flange
	% From mean UR3
	meanTwist = (dataUR3.ur3_up + dataUR3.ur3_dn + dataUR3.ur3_ri + dataUR3.ur3_lf) ./4;
	[meanTwistTip] = FsClass.finalGuessOfOverallTwistAtTip(dataUR3.zAdim, meanTwist .* (180/pi));
	twistTipAbaqus_from_UR3_mean(i_study) = meanTwistTip;
	
	fprintf(['Executing analytical model...' '\n']) %Show progress

	fprintf(['-> End of iteration ' num2str(i_study) '\n' '\n']) %Show progress

end

figure('Units', 'normalized', 'Position', [0.15 0.1 0.7 0.75])
set(gcf, 'Name', 'Twist for different values E_1/E_2')

ax = gca;

xlabel('E_1/E_2')
ylabel('\theta_{tip} [deg]')

% set(ax,'XScale','log');

hold on

y1 = plot(ax, study.E1overE2, twistTipAbaqus_from_U2, '.k', 'MarkerSize', plotSettings.MarkerSize);
y2 = plot(ax, study.E1overE2, twistTipAbaqus_from_UR3, '.r', 'MarkerSize', plotSettings.MarkerSize);
y3 = plot(ax, study.E1overE2, twistTipAnalytical, 'b', 'LineWidth', plotSettings.LineWidth);
y4 = plot(ax, study.E1overE2, twistTipAbaqus_from_UR3_mean, 'ob', 'LineWidth', plotSettings.LineWidth);

legend(ax, [y1 y2 y4 y3], 'FEM-U2','FEM-UR3_{up}', 'FEM-UR3_{mean}', 'Analytical', 'location','Best')

FsClass.SetAxisProp(ax, plotSettings);

if plotSettings.savePlot
    saveas(gcf, ['twist-E1OverE2.png'])
end

for j=1:length(operCell)
    
    j
    operCell{j}.Phi_y
    operCell{j}.y_sc_closed
    operCell{j}.torStiff
    
end

end