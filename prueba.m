% study.E1overE2 = [1, 20, 50];
study.E1overE2 = linspace(1, 100, 10);

figure('Units', 'normalized', 'Position', [0.15 0.1 0.7 0.75])
set(gcf, 'Name', 'Comparison of results')
ax = gca;
hold on
subPlotHandle = cell(1, 4);
for e=1:length(subPlotHandle)
    
    subPlotHandle{e} = subplot(length(subPlotHandle), 1, e);
    hold on
    FsClass.SetAxisProp(subPlotHandle{e}, plotSettings);
end
for i_study = 1:length(study.E1overE2)

    %Update variable
    mat.E2 = mat.E1 / study.E1overE2(i_study);
    mat.G2 = mat.E2 / (2*(0.3269 + 1) ); %N/mm2,
	mainBeam %Execute analytical model script
    
    scatter(subPlotHandle{1}, study.E1overE2(i_study), oper.torStiff, 'b')
    scatter(subPlotHandle{2}, study.E1overE2(i_study), oper.y_sc_closed, 'r')
    momentAnalytical = (loadCase.q_z .* geom.L) .* (y_load - oper.y_sc_closed);
    scatter(subPlotHandle{3}, study.E1overE2(i_study), momentAnalytical, 'sk')
    
    %For abaqus model
    forceMagnitude = -4000; %Total force
    forceRelXPos = 0.5;
    x_load = -((B/2) - (forceRelXPos*B)); %The abaqus and the analytical model have opposite signs for the axes
    x_sc = -oper.y_sc_closed;
    x_moment = x_load - x_sc;
    moment = x_moment * forceMagnitude;
    
    scatter(subPlotHandle{3}, study.E1overE2(i_study), moment, 'r')
    
    %Twist from analytical model+
    twistAnalytical = (momentAnalytical/oper.torStiff) * geom.L;
    scatter(subPlotHandle{4}, study.E1overE2(i_study), twistAnalytical * (180/pi), 'k')
end

xlabel('E_1/E_2')

ylabel(subPlotHandle{1}, 'G I_t')
ylabel(subPlotHandle{2}, 'y_{sc, closed}')
ylabel(subPlotHandle{3}, 'M_{t}')
ylabel(subPlotHandle{4}, '\theta [deg]')

linkaxes([subPlotHandle{1},subPlotHandle{2},subPlotHandle{3}, subPlotHandle{4}],'x')
