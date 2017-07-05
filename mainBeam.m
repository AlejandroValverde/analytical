% MAIN SCRIPT
% Analytical calculation of the displacement of the beam
%
%   Original code by Alejandro Valverde: 18 Mar. 2017
%   Last modified by Alejandro Valverde: 28 Mar. 2017
%
% close all
% clear all

%This parameters are assumed to be entered in the model
%L, B, t1, t2, nPointsPerSection, E1, E2, G1, G2, Q_z_total, posForceAdim

%Pre-calculations
A = geom.H*geom.B;
q_z = loadCase.Q_z_total/geom.L; %N/mm
y_load = (geom.B/2) - (loadCase.posForceAdim*geom.B); %mm

loadCase.q_z = q_z;
loadCase.y_load = y_load;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Operations
t1 = geom.t1;
t2 = geom.t2;
H = geom.H;
B = geom.B;
nPointsPerSection = geom.nPointsPerSection;
E1 = mat.E1;
E2 = mat.E2;
G1 = mat.G1;
G2 = mat.G2;
Q_z = 1000; %The value inserted here does not play a role
y_load = loadCase.y_load;

A = H*B;

%Phi_y
Phi_y_1And5 = E1*(1/12)*t1*H^3;
Phi_y_2 = E1*( ((1/12)*B*t1^3) + ((B*t1*(H/2)^2)) );
Phi_y_3 = E2*(1/12)*t2*H^3;

Phi_y = Phi_y_1And5 + (2*Phi_y_2) + Phi_y_3;

oper.Phi_y = Phi_y;

%Shear flow functions
q1 = @(s) - (Q_z./Phi_y) .* E1 .* t1 .* (s.^2./2);
q2 = @(s) - (Q_z./Phi_y) .* E1 .* t1 .* ((H.^2./8) + (s.*H./2));
q3 = @(s) - (Q_z./Phi_y) .* ( (E1 .* t1 .* ((H.^2./8) + (B.*H./2))) + (E2 .* t2 .* ((s.*H./2) - (s.^2./2))) );
q4 = @(s) - (Q_z./Phi_y) .* E1 .* t1 .* ((H.^2./8) + (B.*H./2) - (s.*H./2));
q5 = @(s) - (Q_z./Phi_y) .* E1 .* t1 .* ((H.^2./8) + (s.^2./2) - (s.*H./2));

%Local s parameter as a function of y, z
s1 = @(z) -z;
s2 = @(y) (B/2) - y;
s3 = @(z) (H/2) + z;
s4 = @(y) (B/2) + y;
s5 = @(z) (H/2) - z;

%Geometry
websPos = cell(1, 5);
websPos{1}.y = B/2 * ones(1, nPointsPerSection);
websPos{1}.z = linspace(0, -H/2, nPointsPerSection);

websPos{2}.y = linspace(B/2, -B/2, nPointsPerSection);
websPos{2}.z = -H/2 * ones(1, nPointsPerSection);

websPos{3}.y = -B/2 * ones(1, nPointsPerSection);
websPos{3}.z = linspace(-H/2, H/2, nPointsPerSection);

websPos{4}.y = linspace(-B/2, B/2, nPointsPerSection);
websPos{4}.z = H/2 * ones(1, nPointsPerSection);

websPos{5}.y = B/2 * ones(1, nPointsPerSection);
websPos{5}.z = linspace(H/2, 0, nPointsPerSection);

%Plotting
shearFlowPlottingFactor = 0.5;

if plotSettings.plotAnalytical

    figure('Units', 'normalized', 'Position', [0.15 0.1 0.7 0.75])
    set(gcf, 'Name', 'Shear flow distribution')

    ax = gca;

    hold on
    for p=1:5

        plot(websPos{p}.y, websPos{p}.z, 'b');

    end

    %Segment 1
    plot_open = plot(websPos{1}.y - shearFlowPlottingFactor*(q1(s1(websPos{1}.z))), websPos{1}.z, 'r');

    %Segment 2
    plot(websPos{2}.y, websPos{2}.z + shearFlowPlottingFactor*(q2(s2(websPos{2}.y))), 'r')

    %Segment 3
    plot(websPos{3}.y + shearFlowPlottingFactor*(q3(s3(websPos{3}.z))), websPos{3}.z, 'r')

    %Segment 4
    plot(websPos{4}.y, websPos{4}.z - shearFlowPlottingFactor*(q4(s4(websPos{4}.y))), 'r')

    %Segment 5
    plot(websPos{5}.y - shearFlowPlottingFactor*(q5(s5(websPos{5}.z))), websPos{5}.z, 'r')

    set(ax, 'Ydir', 'reverse')
    set(ax, 'Xdir', 'reverse')
    set(ax, 'XAxisLocation', 'origin')
    set(ax, 'YAxisLocation', 'origin')

    ylabel('z')
    xlabel('y')

    grid minor

    %Check shear flow continuity
    assert(round(abs(q1(s1(websPos{1}.z(end)))), 4) == round(abs(q2(s2(websPos{2}.y(1)))), 4), 'Error: Boundary condition for the shear flow error');

    assert(round(abs(q2(s2(websPos{2}.y(end)))), 4) == round(abs(q3(s3(websPos{3}.z(1)))), 4), 'Error: Boundary condition for the shear flow error');

end

%Shear center open section
y_sc_open = ((E1*t1)/Phi_y) * (-1) * ( (5/24*H^3*B) + (1/2*H^2*B^2) + ((E2*t2)/(E1*t1)*H^3*B/24) );

%y_sc_open = -54.099;
if plotSettings.plotAnalytical
    plot_y_sc_open = plot(y_sc_open,0,'+g');
end

%% Constant shear flow
term_A = ( (2*B+H)/(G1*t1) ) + (H/(G2*t2));
term_B = (E1/G1)*( (-H^3/6) - ((3/4)*H^2*B) - (0.5*H*B^2) -( (E2*t2/(E1*t1))*H^3/12 ) );

q0 = (Q_z/Phi_y)*(term_B/term_A);

%Shear center closed section
y_sc_closed = - y_sc_open + ( (2*q0*A)/Q_z );
if plotSettings.plotAnalytical
    plot_y_sc_closed = plot(y_sc_closed,0,'+r');
end

oper.y_sc_closed = y_sc_closed;

%% Total shear flow.

if plotSettings.plotAnalytical
    %Segment 1
    plot_closed = plot(websPos{1}.y - shearFlowPlottingFactor*((q1(s1(websPos{1}.z)) - q0)), websPos{1}.z, 'k');

    %Segment 2
    plot(websPos{2}.y, websPos{2}.z + shearFlowPlottingFactor*((q2(s2(websPos{2}.y)) - q0)), 'k')

    %Segment 3
    plot(websPos{3}.y + shearFlowPlottingFactor*((q3(s3(websPos{3}.z)) - q0)), websPos{3}.z, 'k')

    %Segment 4
    plot(websPos{4}.y, websPos{4}.z - shearFlowPlottingFactor*((q4(s4(websPos{4}.y)) - q0)), 'k')

    %Segment 5
    plot(websPos{5}.y - shearFlowPlottingFactor*((q5(s5(websPos{5}.z)) - q0)), websPos{5}.z, 'k')

end


%Shear flow due to shear force displacement from the shear center
M_t = Q_z * (y_load - y_sc_closed);

q_m = M_t / (2*A);

%Total shear flow
q_total =cell(1, 5);
for i=1:5

    q_total{i} = zeros(1, nPointsPerSection);

    switch i
        case 1
            q_total{1} = -q_m + (q1(s1(websPos{1}.z)) - q0);
            if plotSettings.plotAnalytical
                plot_total = plot(websPos{1}.y + shearFlowPlottingFactor*q_total{1}, websPos{1}.z, 'g');
            end

        case 2
            q_total{2} = -q_m + (q2(s2(websPos{2}.y)) - q0);
            if plotSettings.plotAnalytical
                plot(websPos{2}.y, websPos{2}.z - shearFlowPlottingFactor*q_total{2}, 'g')
            end

        case 3
            q_total{3} = -q_m + (q3(s3(websPos{3}.z)) - q0);
            if plotSettings.plotAnalytical
                plot(websPos{3}.y - shearFlowPlottingFactor*q_total{3}, websPos{3}.z, 'g')
            end

        case 4
            q_total{4} = -q_m + (q4(s4(websPos{4}.y)) - q0);
            if plotSettings.plotAnalytical
                plot(websPos{4}.y, websPos{4}.z + shearFlowPlottingFactor*q_total{4}, 'g')
            end

        case 5
            q_total{5} = -q_m + (q5(s5(websPos{5}.z)) - q0);
            if plotSettings.plotAnalytical
                plot(websPos{5}.y + shearFlowPlottingFactor*q_total{5}, websPos{5}.z, 'g')
            end
    end

end


if plotSettings.plotAnalytical
    legend([plot_open plot_closed plot_total plot_y_sc_open plot_y_sc_closed], 'Shear flow for open section', 'Shear flow for closed section', 'Total shear flow', 'Shear center for open section', 'Shear center for closed section', 'location','best');
end

%Get out total shear flow
oper.q_total = q_total;

%Torsional stiffness
oper.torStiff = (4*A^2) / term_A;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotSettings.plotAnalytical && plotSettings.plotDistributedLoad
	figure('Units', 'normalized', 'Position', [0.15 0.1 0.7 0.75])
	set(gcf, 'Name', 'Twist due to distributed load')
	ax_distributed = gca;
	xlabel('z/L')
	ylabel('\phi [deg]')
end

if plotSettings.plotAnalytical && ~ plotSettings.plotDistributedLoad
    figure('Units', 'normalized', 'Position', [0.15 0.1 0.7 0.75])
    set(gcf, 'Name', 'Twist due to distributed load')
    ax_concentrated = gca;
    xlabel('z/L')
    ylabel('\phi [deg]')
end

xAdimSec = linspace(0, 1, 100);

F_x_distributed = @(x) (loadCase.q_z .* geom.L) + (loadCase.q_z .* x);

M_t_distributed = @(x) F_x_distributed(x) .* (y_load - oper.y_sc_closed);

specific_twist_fun = @(x) (M_t_distributed(x) ./ oper.torStiff);

twist_fun = @(x) specific_twist_fun(x) .* x;

<<<<<<< HEAD
% twist_concentratedLoad = (((loadCase.q_z .* geom.L) .* (y_load - oper.y_sc_closed)) ./ oper.torStiff) .* geom.L;

% twist_concentratedLoad = ((-150000) ./ oper.torStiff) .* geom.L;
=======
twist_concentratedLoad = (((loadCase.q_z .* geom.L) .* (y_load - oper.y_sc_closed)) ./ oper.torStiff) .* geom.L;
>>>>>>> parent of ded8c33... Added code to make parameters plots. Also added constant moment.

twist_concentratedLoad = ((loadCase.Q_z_total .* (y_load - oper.y_sc_closed)) ./ oper.torStiff) .* xAdimSec .* geom.L;

if plotSettings.plotAnalytical && plotSettings.plotDistributedLoad
    plot(ax_distributed, xAdimSec .* geom.L, twist_fun(xAdimSec .* geom.L) .* (180/pi))
    FsClass.SetAxisProp(ax_distributed, plotSettings);
end

if plotSettings.plotAnalytical && ~ plotSettings.plotDistributedLoad
    plot(ax_concentrated, xAdimSec .* geom.L, twist_concentratedLoad .* (180/pi))
    FsClass.SetAxisProp(ax_concentrated, plotSettings);
end 