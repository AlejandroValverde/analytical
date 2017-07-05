classdef FsClass
    %% Auxiliary collection of functions
    %
    %   Original code by Alejandro Valverde: 19 May. 2017
    %   Last modified by Alejandro Valverde: 19 May. 2017
    %

    properties
       %
    end
    
    methods (Static)

        function [] = writeToFile(filename, dirWork, parameters, values, Iter)

            if isunix
                fileID = fopen([dirWork.abaqus '/' filename], 'w');
            elseif ispc
                fileID = fopen([dirWork.abaqus '\' filename], 'w');
            end

            fprintf(fileID, 'Iter\n');

            fprintf(fileID, [num2str(Iter) '\n']);

            for t = 1:length(parameters)

                fprintf(fileID, [parameters{t} '\n']);
                fprintf(fileID, [num2str(values{t}) '\n'] );

            end

            fclose(fileID);
        end

        function [dirWork] = organizeFolders(folderPostproc)

            if nargin<1
                folderPostproc = '.';
            end

            if exist('figures', 'dir') == 0
                mkdir('figures'); 
            end

            if isunix
                dirWork.figures = ['./figures/'];
            elseif ispc
                dirWork.figures = ['.\figures\'];
            end


            dirWork.main = pwd;

            if isunix
                dirWork.abaqus = ['../../workfolder'];
            elseif ispc
                dirWork.abaqus = ['..\..\workfolder'];
            end

            if isunix
                dirWork.postProc = ['../../workfolder/postProc/' folderPostproc];
            elseif ispc
                dirWork.postProc = ['..\..\workfolder\postProc\' folderPostproc];
            end

        end

        function [] = SetAxisProp(axesHandle, plotSet)

            set(axesHandle, 'GridAlpha', plotSet.axGridAlpha);
            set(axesHandle, 'FontSize', plotSet.axFontSize);
            set(axesHandle, 'LineWidth', plotSet.axLineWidth);
            set(axesHandle, 'TitleFontSizeMultiplier', plotSet.TitleFontSizeMultiplier);
            grid minor

        end

        function [twist] = finalGuessOfOverallTwistAtTip(xFit, yFit)

            %% Fit: 'untitled fit 1'.
            [xData, yData] = prepareCurveData( xFit, yFit );

            % Set up fittype and options.
            ft = fittype( 'poly1' );

            % Fit model to data.
            %Type of polyominal fitting algorithm, 'poly1': f(x) = p1*x + p2
            [fitresult, ~] = fit( xData, yData, ft );

            confi = confint(fitresult);

            %Build output variables
            p1 = [fitresult.p1, confi(1, 1), confi(2, 1)];
            p2 = [fitresult.p2, confi(1, 2), confi(2, 2)];

            %Twist at tip, x=1
            twist = p1(1) + p2(2);

        end

        function [oper] = shearFlowCalculations(geom, mat, loadCase, plotSettings)

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
        end

        function [xPos, dataU2, dataUR3] = loadAbaqus(dirWork, IterStr)

            %Move to the result folder
            cd(dirWork.postProc)

            if isunix
                abaqusFileInfo = dir([['.'], '/XYData_u2*']);
            elseif ispc
                abaqusFileInfo = dir([['.'], '\XYData_u2*']);
            end

            %Count number of files
            nFiles = length(abaqusFileInfo(not([abaqusFileInfo.isdir])));

            dataU2 = cell(1, nFiles);

            xPos = zeros(1, nFiles);

            for iFile=1:(nFiles)
                [u2, xPos_temp] = FsClass.importU2dataForXpos(abaqusFileInfo(iFile).name, dirWork, IterStr);

                dataU2{iFile} = u2;

                xPos(iFile) = xPos_temp;
            end

            if isunix
                abaqusFileInfo = dir([['.'], '/XYData_ur3*']);
            elseif ispc
                abaqusFileInfo = dir([['.'], '\XYData_ur3*']);
            end

            %Count number of files
            nFiles = length(abaqusFileInfo(not([abaqusFileInfo.isdir])));

            assert(nFiles == 4, 'There is more than four file for UR3')

            dataUR3 = [];

            for iFile=1:(nFiles)
                [dataUR3] = FsClass.importUR3dataForXpos(abaqusFileInfo(iFile).name, dirWork, IterStr, dataUR3);
            end

            %Return to original folder
            cd(dirWork.main)
        end

        function [u2, xPos] = importU2dataForXpos(filename, dirWork, IterStr)

            %% Initialize variables.
            delimiter = ' ';
            startRow = 3;
            endRow = 33;

            %% Create full path to file - Not neccessary  
            if isunix
                fullfilename = [dirWork.postProc '/' IterStr '/' filename];
            elseif ispc
                fullfilename = [dirWork.postProc '\' IterStr '\' filename];
            end

            %Obtain position
            filename1 = strrep(filename, 'XYData_u2_','');
            filename2 = strrep(filename1, '.rpt','');
            xPos = str2num(filename2);

            %% Format for each line of text:
            %   column2: double (%f)
            % For more information, see the TEXTSCAN documentation.
            formatSpec = '%*s%f%[^\n\r]';

            %% Open the text file.
            fileID = fopen(filename,'r');

            %% Read columns of data according to the format.
            % This call is based on the structure of the file used to generate this
            % code. If an error occurs for a different file, try regenerating the code
            % from the Import Tool.
            textscan(fileID, '%[^\n\r]', startRow(1)-1, 'WhiteSpace', '', 'ReturnOnError', false);
            dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            for block=2:length(startRow)
                frewind(fileID);
                textscan(fileID, '%[^\n\r]', startRow(block)-1, 'WhiteSpace', '', 'ReturnOnError', false);
                dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false, 'EndOfLine', '\r\n');
                dataArray{1} = [dataArray{1};dataArrayBlock{1}];
            end

            %% Close the text file.
            fclose(fileID);

            %% Allocate imported array to column variable names
            u2 = dataArray{:, 1};

        end

        function [dataUR3] = importUR3dataForXpos(filename, dirWork, IterStr, previousTable)

            %% Initialize variables.
            delimiter = ' ';
            startRow = 2;
            endRow = 32;

            %% Create full path to file - Not necessary 
            if isunix
                fullfilename = [dirWork.postProc '/' IterStr '/' filename];
            elseif ispc
                fullfilename = [dirWork.postProc '\' IterStr '\' filename];
            end

            %Obtain location
            filename1 = strrep(filename, 'XYData_ur3_','');
            filename2 = strrep(filename1, '.rpt','');
            loc = filename2;

            %% Format for each line of text:
            %   column2: double (%f)
            % For more information, see the TEXTSCAN documentation.
            formatSpec = '%f%f%[^\n\r]';

            %% Open the text file.
            fileID = fopen(filename,'r');

            %% Read columns of data according to the format.
            % This call is based on the structure of the file used to generate this
            % code. If an error occurs for a different file, try regenerating the code
            % from the Import Tool.
            textscan(fileID, '%[^\n\r]', startRow(1)-1, 'WhiteSpace', '', 'ReturnOnError', false);
            dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            for block=2:length(startRow)
                frewind(fileID);
                textscan(fileID, '%[^\n\r]', startRow(block)-1, 'WhiteSpace', '', 'ReturnOnError', false);
                dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false, 'EndOfLine', '\r\n');
                for col=1:length(dataArray)
                    dataArray{col} = [dataArray{col};dataArrayBlock{col}];
                end
            end

            %% Close the text file.
            fclose(fileID);

            %% Allocate imported array to column variable names
            if isempty(previousTable)
                dataUR3 = table(dataArray{1:end-1}, 'VariableNames', {'zAdim',['ur3_' loc]});
            else
                data_temp = table(dataArray{1:end-1}, 'VariableNames', {'zAdim','ur3'});
                previousTable.(['ur3_' loc]) = data_temp.ur3;
                dataUR3 = previousTable;
            end

        end

    end

end