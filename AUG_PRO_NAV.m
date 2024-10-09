function AUG_PRO_NAV
    % Nonlinear Augmented Proportional Navigation (APNG) Engagement Simulator
    % This GUI allows the user to simulate a missile-target engagement using
    % nonlinear APNG guidance law with optional noise and command limits.
    %
    % ... Created by: Islam Elnady - islamelnady@yahoo.com
  
    % Create the main figure window
    set(0, 'DefaultFigureWindowStyle', 'default');
    set(0, 'DefaultTextInterpreter', 'tex');
    set(0, 'DefaultLegendInterpreter', 'tex');
    set(0, 'DefaultAxesTickLabelInterpreter', 'tex');

    hFig = figure('Name', 'Nonlinear APNG Engagement Simulator', 'NumberTitle', 'off', 'Position', [100, 100, 900, 800], ...
        'WindowState','maximized','Units','normalized');

    % Load the background image from the current working directory
    imageFileName = 'apng.png';  % Replace with your image filename
    imagePath = fullfile(pwd, imageFileName);  % Construct full path to the image
    
    bgImage = imread(imagePath);  % Load the image
    % Create an invisible axes for background
    hBackground = axes('Parent', hFig, 'Position', [0.05, 0.36, 0.2, 1], 'Visible', 'off');
    imshow(bgImage, 'Parent', hBackground);
    % Set the axes to be the bottom layer
    uistack(hBackground, 'bottom');


    %% Input Fields for Initial Conditions
    x1 = 10; y1 = 550; w1 = 150; w2 = 50; 
    uicontrol('Style', 'text', 'Position', [x1 y1, w1, 20], 'String', 'Target Initial Pos. [Rt] (m):');
    hRt = uicontrol('Style', 'edit', 'Position', [x1+130, y1, w2, 20], 'String', '5000; 5000');
    
    uicontrol('Style', 'text', 'Position', [x1, y1-30, w1, 20], 'String', 'Missile Initial Pos. [Rm] (m):');
    hRm = uicontrol('Style', 'edit', 'Position', [x1+130, y1-30, w2, 20], 'String', '0; 5000');
    
    uicontrol('Style', 'text', 'Position', [x1, y1-60, w1, 20], 'String', 'Target Velocity [VT] (m/s):');
    hVT = uicontrol('Style', 'edit', 'Position', [x1+130, y1-60, w2, 20], 'String', '150');
    
    uicontrol('Style', 'text', 'Position', [x1, y1-90, w1, 20], 'String', 'Missile Velocity [VM] (m/s):');
    hVM = uicontrol('Style', 'edit', 'Position', [x1+130,  y1-90, w2, 20], 'String', '300');
    
    uicontrol('Style', 'text', 'Position', [x1, y1-120, w1, 20], 'String', 'Missile Heading Angle [HE] (rad):');
    hHE = uicontrol('Style', 'edit', 'Position', [x1+130, y1-120, w2, 20], 'String', '0');
    
    uicontrol('Style', 'text', 'Position', [x1, y1-150, w1, 20], 'String', 'Target Acceleration [nT] (g):');
    hNT = uicontrol('Style', 'edit', 'Position', [x1+130, y1-150, w2, 20], 'String', '2'); % Example: 2g
    
    uicontrol('Style', 'text', 'Position', [x1, y1-180, w1, 20], 'String', 'Navigation Constant [N]:');
    hN = uicontrol('Style', 'edit', 'Position', [x1+130, y1-180, w2, 20], 'String', '3');

    %% Input Fields for Noise (Optional)
    uicontrol('Style', 'text', 'Position', [x1, y1-210, w1, 20], 'String', 'Noise Sigma LambdaDot:');
    hSigmaLambdaDot = uicontrol('Style', 'edit', 'Position', [x1+130, y1-210, w2, 20], 'String', '0');
    
    uicontrol('Style', 'text', 'Position', [x1, y1-240, w1, 20], 'String', 'Noise Sigma Vc:');
    hSigmaVc = uicontrol('Style', 'edit', 'Position', [x1+130, y1-240, w2, 20], 'String', '0');
    
    uicontrol('Style', 'text', 'Position', [x1, y1-270, w1, 20], 'String', 'Noise Sigma nT:');
    hSigmanT = uicontrol('Style', 'edit', 'Position', [x1+130, y1-270, w2, 20], 'String', '0');


    %% Input Fields for Guidance Command Limits (Optional)
    uicontrol('Style', 'text', 'Position', [x1, y1-300, w1, 20], 'String', 'Upper Limit for nc:');
    hUpperLimit = uicontrol('Style', 'edit', 'Position', [x1+130, y1-300, w2, 20], 'String', 'inf');
    
    uicontrol('Style', 'text', 'Position', [x1, y1-330, w1, 20], 'String', 'Lower Limit for nc:');
    hLowerLimit = uicontrol('Style', 'edit', 'Position', [x1+130,  y1-330, w2, 20], 'String', '-inf');

    %% Simulation and Reset Buttons
    uicontrol('Style', 'text','Units','normalized', 'Position', [0.037, 0.222, 0.093, 0.024], 'String', 'Simulation Time Step (s):', ...
              'BackgroundColor',[0.6 0.9 0.9]);
    hdt = uicontrol('Style', 'edit', 'Position', [x1+130, y1-370, w2, 18], 'String', '1e-3');

    uicontrol('Style', 'pushbutton', 'Position', [50, y1-410, w2, 25], 'String', 'Simulate', ...
             'Callback', @runSimulation,'BackgroundColor','g','FontSize',11,'FontWeight','bold');
    uicontrol('Style', 'pushbutton', 'Position', [x1+130, y1-410, w2, 25], 'String', 'Reset', ...
             'Callback', @resetFields,'BackgroundColor','r','FontSize',11,'FontWeight','bold');

    %% Label for Miss Distance
    hMissDistanceLabel = uicontrol('Style', 'text', 'Position', [x1+20, y1-455, w1+20, 30], ...
                                    'String', 'Final Miss Distance: ...', ...
                                    'FontWeight', 'bold', 'FontSize', 12, ...
                                    'ForegroundColor', 'red');

    %% Credit
    uicontrol('Style', 'text','Units','normalized', 'Position', [0.01, 0.013, 0.25, 0.025], ...
        'String', 'Created by: Islam Elnady * islamelnady@yahoo.com*', ...
             'BackgroundColor',[0.95 0.95 0.0],'FontSize',10,'HorizontalAlignment','center');

    % Normalize units for all objects
    set(hFig.Children, 'Units', 'normalized');
    %% Axes for plotting results
    hTrajPlot = axes('Parent', hFig, 'Position',      [0.35, 0.65, 0.63, 0.31]); box on % Engagement Trajectories
    hPlotNc = axes('Parent', hFig, 'Position',        [0.35, 0.37, 0.295, 0.2]); box on% nc plot
    hPlotLambda = axes('Parent', hFig, 'Position',    [0.69, 0.37, 0.295, 0.2]); box on % lambda plot
    hPlotLambdaDot = axes('Parent', hFig, 'Position', [0.35, 0.09, 0.295, 0.2]); box on % lambdaDot plot
    hPlotVc = axes('Parent', hFig, 'Position',        [0.69, 0.09, 0.295, 0.2]); box on % Vc plot

    % Create a text label for the message
    hMessageLabel = uicontrol('Style', 'text', 'Units', 'normalized', ...
        'Position', [0.76, 0.012, 0.2, 0.04], 'String', '', ...  % Initial empty string
        'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [0.95, 0.95, 0.95]);
    
    % Hold on button
    uicontrol('Style', 'pushbutton', 'Units', 'normalized', ... % Use normalized units
        'Position', [0.58, 0.012, 0.16, 0.04], 'String', 'Hold on plots to compare!', ...
        'Callback', @(src, event) hold_on_button_callback([hTrajPlot, hPlotNc, hPlotLambda, hPlotLambdaDot, hPlotVc], hMessageLabel), ...
        'BackgroundColor', [0.4, 0.9, 0.7], 'FontSize', 10, 'FontWeight', 'bold');

    %% Callback Function for Running the Simulation
    function runSimulation(~, ~)
        % Get inputs from GUI
        Rt = str2num(get(hRt, 'String'));
        Rm = str2num(get(hRm, 'String'));
        VT = str2double(get(hVT, 'String'));
        VM = str2double(get(hVM, 'String'));
        HE = deg2rad(str2double(get(hHE, 'String')));
        
        % Convert nT from g to m/s^2
        g = 9.81; % acceleration due to gravity in m/s^2
        nT_g = str2double(get(hNT, 'String'));
        nT = nT_g * g; % convert to m/s^2
        
        N = str2double(get(hN, 'String'));
        sigma_LambdaDot = str2double(get(hSigmaLambdaDot, 'String'));
        sigma_nT = str2double(get(hSigmanT, 'String'));
        sigma_Vc = str2double(get(hSigmaVc, 'String'));
        upperLimit = str2double(get(hUpperLimit, 'String'));
        lowerLimit = str2double(get(hLowerLimit, 'String'));

        % Initial Conditions and Parameters
        dt = str2double(get(hdt, 'String'));
        t_max = 1000;
        beta = deg2rad(0);
        Vt = [-VT * cos(beta); VT * sin(beta)];
        nc = 0;
        Xr = Rt(1) - Rm(1);
        Yr = Rt(2) - Rm(2);
        lambda = atan2(Yr, Xr);
        L = asin(VT * sin(beta + lambda) / VM);
        Vm = [VM * cos(lambda + L + HE); VM * sin(lambda + L + HE)];
        Am = [0; 0];
        missile_positions = [];
        lm_array = [];
        target_positions = [];
        nc_array = [];
        lambda_array = [];
        lambdaDot_array = [];
        Vc_array = [];
        time_array = [];
        
        t = 0;
        R = norm(Rt - Rm);
        Vc = 1000;
               
        % Simulation loop
        while R >= 0.5 && Vc >= 0 && t <= t_max
            Rr = Rt - Rm;
            R = norm(Rr);
            Vr = Vt - Vm;
            lambda = atan2(Rr(2), Rr(1));
            lambdaD = (Rr(1) * Vr(2) - Rr(2) * Vr(1)) / (R^2);

            Vc = -(Rr(1) * Vr(1) + Rr(2) * Vr(2)) / R;

            lambdaD_n = lambdaD + sigma_LambdaDot * randn();
            Vc_n = Vc + sigma_Vc * randn();
            nT_n = nT + sigma_nT * randn();
            % Nonlinear APNG Command with limits
            nc = N * Vc_n * lambdaD_n + N * nT_n / 2;

            nc = max(min(nc, upperLimit), lowerLimit); % Apply limits to nc
            Am = [-nc * sin(lambda); nc * cos(lambda)];
            
            % RK4 for missile
            [Rm, Vm] = rk4_missile(Rm, Vm, Am, dt);
            gamma =   atan2(Vm(2),Vm(1));
            Lm = gamma - lambda;
            VM = norm(Vm);

            % RK4 for target
            betaD = nT / norm(Vt);
            beta = beta + betaD * dt;
            At = [nT * sin(beta); nT * cos(beta)];
            [Rt, Vt] = rk4_target(Rt, Vt, At, dt);
            
            VT = norm(Vt);
            L = asin(VT * sin(beta + lambda) / VM);

            % Store data
            missile_positions = [missile_positions, Rm];
            lm_array = [lm_array Lm];
            target_positions = [target_positions, Rt];
            nc_array = [nc_array, nc];
            lambda_array = [lambda_array, lambda];
            lambdaDot_array = [lambdaDot_array, lambdaD];
            Vc_array = [Vc_array, Vc];
            time_array = [time_array, t];
            
            % Time update
            t = t + dt;
        end
        
        % Calculate final miss distance
        final_miss_distance = norm(Rt - Rm);
       set(hMissDistanceLabel, 'String', ['Final Miss Distance: ', sprintf('%.3f', final_miss_distance), ' m']);
        
        % Check if missile intercepted the target
        if R < 1
            disp(['Missile intercepted the target at time t = ', num2str(t), ' seconds.']);
        else
            disp('Missile failed to intercept the target.');
        end
        
        % Plot the results
        plotTrajectories(missile_positions, target_positions, hTrajPlot);
        plotGuidanceCommands(time_array, nc_array./g, rad2deg(lambda_array),rad2deg(lm_array), ...
                             rad2deg(lambdaDot_array), Vc_array, hPlotNc, hPlotLambda, hPlotLambdaDot, hPlotVc);
    end

    %% RK4 Integration for Missile and Target
    function [R_next, V_next] = rk4_missile(R, V, A, dt)
        k1_V = A * dt;
        k1_R = V * dt;

        k2_V = (A + k1_V / 2) * dt;
        k2_R = (V + k1_R / 2) * dt;

        k3_V = (A + k2_V / 2) * dt;
        k3_R = (V + k2_R / 2) * dt;

        k4_V = (A + k3_V) * dt;
        k4_R = (V + k3_R) * dt;

        V_next = V + (k1_V + 2*k2_V + 2*k3_V + k4_V) / 6;
        R_next = R + (k1_R + 2*k2_R + 2*k3_R + k4_R) / 6;
    end

    function [R_next, V_next] = rk4_target(R, V, A, dt)
        k1_V = A * dt;
        k1_R = V * dt;

        k2_V = (A + k1_V / 2) * dt;
        k2_R = (V + k1_R / 2) * dt;

        k3_V = (A + k2_V / 2) * dt;
        k3_R = (V + k2_R / 2) * dt;

        k4_V = (A + k3_V) * dt;
        k4_R = (V + k3_R) * dt;

        V_next = V + (k1_V + 2*k2_V + 2*k3_V + k4_V) / 6;
        R_next = R + (k1_R + 2*k2_R + 2*k3_R + k4_R) / 6;
    end

    %% Plot Functions
    fontSize = 12;
    held = 0;
    runCount = 1;
    legtraj = {};  legnc = {};   legl = {};  legld = {};   legvc = {};

    function plotTrajectories(missile_positions, target_positions, ax)
    axes(ax);
    plot(missile_positions(1,:), missile_positions(2,:), 'LineWidth',2);
    hold on;
    plot(target_positions(1,:), target_positions(2,:), 'LineWidth', 2);
    xlabel('X Position (m)', 'FontSize', fontSize);
    ylabel('Y Position (m)', 'FontSize', fontSize);
    title('Engagement Trajectories', 'FontSize', fontSize);
    hold off;
    legtraj{end+1} = [' Missile - Run ', num2str(runCount)];
    legtraj{end+1} = [' Target - Run ', num2str(runCount)];

    if held == 1  % Update the legend with all runs
        legend(ax, legtraj);
    else
        % Show only the default legend for the first run
        legend(' Missile', ' Target');
    end 
    grid on;    axis equal;   box on;  % Add box around the plot
    end


function plotGuidanceCommands(time_array, nc_array, lambda_array, lm_array, lambdaDot_array, Vc_array, ...
                              axNc, axLambda, axLambdaDot, axVc)
    
    % Plot nc
    axes(axNc);
    plot(time_array, nc_array, 'LineWidth', 2);
    xlabel('Time (s)', 'FontSize', fontSize);
    ylabel('Guidance Command, n_c  (g)', 'FontSize', fontSize);
    title('Guidance Command vs Time', 'FontSize', fontSize);
    legnc{end+1} = [' n_c - Run ', num2str(runCount)];

    if held == 1  % Update the legend with all runs
        legend(axNc, legnc);
    else
        legend off;
    end
    grid on; box on;
    
    
    % Plot lambda and leading angle
    axes(axLambda);
    plot(time_array, lambda_array, 'LineWidth', 2, 'DisplayName', [' \lambda - Run ', num2str(runCount)]);
    hold on;
    plot(time_array, lm_array, 'LineWidth', 2, 'DisplayName', ['\beta_m - Run ', num2str(runCount)]);
    xlabel('Time (s)', 'FontSize', fontSize);
    ylabel('\lambda_L_O_S, \beta_m (^o)', 'FontSize', fontSize);
    title('LOS Angle and Leading Angle vs Time', 'FontSize', fontSize);
    hold off;

    % Append the new legend entries
      legl{end+1} = [' \lambda - Run ', num2str(runCount)];
      legl{end+1} = [' \beta_m - Run ', num2str(runCount)];
    
    if held == 1  % Update the legend with all runs
        legend(axLambda, legl);
    else
        % Show only the default legend for the first run
        legend('\lambda', '\beta_m');
    end 
    grid on; box on;


    % Plot lambdaDot
    axes(axLambdaDot);
    plot(time_array, lambdaDot_array, 'LineWidth', 2);
    xlabel('Time (s)', 'FontSize', fontSize);
    ylabel('LOS Angle Rate, $\dot{\lambda}$ ($^o$/s)', 'Interpreter', 'latex', 'FontSize', fontSize);
    title('Rate of Change of LOS Angle vs Time', 'FontSize', fontSize);
    
    legld{end+1} = [' \lambdadot - Run ', num2str(runCount)];
    if held == 1  % Update the legend with all runs
        legend(axLambdaDot, legld);
    else
        % Show only the default legend for the first run
        legend off;
    end  
    grid on; box on;

    % Plot Vc
    axes(axVc);
    plot(time_array, Vc_array, 'LineWidth', 2);
    xlabel('Time (s)', 'FontSize', fontSize);
    ylabel('Closing Velocity, Vc (m/s)', 'FontSize', fontSize);
    title('Closing Velocity vs Time', 'FontSize', fontSize);

    legvc{end+1} =  [' V_c Run ', num2str(runCount)];
    if held == 1  % Update the legend with all runs
        legend(axVc, legvc);
    else
        % Show only the default legend for the first run
        legend off;
    end    
    grid on; box on;

    runCount = runCount + 1 ; % Increment run count after each simulation
  
end

    % Callback function for the button
    function hold_on_button_callback(axHandles, messageLabel)
        held = 1;
        hold(axHandles, 'on');  % Hold on for each axis
        set(messageLabel, 'String', 'Press again after every run!');  % Update the message label
        disp('Holding on to the plots for comparison! Press again after every run!');
    end

  %% Reset All
    function resetFields(~, ~)
        set(hRt, 'String', '5000; 5000');
        set(hRm, 'String', '0; 5000');
        set(hVT, 'String', '150');
        set(hVM, 'String', '300');
        set(hHE, 'String', '0');
        set(hNT, 'String', '2'); % Reset to 2g
        set(hN, 'String', '3');
        set(hSigmaLambdaDot, 'String', '0');
        set(hSigmanT, 'String', '0');
        set(hSigmaVc, 'String', '0');
        set(hUpperLimit, 'String', 'inf');
        set(hLowerLimit, 'String', '-inf');
        set(hdt, 'String', '1e-3');
        set(hMissDistanceLabel, 'String', 'Final Miss Distance: ...');
        
        % Clear all plots
        set(hMessageLabel, 'String', '');  % Clear the message
        hold([hTrajPlot,hPlotNc,hPlotLambda,hPlotLambdaDot,hPlotVc],'off')
        cla(hTrajPlot);
        cla(hPlotNc);
        cla(hPlotLambda);
        cla(hPlotLambdaDot);
        cla(hPlotVc); 
        runCount = 1;
        held = 0;
        legtraj = {};  legnc = {};   legl = {};  legld = {};   legvc = {};
    end
end
