function AUG_PRO_NAV
    % Nonlinear Augmented Proportional Navigation (APNG) Engagement Simulator
    % This GUI allows the user to simulate a missile-target engagement using
    % nonlinear APNG guidance law with optional noise and command limits.
    %
    % ... Created by: Islam Elnady - islamelnady@yahoo.com
  
    %% Setup default figure properties
    set(0, 'DefaultFigureWindowStyle', 'default');
    set(0, 'DefaultTextInterpreter', 'tex');
    set(0, 'DefaultLegendInterpreter', 'tex');
    set(0, 'DefaultAxesTickLabelInterpreter', 'tex');

    % GLobal Variables
    global stopSimulation guid_type 
    guid_type = 'True';
    fontSize = 11;
    held = 0;
    runCount = 1;
    legtraj = {};  legl = {};  legld = {};   legvc = {}; dmiss = {}; lega = {};
    LatD = {}; dmissHold = {}; LatDHold = {};

    % Create the main figure window
    hFig = figure('Name', 'Nonlinear APNG Engagement Simulator', 'NumberTitle', 'off', 'Position', [100, 100, 900, 800], ...
        'WindowState','maximized','Units','normalized');

    % Load the background image from the current working directory
    imageFileName = 'apng.png';  % Replace with your image filename
    imagePath = fullfile(pwd, imageFileName);  % Construct full path to the image
    
    bgImage = imread(imagePath);  % Load the image
    % Create an invisible axes for background
    hBackground = axes('Parent', hFig, 'Position', [0.02, 0.36, 0.19, 1], 'Visible', 'off');
    imshow(bgImage, 'Parent', hBackground);
    % Set the axes to be the bottom layer
    uistack(hBackground, 'bottom');


    %% Input Fields for Initial Conditions
    x1 = 10; y1 = 560; w1 = 150; w2 = 50; 
    uicontrol('Style', 'text', 'Position', [x1 y1, w1, 20], 'String', 'Target Initial Pos. [Rt] (m):');
    hRt = uicontrol('Style', 'edit', 'Position', [x1+130, y1, w2, 20], 'String', '5000; 5000');
    
    uicontrol('Style', 'text', 'Position', [x1 y1-30, w1, 20], 'String', 'Target Initial Heading. [beta] (deg) :');
    hbt = uicontrol('Style', 'edit', 'Position', [x1+130, y1-30, w2, 20], 'String', '0');
  
    %%  Pseudo Pilot Dynamics controls
    hTauPanel = uipanel('Parent', hFig,'Units', 'normalized', 'Position', [0.22 0.60 0.09 0.08] ...
        ,'BackgroundColor',[0.8 0.8 0.9]);  % Position and size of the frame
  
    % Add the button inside the panel
    uicontrol('Parent', hTauPanel, 'Style', 'text', 'String', '1st Order Dynamics', ...
              'Units', 'normalized', 'Position', [0.05 0.65 0.90 0.30],'BackgroundColor',[0.8 0.8 0.9]);

    % Add a label for the tau input field inside the panel
    uicontrol('Parent', hTauPanel, 'Style', 'text', 'String', 'Enter tau (s):', ...
             'Units', 'normalized', 'Position', [0.08 0.2 0.60 0.30]);
    hTauInput = uicontrol('Parent', hTauPanel, 'Style', 'edit', 'String', '0', ...
        'Units', 'normalized', 'Position', [0.71 0.2 0.25 0.30]);
    uicontrol('Style', 'text', 'Position', [x1, y1-60, w1, 20], 'String', 'Missile Initial Pos. [Rm] (m):');

    %%
    
    hRm = uicontrol('Style', 'edit', 'Position', [x1+130, y1-60, w2, 20], 'String', '0; 5000');
    
    uicontrol('Style', 'text', 'Position', [x1, y1-90, w1, 20], 'String', 'Target Velocity [VT] (m/s):');
    hVT = uicontrol('Style', 'edit', 'Position', [x1+130, y1-90, w2, 20], 'String', '150');
    
    uicontrol('Style', 'text', 'Position', [x1, y1-120, w1, 20], 'String', 'Missile Velocity [VM] (m/s):');
    hVM = uicontrol('Style', 'edit', 'Position', [x1+130,  y1-120, w2, 20], 'String', '300');
    
    uicontrol('Style', 'text', 'Position', [x1, y1-150, w1, 20], 'String', 'Missile Heading Error [HE] (deg):');
    hHE = uicontrol('Style', 'edit', 'Position', [x1+130, y1-150, w2, 20], 'String', '0');
    
    uicontrol('Style', 'text', 'Position', [x1, y1-180, w1, 20], 'String', 'Target Acceleration [nT] (g):');
    hNT = uicontrol('Style', 'edit', 'Position', [x1+130, y1-180, w2, 20], 'String', '2'); % Example: 2g
    
    uicontrol('Style', 'text', 'Position', [x1, y1-210, w1, 20], 'String', 'Navigation Constant [N]:');
    hN = uicontrol('Style', 'edit', 'Position', [x1+130, y1-210, w2, 20], 'String', '3');

    %% Input Fields for Noise (Optional)
    uicontrol('Style', 'text', 'Position', [x1, y1-240, w1, 20], 'String', 'Noise Sigma LambdaDot:');
    hSigmaLambdaDot = uicontrol('Style', 'edit', 'Position', [x1+130, y1-240, w2, 20], 'String', '0');
    
    uicontrol('Style', 'text', 'Position', [x1, y1-270, w1, 20], 'String', 'Noise Sigma Vc:');
    hSigmaVc = uicontrol('Style', 'edit', 'Position', [x1+130, y1-270, w2, 20], 'String', '0');
    
    uicontrol('Style', 'text', 'Position', [x1, y1-300, w1, 20], 'String', 'Noise Sigma nT:');
    hSigmanT = uicontrol('Style', 'edit', 'Position', [x1+130, y1-300, w2, 20], 'String', '0');


    %% Input Fields for Guidance Command Limits (Optional)
    uicontrol('Style', 'text', 'Position', [x1, y1-330, w1, 20], 'String', 'Upper Limit for nc (g):');
    hUpperLimit = uicontrol('Style', 'edit', 'Position', [x1+130, y1-330, w2, 20], 'String', 'inf');
    
    uicontrol('Style', 'text', 'Position', [x1, y1-360, w1, 20], 'String', 'Lower Limit for nc (g):');
    hLowerLimit = uicontrol('Style', 'edit', 'Position', [x1+130,  y1-360, w2, 20], 'String', '-inf');

    uicontrol('Style', 'text','Units','normalized', 'Position', [0.037, 0.21, 0.093, 0.024], 'String', 'Simulation Time Step (s):', ...
              'BackgroundColor',[0.6 0.9 0.9]);
    hdt = uicontrol('Style', 'edit', 'Position', [x1+130, y1-390, w2, 18], 'String', '1e-3');

%% TRUE OR PURE PRO NAV
    % Create a label
    uicontrol('Style', 'text', 'Position', [20, y1-420, 100, 20], 'String', 'Choose guidance type:', ...
              'FontSize', 10,'ForegroundColor','b');
    
    % Create checkboxes with "True" checked by default
    trueBox = uicontrol('Style', 'checkbox', 'String', 'True', 'Position', [130, y1-420, 50, 20], 'Value', 1, ...
              'Callback', @(src, event) setGuidType('True'));
    pureBox = uicontrol('Style', 'checkbox', 'String', 'Pure', 'Position', [170, y1-420, 50, 20], 'Value', 0, ...
              'Callback', @(src, event) setGuidType('Pure'));
    
    % Nested function to update the guid_type variable and uncheck the other option
    function setGuidType(option)
        if strcmp(option, 'True')
            set(trueBox, 'Value', 1);  % Check 'True' box
            set(pureBox, 'Value', 0);  % Uncheck 'Pure' box
        elseif strcmp(option, 'Pure')
            set(pureBox, 'Value', 1);  % Check 'Pure' box
            set(trueBox, 'Value', 0);  % Uncheck 'True' box
        end
        guid_type = option;
        % Display selected option
        disp(['Selected guidance type: ', option]);
    end

    %% Simulation and Reset Controls

    uicontrol('Style', 'pushbutton', 'Position', [50, y1-460, w2, 30], 'String', 'Simulate', ...
             'Callback', @runSimulation,'BackgroundColor','g','FontSize',11,'FontWeight','bold');
    uicontrol('Style', 'pushbutton', 'Position', [x1+130, y1-460, w2, 30], 'String', 'Reset', ...
             'Callback', @resetFields,'BackgroundColor','r','FontSize',11,'FontWeight','bold');

%% HELP
    uicontrol('Style', 'pushbutton','Units','Normalized', 'Position', [0.96, 0.01, 0.03, 0.03], ...
              'String', 'Help', 'Callback', @showHelp, 'BackgroundColor', 'k', ...
              'FontSize', 9, 'FontWeight', 'bold','ForegroundColor','w');
    
    % Function for showing help information
    function showHelp(~, ~)
    helpMessage = sprintf(['%% Help for Nonlinear APNG Engagement Simulator\n\n' ...
                           'This GUI allows the user to simulate a missile-target engagement using nonlinear APNG guidance law with optional noise and command limits.\n\n' ...
                           '%% Instructions:\n' ...
                           '1. Input initial conditions for the target and missile in the provided fields.\n' ...
                           '2. Adjust optional parameters such as noise and guidance command limits as needed.\n' ...
                           '3. Select the type of guidance: True or Pure.\n' ...
                           '4. Click "Simulate" to run the simulation and view the results.\n' ...
                           '5. Use "Reset" to clear all fields and plots.\n' ...
                           '6. **Important:** Press "Hold on plots to compare!" after each run to keep the previous results visible.\n\n' ...
                           '%% List of Labels and Their Rules:\n' ...
                           '- Target Initial Pos. [Rt]: Initial position of the target in meters.\n' ...
                           '- Target Initial Heading. [beta]: Initial heading of the target in degrees.\n' ...
                           '- Missile Initial Pos. [Rm]: Initial position of the missile in meters.\n' ...
                           '- Target Velocity [VT]: Velocity of the target in meters per second.\n' ...
                           '- Missile Velocity [VM]: Velocity of the missile in meters per second.\n' ...
                           '- Noise Sigma LambdaDot: Standard deviation for noise in LambdaDot.\n' ...
                           '- Noise Sigma Vc: Standard deviation for noise in closing velocity (Vc).\n' ...
                           '- Noise Sigma nT: Standard deviation for noise in target acceleration.\n' ...
                           '- Upper/Lower Limit for nc: Guidance command limits in g.\n' ...
                           '- Simulation Time Step (s): Time step for the simulation.\n\n' ...
                           '%% Description of 1st Order Dynamics and Delay:\n' ...
                           'The 1st order dynamics is represented by a time constant (tau) that models the response delay of the missile guidance system.\n' ...
                           'This delay could represent the autopilot dynamics.\n' ...
                           'The guidance command (n_c) is calculated with a time delay, making it more realistic in real-world applications.\n\n' ...
                           '%% Main Features:\n' ...
                           '- Visualization of engagement trajectories and key parameters.\n' ...
                           '- Real-time simulation with feedback.\n' ...
                           '- Option to hold on plots for comparison.\n' ...
                           '- Animation of missile and target engagement.']);
                       
    msgbox(helpMessage, 'Help', 'help');
    end



    %% Miss Distance
    hMissDistanceLabel = uicontrol('Style', 'text', 'Position', [x1+20, y1-495, w1+70, 25], ...
                                    'String', 'Final Miss Distance: ...', ...
                                    'FontWeight', 'bold', 'FontSize', 10,...
                                    'ForegroundColor', 'r','HorizontalAlignment','left');
    %% Lateral Divert
    hLateralD = uicontrol('Style', 'text', 'Position', [x1+20, y1-515, w1+70, 20], ...
                                    'String', 'Total Lateral Divert: ...', ...
                                    'FontWeight', 'bold', 'FontSize', 10, ...
                                    'ForegroundColor', 'b','HorizontalAlignment','left');
    %% Credit
    uicontrol('Style', 'text','Units','normalized', 'Position', [0.01, 0.013, 0.22, 0.025], ...
        'String', 'Created by: Islam Elnady * islamelnady@yahoo.com*', ...
             'BackgroundColor',[0.95 0.95 0.0],'FontSize',9,'HorizontalAlignment','center');


    set(hFig.Children, 'Units', 'normalized');
    
    %% Axes for plotting results
    hTrajPlot = axes('Parent', hFig,'Units','normalized','Position', [0.35,  0.66, 0.295, 0.3]); box on % Engagement Trajectories
    hAnimate = axes('Parent', hFig,'Units','normalized', 'Position', [0.69,  0.66, 0.295, 0.3]); box on % Animate Engagement Trajectories

    % Create a button to start the animation
    hButton = uicontrol('Style', 'pushbutton', 'String', 'Start Animation','Units','normalized', ...
              'Position', [0.8, 0.8, 0.083, 0.065],'BackgroundColor',[0.9 0.9 0.6],'FontSize',9, ...
              'Callback', @(src, event) startanimation(hAnimate));

    hPlotNc = axes('Parent', hFig, 'Position',        [0.35, 0.37, 0.295, 0.2]); box on% nc plot
    hPlotLambda = axes('Parent', hFig, 'Position',    [0.69, 0.37, 0.295, 0.2]); box on % lambda plot
    hPlotLambdaDot = axes('Parent', hFig, 'Position', [0.35, 0.09, 0.295, 0.2]); box on % lambdaDot plot
    hPlotVc = axes('Parent', hFig, 'Position',        [0.69, 0.09, 0.295, 0.2]); box on % Vc plot

    % Create a text label for the message
    hMessageLabel = uicontrol('Style', 'text', 'Units', 'normalized', ...
        'Position', [0.71, 0.012, 0.2, 0.03], 'String', '', ...  % Initial empty string
        'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [0.95, 0.95, 0.95],'ForegroundColor', 'red');
    
    % Hold on button
    uicontrol('Style', 'pushbutton', 'Units', 'normalized', ... % Use normalized units
        'Position', [0.58, 0.012, 0.16, 0.04], 'String', 'Hold on plots to compare!', ...
        'Callback', @(src, event) hold_on_button_callback([hTrajPlot,hAnimate, hPlotNc, hPlotLambda, hPlotLambdaDot, hPlotVc], hMessageLabel), ...
        'BackgroundColor', [0.4, 0.9, 0.7], 'FontSize', 10, 'FontWeight', 'bold');

    %% Callback Function for Running the Simulation
    function runSimulation(~, ~)
        set(hButton, 'Visible', 'on');
        if held ~= 1 
            cla(hAnimate);  
        end

        % Get inputs from GUI
        Rt = str2num(get(hRt, 'String'));
        Rm = str2num(get(hRm, 'String'));
        VT = str2double(get(hVT, 'String'));
        VM = str2double(get(hVM, 'String'));
        HE = deg2rad(str2double(get(hHE, 'String')));
        dt = str2double(get(hdt, 'String'));
        tau = str2double(get(hTauInput, 'String'));  % Get the tau value from the input
        if tau == 0; disp('Zero-Lag Guidance!')
        else; disp(['1st Order Lag Gudiance with time constant: ' num2str(tau) ' sec'])
        end
        
        alpha = dt / (tau + dt); % Calculate the alpha coefficient for the delay
        nc_d = 0;

        t_max = 1000;
        g = 9.81; % acceleration due to gravity in m/s^2
        nT_g = str2double(get(hNT, 'String'));
        nT = nT_g * g; % convert to m/s^2
        beta = deg2rad(str2double(get(hbt, 'String'))); % get target's heading in rads

        N = str2double(get(hN, 'String'));
        sigma_LambdaDot = str2double(get(hSigmaLambdaDot, 'String'));
        sigma_nT = str2double(get(hSigmanT, 'String'));
        sigma_Vc = str2double(get(hSigmaVc, 'String'));
        upperLimit = str2double(get(hUpperLimit, 'String')) * g;
        lowerLimit = str2double(get(hLowerLimit, 'String')) * g;

        % Initial Conditions and Parameters

        Vt = [-VT * cos(beta); VT * sin(beta)];
        nc = 0;
        Xr = Rt(1) - Rm(1);
        Yr = Rt(2) - Rm(2);
        lambda = atan2(Yr, Xr);
        L = asin(VT * sin(beta + lambda) / VM);
        Vm = [VM * cos(lambda + L + HE); VM * sin(lambda + L + HE)];
        Am = [0; 0];
        gamma = atan2(Vm(2),Vm(1));

        missile_pos = [];
        missile_vel = [];
        missile_acc = [];
        gamma_array = [];
        lm_array = [];
        
        target_pos = [];
        Rr_array = [];

        nc_array = [];
        AM_array = [];
        lambda_array = [];
        lambdaDot_array = [];
        Vc_array = [];
        time_array = [];
        
        t = 0;
        R = norm(Rt - Rm);
        Vc = 1000;
               
        % Simulation loop
        while R >= 0.1 && Vc >= 0 && t <= t_max
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
            if strcmp(guid_type, 'True')
               nc = N * Vc_n * lambdaD_n + N * nT_n / 2;
               nc = max(min(nc, upperLimit), lowerLimit); % Apply limits to nc
            % Apply first-order time delay if exists
               nc_d = alpha * nc + (1 - alpha) * nc_d; % Update delayed value
               AM = nc_d;
               Am = [-nc_d * sin(lambda); nc_d * cos(lambda)];

            elseif strcmp(guid_type, 'Pure')
               nc = N * VM * lambdaD_n + N * nT_n / 2;
               nc = max(min(nc, upperLimit), lowerLimit); % Apply limits to nc
            % Apply first-order time delay if exists
               nc_d = alpha * nc + (1 - alpha) * nc_d; % Update delayed value
               AM = nc_d;
               Am = [-nc_d * sin(gamma); nc_d * cos(gamma)];
            else
               disp('Invalid Guidance Type!');
            end
           
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
            missile_pos = [missile_pos, Rm];
            missile_vel = [missile_vel, Vm];
            missile_acc = [missile_acc, Am];
            gamma_array = [gamma_array gamma];
            lm_array = [lm_array Lm];

            target_pos = [target_pos, Rt];
            
            Rr_array = [Rr_array, Rr];


            nc_array = [nc_array, nc];
            AM_array = [AM_array AM];
            lambda_array = [lambda_array, lambda];
            lambdaDot_array = [lambdaDot_array, lambdaD];
            Vc_array = [Vc_array, Vc];
            time_array = [time_array, t];
            
            % Time update
            t = t + dt;
        end
                    
            % Calculate final miss distance
            md = norm(Rt - Rm);
            dmiss{end+1} = sprintf('%.2f m', md);  

            if held ~= 1
                set(hMissDistanceLabel, 'String', ['Final Miss Distance:  ', sprintf('%.3f', md), ' m']);
            elseif held == 1
                % Hold is active
                if isempty(dmissHold) && numel(dmiss) > 1  % Use a separate variable to store data when hold is active
                    % Initialize and start with the current and previous value
                    dmissHold = dmiss(end-1:end);  
                else
                    % Append new values as additional runs are done while hold is active
                    dmissHold{end+1} = sprintf('%.2f m', md);
                end
                % Concatenate all miss distances after hold
                allmd = strjoin(dmissHold, ', '); 
                set(hMissDistanceLabel, 'Backgroundcolor', [0.9 0.9 0.9], 'String', ['Final Miss Distances: ', allmd]);
            end




        % Check if missile intercepted the target
        if R < 1
            disp(['Missile intercepted the target at time t = ', num2str(t), ' seconds.']);
        else
            disp('Missile failed to intercept the target.');
        end 
    
        stopSimulation = false;

        % Plot the results

        % Make the handle accessible for the callback
        guidata(hFig, struct('missile_pos', missile_pos, 'missile_vel', missile_vel, ...
                         'missile_acc', missile_acc, 'target_pos', target_pos));

        plotTrajectories(missile_pos, target_pos, hTrajPlot);

        if length(nc_array) > 100
           % Apply median filter to the last 100 indices to suppress spikes
           time_array(end-99:end)      = medfilt1(time_array(end-99:end), 5);
           nc_array(end-99:end)        = medfilt1(nc_array(end-99:end), 5);
           lambda_array(end-99:end)    = medfilt1(lambda_array(end-99:end), 5);
           lm_array(end-99:end)        = medfilt1(lm_array(end-99:end), 5);
           lambdaDot_array(end-99:end) = medfilt1(lambdaDot_array(end-99:end), 5);
           Vc_array(end-99:end)        = medfilt1(Vc_array(end-99:end), 5);
        end

        % Get lateral divert
        LD = trapz(time_array,abs(nc_array)./9.81);
        LatD{end+1} = sprintf('%.2f g', LD);  % Curly braces {} for cell assignment 
        if held ~= 1
           set(hLateralD, 'String', ['Total Lateral Divert:  ', sprintf('%.3f', LD), ' g']);      
        elseif held == 1
                % Hold is active
                if isempty(LatDHold) && numel(LatD) > 1 % Use a separate variable to store data when hold is active
                    % Initialize and start with the current and previous value
                    LatDHold = LatD(end-1:end);  
                else
                    % Append new values as additional runs are done while hold is active
                    LatDHold{end+1} = sprintf('%.2f g', LD);
                end
            allLD = strjoin(LatDHold, ', '); % Concatenate with a comma and space
            set(hLateralD,'Backgroundcolor',[0.9 0.9 0.9], 'String', ['Total Lateral Divert: ', allLD]);
        end


        plotGuidanceCommands(time_array, nc_array, AM_array, rad2deg(lambda_array), rad2deg(gamma_array), rad2deg(lm_array), ...
                             rad2deg(lambdaDot_array), Vc_array, hPlotNc, hPlotLambda, hPlotLambdaDot, hPlotVc);
                
        runCount = runCount + 1 ; % Increment run count after each simulation

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
    function startanimation(ax)
        data = guidata(gca); % Get the data stored in the figure
        missile_pos = data.missile_pos;
        missile_vel = data.missile_vel;
        missile_acc = data.missile_acc;
        target_pos = data.target_pos;
        set(hButton, 'Visible', 'off');
        animateTraj(missile_pos, missile_vel, missile_acc, target_pos, ax)

    end

    function plotTrajectories(missile_pos, target_pos, ax)
        axes(ax);
    
        plot(missile_pos(1,:), missile_pos(2,:), 'LineWidth',2);
        hold on;
        plot(target_pos(1,:), target_pos(2,:), 'LineWidth', 2);
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

%% Animate Trajectories
   function animateTraj(missile_pos, missile_vel, missile_acc, target_pos, ax)
        axes(ax);  axis equal;  box on; hold on
        % Add box around the plot
        xlabel(ax,'X Position (m)', 'FontSize', fontSize);
        ylabel(ax,'Y Position (m)', 'FontSize', fontSize);
        title(ax,'Engagement Trajectories', 'FontSize', fontSize);
        htgt = plot(target_pos(:,1),  'r.', 'MarkerSize', 20); % Initialize red dot
        hmis = plot(missile_pos(:,1), 'b.', 'MarkerSize', 20); % Initialize red dot

        d =  round(length(missile_pos)/(150)); 
        for i = 1:d:length(missile_pos)  
     
        % Define the zoom range around the current position
            zoom_range = 1000; % You can adjust this range as needed
            % Check if the position index is within valid bounds
            if norm(target_pos(:,i) - missile_pos(:,i)) <= 1500 && i <= length(missile_pos)
                % Calculate the limits for X and Y
                xl = [missile_pos(1, i) - zoom_range, missile_pos(1, i) + zoom_range];
                yl = [missile_pos(2, i) - zoom_range, missile_pos(2, i) + zoom_range];
                
                % Update the plot limits
                set(gca, 'XLim', xl);
                set(gca, 'YLim', yl);
            end
    
        if stopSimulation
        disp('Simulation stopped.');
        break; % Exit the loop
        end
    
        plot(missile_pos(1,i), missile_pos(2,i),'ko', 'LineWidth', 1,'MarkerSize', 1); 
        set(hmis, 'XData', missile_pos(1,i), 'YData', missile_pos(2,i));

        plot(target_pos(1,i), target_pos(2,i),'ko', 'LineWidth', 1, 'MarkerSize', 1);
        set(htgt, 'XData', target_pos(1,i), 'YData', target_pos(2,i));

        if i > 1 % Ensure quivers are only drawn after the first iteration
           delete(qv_v);
           delete(qv_a);
           delete(hr);
        end
    
        hr =  plot([missile_pos(1,i) target_pos(1,i)], [missile_pos(2,i) target_pos(2,i)],'k--','LineWidth', 1);
    
        % Plot velocity vector (Vm)
        qv_v = quiver(missile_pos(1,i), missile_pos(2,i), missile_vel(1,i), missile_vel(2,i),0, ...
                     'g','LineWidth', 2,'MaxHeadSize',10);    
        % Plot guidance command vector (nc)
        qv_a =  quiver(missile_pos(1,i), missile_pos(2,i), missile_acc(1,i)*8, missile_acc(2,i)*8,0, ...
                      'r', 'LineWidth', 2,'MaxHeadSize',10);    
        legend([qv_v, qv_a], {' V_m', ' A_m'});
        drawnow
        pause(0.001);
        grid on;
        end 
        xlim('auto');  ylim('auto'); 
   end


   function plotGuidanceCommands(time_array, nc_array, AM_array, lambda_array, gamma_array, lm_array, lambdaDot_array, ...
                                 Vc_array,axNc, axLambda, axLambdaDot, axVc)
    
     % Plot nc
        axes(axNc);
        plot(time_array, nc_array./9.81, 'LineWidth', 2);
        hold on;
        plot(time_array, AM_array./9.81, 'LineWidth', 2);
        ylabel('n_c,  A_m (g)', 'FontSize', fontSize,'FontWeight','bold');
        title('Guidance Command vs Achieved Acceleration', 'FontSize', fontSize); 
        hold off;
        lega{end+1} = [' n_c - Run ', num2str(runCount)];
        lega{end+1} = [' A_m - Run ', num2str(runCount)];  % Update the legend with all runs
        if held == 1 
           legend(axNc, lega);
        else
           legend('n_c', 'A_m');
        end 
        grid on; box on;
        
        
     % Plot lambda and gamma
        axes(axLambda);
        plot(time_array, lambda_array, 'LineWidth', 2, 'DisplayName', [' \lambda - Run ', num2str(runCount)]);
        hold on;
        plot(time_array, gamma_array, 'LineWidth', 2, 'DisplayName', ['\beta_m - Run ', num2str(runCount)]);
        xlabel('Time (s)', 'FontSize', fontSize);
        ylabel('\lambda_{los},  \gamma_m (^o)', 'FontSize', fontSize,'FontWeight','bold');
        title('LOS Angle and Flight Path Angle vs Time', 'FontSize', fontSize);
        hold off;
        % Append the new legend entries
          legl{end+1} = [' \lambda - Run ', num2str(runCount)];
          legl{end+1} = [' \gamma_m - Run ', num2str(runCount)];
        
        if held == 1  % Update the legend with all runs
            legend(axLambda, legl);
        else
            % Show only the default legend for the first run
            legend('\lambda', '\gamma_m');
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
        stopSimulation = true;
        set(hRt, 'String', '5000; 5000');
        set(hbt, 'String', '0');
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
        set(hMissDistanceLabel,'BackgroundColor',[0.95 0.95 0.95], 'String', 'Final Miss Distance: ...');
        set(hLateralD,'BackgroundColor',[0.95 0.95 0.95], 'String', 'Total Lateral Divert: ...');
        set(hTauInput, 'String', '0');
        set(trueBox, 'Value', 1);  % Check 'True' box
        set(pureBox, 'Value', 0);  % Check 'True' box
        guid_type = 'True';

        
        % Clear all plots
        set(hMessageLabel, 'String', '');  % Clear the message
        % Clear all plots
        plotHandles = [hTrajPlot, hPlotNc, hPlotLambda, hPlotLambdaDot, hPlotVc, hAnimate];   
        for h = plotHandles
            cla(h);   title(h, '');  xlabel(h, '');   ylabel(h, '');   hold(h,'off')
        end
        runCount = 1;
        held = 0;
        legtraj = {};  legl = {};     legld = {};      legvc = {}; lega = {};
        dmiss   = {};  LatD = {}; dmissHold = {};   LatDHold = {};
    end 
end
