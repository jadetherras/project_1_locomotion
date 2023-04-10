
%% gait parameters file - PCA

% In this file, we calculate some gait parameters

%% Loading the data

%dataset sain 
data_healthy=load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/4_AML01_3kmh.mat");

% dataset SCI Human
data_SCI=load("SCI Human/DM002_TDM_08_1kmh.mat");

%% calculation of all parameters

% we can do all at one, just for optimisation

% test plot
plot(data_healthy.data.LTOE(27:191,2))
    hold on
    plot(data_healthy.data.RTOE(27:191,2))
    plot(filtering(data_healthy.data.LTOE(27:191,2)))
    plot(filtering(data_healthy.data.RTOE(27:191,2)))
    legend(["1","2","3","4"])
    hold off

% calculate parameters
parameters_healthy = calculated_parameters(data_healthy.data);
parameters_SCI = calculated_parameters(data_SCI.data);

%merge
parameters = [parameters_healthy,parameters_SCI];

%% PCA


coef = pca(parameters);



%% gate parameters functions

function parameters = calculated_parameters(data)

cut = cut_gate(data);

parameters = [];

    for i = 1:(length(cut)-1) 
        parameters = [parameters,struct2cell(gate_parameters(data, cut(i),cut(i+1)-1))];
    end

end

function parameters = gate_parameters(data, start, stop)
    
    % preli

    T = 1/120;
    speed = 3000/3600*1000;

    % jade 

    [off_distL, off_indexL] = max(data.LTOE(start:stop,2));
    [strike_distL,strike_indexL] = min(data.LTOE(start:stop,2));
    
    [off_distR, off_indexR] = max(data.RTOE(start:stop,2));
    [strike_distR,strike_indexR] = min(data.RTOE(start:stop,2));

    parameters.gate_duration_sec = (stop-start)*T;

    parameters.swing_duration_left_sec = (strike_indexL-off_indexL)*T;
    
    parameters.swing_duration_right_sec = (strike_indexR-off_indexR)*T;

    parameters.stance_percentage_left = (parameters.gate_duration_sec - parameters.swing_duration_left_sec)/parameters.gate_duration_sec;

    parameters.stance_percentage_right = (parameters.gate_duration_sec - parameters.swing_duration_right_sec)/parameters.gate_duration_sec;
    
    parameters.step_height_left_mm = max(data.LANK(start:stop,3))-min(data.LANK(start:stop,3));

    parameters.step_length_left_mm = abs(strike_distL -(off_distL+parameters.swing_duration_left_sec*speed));

    parameters.step_height_right_mm = max(data.RANK(start:stop,3))-min(data.RANK(start:stop,3));

    parameters.step_length_right_mm = abs(strike_distR -(off_distR +parameters.swing_duration_right_sec*speed));

    % lena

    x1L = data.LKNE(start:stop,2);
    x3L = data.LHIP(start:stop,2);
    x2L = data.LANK(start:stop,2);
    
    %z-space
    y1L = data.LKNE(start:stop,3);
    y3L = data.LHIP(start:stop,3);
    y2L = data.LANK(start:stop,3);

    x1R = data.RKNE(start:stop,2);
    x3R = data.RHIP(start:stop,2);
    x2R = data.RANK(start:stop,2);
    
    %z-space
    y1R = data.RKNE(start:stop,3);
    y3R = data.RHIP(start:stop,3);
    y2R = data.RANK(start:stop,3);


    x_hip_latL = data.LHIP(start:stop,1);
    x_hip_latR = data.RHIP(start:stop,1);

    angleL = atan2(y3L - y1L, x3L - x1L) - atan2(y2L - y1L, x2L - x1L);
    angle_radL = angleL*(180/pi);
    vel_angleL = gradient(angle_radL);

    angleR = atan2(y3R - y1R, x3R - x1R) - atan2(y2R - y1R, x2R - x1R);
    angle_radR = angleR*(180/pi);
    vel_angleR = gradient(angle_radR);


    parameters.var_lateral_hip_left_mm = var(x_hip_latL);
    parameters.var_ver_hip_left_mm = var(y3L);
    parameters.max_joint_angle_left_deg = max(angle_radL);
    parameters.min_joint_angle_left_deg = min(angle_radL);
    parameters.max_angle_vel_left_deg = max(abs(vel_angleL));

    parameters.var_lateral_hip_right_mm = var(x_hip_latR);
    parameters.var_ver_hip_right_mm = var(y3R);
    parameters.max_joint_angle_right_deg = max(angle_radR);
    parameters.min_joint_angle_right_deg = min(angle_radR);
    parameters.max_angle_vel_right_deg = max(abs(vel_angleR));

end



%% Basic functions

% filtering : take a signal (ex the y value for a marker) and return a
% clean the signal S -> S_f
function [S_f] = filtering(S)

    %low, high pass filter
    S_f = lowpass(highpass(S,1e-1,1e2),0.6,1e2, 'ImpulseResponse','iir');

end


% cut_gate : cut the date in gate cycle (using a gradient)
% each value represent a foot strike
% data is the dataset
% gate is a array of position 
function gate = cut_gate(data) 

    %filtering
    S_L = filtering(data.LTOE(:,2));

    % calculate the gradient
    G = gradient(S_L);

    %initialisation
    gate = [];

    % calculate the position
    for i = 2:length(G)
        if sign(G(i)) ~= sign(G(i-1)) && sign(G(i)) < 0 
            % if change sign, pass 0 and <0 foot off
                gate = [gate,i];
        end
    end
end