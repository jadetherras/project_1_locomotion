
%% gait parameters file

% In this file, we calculate some gait parameters

%% Loading the data

%dataset sain 
%Jade
data_healthy=load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/4_AML01_3kmh.mat");

% dataset SCI Human
%Jade
data =load("SCI Human/DM002_TDM_08_1kmh.mat");

%Lena
%data_healthy = load( 'DATA/1_AML01_2kmh.mat');
%data = load( 'Dataset SCI Human/SCI Human/DM002_TDM_08_1kmh.mat');

data = data.data;
data_healthy = data_healthy.data;







%% gate parameters for healthy subject

gates = cut_gate(data_healthy);
max_joint_angle = []; %utile
min_joint_angle = []; %utile
avg_angle_vel =[]; 
min_angle_vel =[]; 
max_angle_vel =[]; %utile

var_ver_hip = []; %utile
var_lateral_hip = []; %utile

for i = 1:(length(gates)-1)
    x1 = data_healthy.LKNE(gates(i):gates(i+1)-1,2);
    x3 = data_healthy.LHIP(gates(i):gates(i+1)-1,2);
    x2 = data_healthy.LANK(gates(i):gates(i+1)-1,2);
    
    %z-space
    y1 = data_healthy.LKNE(gates(i):gates(i+1)-1,3);
    y3 = data_healthy.LHIP(gates(i):gates(i+1)-1,3);
    y2 = data_healthy.LANK(gates(i):gates(i+1)-1,3);

    x_hip_lat = data_healthy.LHIP(gates(i):gates(i+1)-1,1);
    var_lateral_hip = [var_lateral_hip, var(x_hip_lat)];
    var_ver_hip = [var_ver_hip,var(y3)];
    
    %angle = atan2(vector2.y, vector2.x) - atan2(vector1.y, vector1.x);
    angle = atan2(y3 - y1, x3 - x1) - atan2(y2 - y1, x2 - x1);
    angle_rad = angle*(180/pi);
    vel_angle = gradient(angle_rad);
    max_joint_angle = [max_joint_angle, max(angle_rad)];
    min_joint_angle = [min_joint_angle, min(angle_rad)];
    %avg_angle_vel =[avg_angle_vel, mean(abs(vel_angle)) ];
    %min_angle_vel =[min_angle_vel, min(abs(vel_angle))]; 
    max_angle_vel =[max_angle_vel, max(abs(vel_angle))];
end

disp(max_joint_angle);


%% gate parameters for SCI patients

gates_SCI = cut_gate(data);
max_joint_angle_SCI = [];
min_joint_angle_SCI = [];
avg_angle_vel_SCI =[];
min_angle_vel_SCI =[];
max_angle_vel_SCI =[];

var_ver_hip_SCI = [];
var_lateral_hip_SCI = [];

for i = 1:(length(gates_SCI)-1)
    x1 = data.LKNE(gates_SCI(i):gates_SCI(i+1)-1,2);
    x3 = data.LHIP(gates_SCI(i):gates_SCI(i+1)-1,2);
    x2 = data.LANK(gates_SCI(i):gates_SCI(i+1)-1,2);
    
    %z-space
    y1 = data.LKNE(gates_SCI(i):gates_SCI(i+1)-1,3);
    y3 = data.LHIP(gates_SCI(i):gates_SCI(i+1)-1,3);
    y2 = data.LANK(gates_SCI(i):gates_SCI(i+1)-1,3);

    x_hip_lat = data.LHIP(gates_SCI(i):gates_SCI(i+1)-1,1);
    var_lateral_hip_SCI = [var_lateral_hip_SCI, var(x_hip_lat)];
    var_ver_hip_SCI = [var_ver_hip_SCI,var(y3)];
    
    %angle = atan2(vector2.y, vector2.x) - atan2(vector1.y, vector1.x);
    angle = atan2(y3 - y1, x3 - x1) - atan2(y2 - y1, x2 - x1);
    angle_rad = angle*(180/pi);
    vel_angle = gradient(angle_rad);
    max_joint_angle_SCI = [max_joint_angle_SCI, max(angle_rad)];
    min_joint_angle_SCI = [min_joint_angle_SCI, min(angle_rad)];
    
    max_angle_vel_SCI =[max_angle_vel_SCI, max(abs(vel_angle))];
end

disp(max_joint_angle_SCI);

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
    S = data.LTOE(:,2);
    S_L = filtering(S);

    % calculate the gradient
    G = gradient(S_L);

    %initialisation
    gate = [];

    % calculate the position
    for i = 2:length(G)
        if sign(G(i)) ~= sign(G(i-1)) && sign(G(i)) < 0 
            % if change sign, pass 0 and <0 foot off
                [~,index] = max(S(i-20:i+20));
                gate = [gate,index+20+i];
        end
    end
end


% cut_gate : cut the date in gate cycle (using a gradient)
% each value represent a foot strike
% data is the dataset
% gate is a array of position 
function gate = old_cut_gate(data) 

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