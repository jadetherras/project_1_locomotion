
%% gait parameters file - PCA

% In this file, we calculate some gait parameters
% We did the PCA and biplot
% we explore the result
clear;
close all;

%% Loading the data jade

%dataset sain 
data_healthy_1_2kmh=load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/1_AML01_2kmh.mat");
data_healthy_2_2kmh=load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/1_AML02_2kmh.mat");

data_healthy_1_1kmh=load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/3_AML01_1kmh.mat");
data_healthy_2_1kmh=load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/3_AML02_1kmh.mat");

data_healthy_1_3kmh=load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/4_AML01_3kmh.mat");
data_healthy_2_3kmh=load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/4_AML02_3kmh.mat");

data_healthy_1_3kmh_inclined=load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/2_AML01_3kmh_inclined.mat");
data_healthy_2_3kmh_inclined=load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/2_AML02_3kmh_inclined.mat");

% dataset SCI Human
data_SCI_1kmh=load("SCI Human/DM002_TDM_08_1kmh.mat");
data_SCI_2kmh=load("SCI Human/DM002_TDM_08_2kmh.mat");

%dataset SCI Human noEES

data_SCI_1kmh_NoEES=load("SCI Human/DM002_TDM_1kmh_NoEES.mat");

% %% Loading the data Lena
% 
% %dataset sain 
% data_healthy_1_2kmh=load("DATA/1_AML01_2kmh.mat");
% data_healthy_2_2kmh=load("DATA/1_AML02_2kmh.mat");
% 
% data_healthy_1_1kmh=load("DATA/3_AML01_1kmh.mat");
% data_healthy_2_1kmh=load("DATA/3_AML02_1kmh.mat");
% 
% data_healthy_1_3kmh=load("DATA/4_AML01_3kmh.mat");
% data_healthy_2_3kmh=load("DATA/4_AML02_3kmh.mat");
% 
% % dataset SCI Human
% data_SCI_1kmh=load("Dataset SCI Human/SCI Human/DM002_TDM_08_1kmh.mat");
% data_SCI_2kmh=load("Dataset SCI Human/SCI Human/DM002_TDM_08_2kmh.mat");
% 
% %dataset SCI Human noEES
% 
% data_SCI_1kmh_NoEES=load("Dataset SCI Human/SCI Human/DM002_TDM_1kmh_NoEES.mat");

% 
% % Loading the data Lucas
% 
% %dataset sain 
% data_healthy_1_2kmh=load("Healthy/1_AML01_2kmh.mat");
% data_healthy_2_2kmh=load("Healthy/1_AML02_2kmh.mat");
% 
% data_healthy_1_1kmh=load("Healthy/3_AML01_1kmh.mat");
% data_healthy_2_1kmh=load("Healthy/3_AML02_1kmh.mat");
% 
% data_healthy_1_3kmh=load("Healthy/4_AML01_3kmh.mat");
% data_healthy_2_3kmh=load("Healthy/4_AML02_3kmh.mat");
% 
% % dataset SCI Human
% data_SCI_1kmh=load("SCI Human/DM002_TDM_08_1kmh.mat");
% data_SCI_2kmh=load("SCI Human/DM002_TDM_08_2kmh.mat");
% 
% %dataset SCI Human noEES
% 
% data_SCI_1kmh_NoEES=load("SCI Human/DM002_TDM_1kmh_NoEES.mat");

%% calculation of all parameters

% we can do all at one, just for optimisation

% % test plot
%  plot(data_healthy.data.LTOE(27:191,2))
%      hold on
%      plot(data_healthy.data.RTOE(27:191,2))
%      plot(filtering(data_healthy.data.LTOE(27:191,2)))
%      plot(filtering(data_healthy.data.RTOE(27:191,2)))
%      legend(["1","2","3","4"])
%      hold off

% calculate parameters

 %Gate = cut_gate(data_healthy_1_3kmh_inclined.data);
 
%   n = size(data_healthy_1_2kmh.data.LTOE(:,2));
%   t = 12166:12348;
%   t = t/100;
%   plot(t,[data_healthy_1_2kmh.data.LTOE(12166:12348,2),filtering(data_healthy_1_2kmh.data.LTOE(12166:12348,2)),data_healthy_1_2kmh.data.RTOE(12166:12348,2),filtering(data_healthy_1_2kmh.data.RTOE(12166:12348,2))]);
%   legend('L','L f','R','R f');
% 
% figure
% t = 2056:2355;
% t = t/100;
% plot(t,[data_healthy_1_1kmh.data.LTOE(2056:2355,2),filtering(data_healthy_1_1kmh.data.LTOE(2056:2355,2)),data_healthy_1_1kmh.data.RTOE(2056:2355,2),filtering(data_healthy_1_1kmh.data.RTOE(2056:2355,2))]);
% legend('L','L f','R','R f');
% 
% gate_parameters(data_healthy_1_1kmh.data, 2056, 2355,1)
% 
% figure
% t = 211:494;
% t = t/100;
% plot(t,[data_healthy_1_1kmh.data.LTOE(211:494,2),filtering(data_healthy_1_1kmh.data.LTOE(211:494,2)),data_healthy_1_1kmh.data.RTOE(211:494,2),filtering(data_healthy_1_1kmh.data.RTOE(211:494,2))]);
% legend('L','L f','R','R f');

parameters_healthy_1_1kmh = calculated_parameters(data_healthy_1_1kmh.data,1);
parameters_healthy_2_1kmh = calculated_parameters(data_healthy_2_1kmh.data,1);
  
N1H = size(parameters_healthy_1_1kmh,2) + size(parameters_healthy_2_1kmh,2);

parameters_healthy_1_2kmh = calculated_parameters(data_healthy_1_2kmh.data,2);
parameters_healthy_2_2kmh = calculated_parameters(data_healthy_2_2kmh.data,2);
  
N2H = size(parameters_healthy_1_2kmh,2) + size(parameters_healthy_2_2kmh,2);

parameters_healthy_1_3kmh = calculated_parameters(data_healthy_1_3kmh.data,3);
parameters_healthy_2_3kmh = calculated_parameters(data_healthy_2_3kmh.data,3);

N3H = size(parameters_healthy_1_3kmh,2) + size(parameters_healthy_2_3kmh,2);

parameters_healthy_1_3kmh_inclined = calculated_parameters(data_healthy_1_3kmh_inclined.data,3);
parameters_healthy_2_3kmh_inclined = calculated_parameters(data_healthy_2_3kmh_inclined.data,3);

N3HInclin = size(parameters_healthy_1_3kmh_inclined,2) + size(parameters_healthy_2_3kmh_inclined,2);

parameters_SCI_1kmh = calculated_parameters(data_SCI_1kmh.data,1);
parameters_SCI_2kmh = calculated_parameters(data_SCI_2kmh.data,2);

N1SCI = size(parameters_SCI_1kmh,2);
N2SCI = size(parameters_SCI_2kmh,2);

parameters_SCI_1kmh_NoEES = calculated_parameters(data_SCI_1kmh_NoEES.data,1);
NSCInoESS = size(parameters_SCI_1kmh_NoEES,2);

NH = N1H + N2H + N3H + N3HInclin;
NSCI = N1SCI + NSCInoESS + N2SCI;

%merge
parameters = [parameters_healthy_1_1kmh,parameters_healthy_2_1kmh,parameters_healthy_1_2kmh, parameters_healthy_2_2kmh,parameters_healthy_1_3kmh,parameters_healthy_2_3kmh,parameters_healthy_1_3kmh_inclined,parameters_healthy_2_3kmh_inclined,parameters_SCI_1kmh,parameters_SCI_2kmh,parameters_SCI_1kmh_NoEES];

%1
%parameters = [parameters_healthy_1_1kmh,parameters_healthy_2_1kmh,parameters_SCI_1kmh,parameters_SCI_1kmh_NoEES];
 
%% processing for PCA and plot

%from structure to matrix for the PCA analysis
prePCA = [];
for i = 1:(length(parameters))
    prePCA = [prePCA,cell2mat(struct2cell(parameters(i)))];
end
prePCA = transpose(prePCA);

val_label = fieldnames(parameters(1));

PCA = normalize(prePCA);

%healthy vs SCI
% color = [transpose(zeros(NH,1)+1), transpose(zeros(NSCI,1)+2)];
% color = transpose(color);

%healthy vs SCI vs stimulation
%color = [transpose(zeros(NH,1)+1), transpose(zeros(N1SCI + N2SCI,1)+2), transpose(zeros(NSCInoESS,1)+3)];
%color = transpose(color);

%VS differents conditions
color = [transpose(zeros(N1H,1)+1), transpose(zeros(N2H,1)+2), transpose(zeros(N3H,1)+3),transpose(zeros(N3HInclin,1)+4), transpose(zeros(N1SCI,1)+5),transpose(zeros(N2SCI,1)+6),transpose(zeros(NSCInoESS,1)+7)];
color = transpose(color);
%1
%color = [transpose(zeros(N1H,1)+1), transpose(zeros(N1SCI,1)+4),transpose(zeros(NSCInoESS,1)+6)];
%color = transpose(color);


%% corellation matrix

matrix = cov(PCA);

figure

imagesc(matrix)
xticks = linspace(1, size(val_label,1),size(val_label,1));
set(gca, 'XTick', xticks, 'XTickLabel',transpose(val_label))

yticks = linspace(1, size(val_label,1),size(val_label,1));
set(gca, 'YTick', yticks, 'YTickLabel',transpose(val_label))
colorbar

%% histogram parameters gate duration

figure

hist(prePCA(:,1),length(prePCA(:,1)))
%hist(prePCA(:,1),100)
title('histogram of gate duration')
xlabel('duration')
ylabel('number')
%% PCA

[coefs,score, ~, ~, explained] = pca(PCA);

figure 
hold on
bar(explained)
plot(1:numel(explained), cumsum(explained), 'o-', 'MarkerFaceColor', 'r')
yyaxis right
h = gca;
h.YAxis(2).Limits = [0 100];
h.YAxis(2).Color = h.YAxis(1).Color;
h.YAxis(2).TickLabel = strcat(h.YAxis(2).TickLabel, '%');

pc1 = score(:,1);
pc2 = score(:,2);
pc3 = score(:,3);

figure
gscatter(pc1,pc2,color);

figure
gscatter(pc1,pc3,color);

figure
gscatter(pc2,pc3,color);

figure
scatter3(pc1,pc2,pc3, 1, color);
title('PCA')
xlabel('pc1')
ylabel('pc2')
zlabel('pc3')

figure
biplot(coefs(:,1:3),'VarLabels',val_label);

figure
biplot(coefs(:,1:2),'VarLabels',val_label);
%biplot(coefs(:,1:2),'Scores',score(:,1:2),'VarLabels',val_label);

%% Basic functions

% filtering : take a signal (ex the y value for a marker) and return a
% clean the signal S -> S_f
function [S_f] = filtering(S)

    %low, high pass filter
    d1 = designfilt("lowpassiir",FilterOrder=2, HalfPowerFrequency=0.03,DesignMethod="butter");
    S_f = filtfilt(d1,S);
end


% cut_gate : cut the date in gate cycle (using a gradient)
% each value represent a foot strike
% data is the dataset
% constrain : do we get rid of the gate that are aberentely long ? (ie. not
% well detected)
% Gate is a structure, each gate is represented but it events
function Gate = cut_gate(data,constraint) 
    T = 1/data.marker_sr;

    %filtering
     S_L = filtering(data.LTOE(:,2));
     S_R = filtering(data.RTOE(:,2));

%     filtering
%     S_L = data.LTOE(:,2);
%     S_R = data.RTOE(:,2);

    % calculate the gradient
    Gl = gradient(S_L);
    Gr = gradient(S_R);

    %initialisation
    Gate = [];
    gate = [];
    
    on_gate = false;

    % calculate the position
    i = 1;
    while i <=(length(Gl)-1)
        i = i+1;
        if sign(Gl(i)) ~= sign(Gl(i-1)) && sign(Gl(i)) < 0 
            
            if isfield(gate,'strikeR')
                gate.offnext = i;
                if not(constraint) || (constraint && (gate.offnext - gate.offL)*T <=4)
                    Gate = [Gate,gate];
                elseif constraint && gate.offnext - gate.offL >4 
                    disp(['remove a gate with duration',num2str(gate.offnext - gate.offL),'in position',num2str(gate.offL)])
                end
            end

            on_gate = true;
            gate = [];
            gate.offL = i;

            while on_gate && i <=(length(Gl)-1)
                i = i+1;
                if sign(Gl(i-1)) ~= sign(Gl(i)) && sign(Gl(i)) > 0 && isfield(gate,'offL')
                    gate.strikeL = i;
                elseif sign(Gr(i-1)) ~= sign(Gr(i)) && sign(Gr(i)) < 0 && isfield(gate,'strikeL')
                    gate.offR = i;
                elseif sign(Gr(i-1)) ~= sign(Gr(i)) && sign(Gr(i)) > 0 && isfield(gate,'offR')
                    gate.strikeR = i;
                    on_gate = false;
                end
            end

        end
    end

end

%% gate parameters functions

function parameters = calculated_parameters(data,S)

    display("start data")
    %cut the data in gate
    Gate = cut_gate(data,true);

    %initiate the structure of parameters
    parameters = [];

    %get the gate dependent parameters
    for i = 1:(length(Gate)) 
        parameters = [parameters,gate_parameters(data,Gate(i),S)];
    end
    
    display("finish single parameters")
    %get the global values (for exemple the mean of the swing duration)
    Global = get_global(parameters);
   
    display("finish global value")
    %get the global values dependent parameters (for exemple the
    %variability)
    parameters = gate_global_parameters(parameters,Global);

    display("finish global parameters")
    %removal 

    parameters = remove_param(parameters);
end

function parameter = gate_parameters(data,gate,S)
    
    %preli

    T = 1/data.marker_sr;
    speed = S*1000/3.6;

    %alone parameters

    parameter.gate_duration_sec = (gate.offnext-gate.offL)*T;

    parameter.swing_duration_left_sec = (gate.strikeL-gate.offL)*T;
    
    parameter.swing_duration_right_sec = (gate.strikeR-gate.offR)*T;

    %parameter.swing_duration_symetry = abs(parameter.swing_duration_right_sec - parameter.swing_duration_left_sec);
    parameter.swing_duration_symetry = 100*(parameter.swing_duration_right_sec - parameter.swing_duration_left_sec)/(0.5*(parameter.swing_duration_right_sec + parameter.swing_duration_left_sec));

    parameter.stance_percentage_left = (parameter.gate_duration_sec - parameter.swing_duration_left_sec)/parameter.gate_duration_sec;

    parameter.stance_percentage_right = (parameter.gate_duration_sec - parameter.swing_duration_right_sec)/parameter.gate_duration_sec;
    
    parameter.double_stance_percentage = 100*max(0,(parameter.gate_duration_sec-(parameter.swing_duration_left_sec+parameter.swing_duration_right_sec))/parameter.gate_duration_sec);
    
    parameter.stance_percentage_total = parameter.stance_percentage_right + parameter.stance_percentage_left - parameter.double_stance_percentage;
    
    parameter.step_height_left_mm = max(data.LANK(gate.offL:gate.offnext,3))-min(data.LANK(gate.offL:gate.offnext,3));
    
    parameter.step_length_left_mm = (abs(data.LTOE(gate.strikeL,2) -data.LTOE(gate.offL,2)) + parameter.swing_duration_left_sec*speed);
    
    parameter.step_height_right_mm = max(data.RANK(gate.offL:gate.offnext,3))-min(data.RANK(gate.offL:gate.offnext,3));
    
    parameter.step_length_right_mm = (abs(data.RTOE(gate.strikeR,2)-(data.RTOE(gate.offR,2))) + parameter.swing_duration_right_sec*speed);
    
    %parameter.step_height_symetry = abs(parameter.step_height_left_mm-parameter.step_height_right_mm);

    parameter.step_height_symetry = 100*(parameter.step_height_left_mm-parameter.step_height_right_mm)/(0.5*(parameter.step_height_left_mm+parameter.step_height_right_mm));

    %parameter.step_length_symetry = abs(parameter.step_length_left_mm-parameter.step_length_right_mm);
    parameter.step_length_symetry = 100*(parameter.step_length_left_mm-parameter.step_length_right_mm)/(0.5*(parameter.step_length_left_mm+parameter.step_length_right_mm));


    parameter.stridewidth_mm = mean(data.LANK(gate.offL:gate.offnext,1) - data.RANK(gate.offL:gate.offnext,1));


    x1L = data.LKNE(gate.offL:gate.offnext,2);
    x3L = data.LHIP(gate.offL:gate.offnext,2);
    x2L = data.LANK(gate.offL:gate.offnext,2);
    
    y1L = data.LKNE(gate.offL:gate.offnext,3);
    y3L = data.LHIP(gate.offL:gate.offnext,3);
    y2L = data.LANK(gate.offL:gate.offnext,3);

    x1R = data.RKNE(gate.offL:gate.offnext,2);
    x3R = data.RHIP(gate.offL:gate.offnext,2);
    x2R = data.RANK(gate.offL:gate.offnext,2);
    
    y1R = data.RKNE(gate.offL:gate.offnext,3);
    y3R = data.RHIP(gate.offL:gate.offnext,3);
    y2R = data.RANK(gate.offL:gate.offnext,3);

    x_hip_latL = data.LHIP(gate.offL:gate.offnext,1);
    mean_x_hip_latL = mean(data.LHIP(:,1));
    x_hip_latR = data.RHIP(gate.offL:gate.offnext,1);
    mean_x_hip_latR = mean(data.RHIP(:,1));

    mean_y3L = mean(data.LHIP(:,3));
    mean_y3R = mean(data.RHIP(:,3));

    angleL = atan2(y3L - y1L, x3L - x1L) - atan2(y2L - y1L, x2L - x1L);
    angle_radL = angleL*(180/pi);
    vel_angleL = gradient(angle_radL);

    angleR = atan2(y3R - y1R, x3R - x1R) - atan2(y2R - y1R, x2R - x1R);
    angle_radR = angleR*(180/pi);
    vel_angleR = gradient(angle_radR);

    % variability is defined as the difference with the mean value
    parameter.var_lateral_hip_left_mm = mean(abs(mean_x_hip_latL-x_hip_latL))/mean_x_hip_latL*100;
    parameter.var_lateral_hip_right_mm = mean(abs(mean_x_hip_latR-x_hip_latR))/mean_x_hip_latR*100;
    parameter.var_ver_hip_left_mm = mean(abs(mean_y3L-y3L))/mean_y3L*100;
    parameter.var_ver_hip_right_mm = mean(abs(mean_y3R-y3R))/mean_y3R*100;

    %joint angle
    parameter.max_joint_angle_left_deg = max(angle_radL);
    parameter.min_joint_angle_left_deg = min(angle_radL);
    parameter.max_angle_vel_left_deg = max(abs(vel_angleL));

    
    parameter.var_ver_hip_right_mm = var(y3R);
    parameter.max_joint_angle_right_deg = max(angle_radR);
    parameter.min_joint_angle_right_deg = min(angle_radR);
    parameter.max_angle_vel_right_deg = max(abs(vel_angleR));

    % EMG

    t_idx_emg = ((gate.offL-1)*(data.EMG_sr/data.marker_sr)+1:gate.offnext*data.EMG_sr/data.marker_sr);
    L_ag  = data.LSol(t_idx_emg);
    L_ant = data.LTA(t_idx_emg);

    R_ag  = data.RSol(t_idx_emg);
    R_ant = data.RTA(t_idx_emg);
    
    L_CI = emgLib.coactivation_index(L_ant,L_ag);
    R_CI = emgLib.coactivation_index(R_ant,R_ag);
    
    parameter.coactivation_index = (L_CI + R_CI)/2;

    
    EMG_filtered = emgLib.filter_emg(data.LSol(t_idx_emg),data.EMG_sr,0);
    %[parameter.amplitude_emg_LSol,parameter.integral_emg_LSol,parameter.rms_emg_LSol] = emgLib.emg_parameters(EMG_filtered,data.EMG_sr);
    [parameter.amplitude_emg_LSol,~, ~] = emgLib.emg_parameters(EMG_filtered,data.EMG_sr);
 
    EMG_filtered = emgLib.filter_emg(data.RSol(t_idx_emg),data.EMG_sr,0);
    %[parameter.amplitude_emg_RSol,parameter.integral_emg_RSol,parameter.rms_emg_RSol] = emgLib.emg_parameters(EMG_filtered,data.EMG_sr);
    [parameter.amplitude_emg_RSol,~, ~] = emgLib.emg_parameters(EMG_filtered,data.EMG_sr);

    EMG_filtered_flexor_R = emgLib.filter_emg(data.RTA(t_idx_emg),data.EMG_sr,0);
    %[parameter.amplitude_emg_RTA,parameter.integral_emg_RTA,parameter.rms_emg_RTA] = emgLib.emg_parameters(EMG_filtered,data.EMG_sr);
    [parameter.amplitude_emg_RTA,~, ~] = emgLib.emg_parameters(EMG_filtered_flexor_R,data.EMG_sr);

    EMG_filtered_flexor_L = emgLib.filter_emg(data.LTA(t_idx_emg),data.EMG_sr,0);
    %[parameter.amplitude_emg_LTA,parameter.integral_emg_LTA,parameter.rms_emg_LTA] = emgLib.emg_parameters(EMG_filtered,data.EMG_sr);
    [parameter.amplitude_emg_LTA,~, ~] = emgLib.emg_parameters(EMG_filtered_flexor_L,data.EMG_sr);

    EMG_filtered_extensor_L = emgLib.filter_emg(data.LMG(t_idx_emg),data.EMG_sr,0);
    [parameter.amplitude_emg_LMG,~,~] = emgLib.emg_parameters(EMG_filtered_extensor_L,data.EMG_sr);
    
    EMG_filtered_extensor_R = emgLib.filter_emg(data.RMG(t_idx_emg),data.EMG_sr,0);
    [parameter.amplitude_emg_RMG,~,~] = emgLib.emg_parameters(EMG_filtered_extensor_R,data.EMG_sr);


    % Burst
    parameter.extensorL_burst_duration = emgLib.calculate_burst_duration(EMG_filtered_extensor_L,data.EMG_sr,0);
    parameter.extensorR_burst_duration = emgLib.calculate_burst_duration(EMG_filtered_extensor_R,data.EMG_sr,0);
    
    parameter.flexorL_burst_duration = emgLib.calculate_burst_duration(EMG_filtered_flexor_L,data.EMG_sr,0);
    parameter.flexorR_burst_duration = emgLib.calculate_burst_duration(EMG_filtered_flexor_R,data.EMG_sr,0);


%     EMG_filtered = emgLib.filter_emg(data.LMG(t_idx_emg),data.EMG_sr,0);
%     [parameter.amplitude_emg_LMG,parameter.integral_emg_LMG,parameter.rms_emg_LMG] = emgLib.emg_parameters(EMG_filtered,data.EMG_sr);
% 
%     
%     EMG_filtered = emgLib.filter_emg(data.LST(t_idx_emg),data.EMG_sr,0);
%     [parameter.amplitude_emg_LST,parameter.integral_emg_LST,parameter.rms_emg_LST] = emgLib.emg_parameters(EMG_filtered,data.EMG_sr);
% 
%     %EMG_filtered = emgLib.filter_emg(data.LVLat(t_idx_emg),data.EMG_sr,0);
%     %[parameter.amplitude_emg_LVLat,parameter.integral_emg_LVLat,parameter.rms_emg_LVLat] = emgLib.emg_parameters(EMG_filtered,data.EMG_sr);
% 
%     EMG_filtered = emgLib.filter_emg(data.LRF(t_idx_emg),data.EMG_sr,0);
%     [parameter.amplitude_emg_LRF,parameter.integral_emg_LRF,parameter.rms_emg_LRF] = emgLib.emg_parameters(EMG_filtered,data.EMG_sr);
% 
%     %EMG_filtered = emgLib.filter_emg(data.LIl(t_idx_emg),data.EMG_sr,0);
%     %[parameter.amplitude_emg_LIl,parameter.integral_emg_LIl,parameter.rms_emg_LIl] = emgLib.emg_parameters(EMG_filtered,data.EMG_sr);

    
end

function Global = get_global(parameters)

    Global = [];

    mean_gate_duration = 0;

    mean_swing_left = 0;
    mean_swing_right = 0;

    mean_step_height_left = 0;
    mean_step_height_right = 0;

    mean_step_length_left = 0;
    mean_step_length_right = 0;

    mean_stridewidth = 0;

    mean_coactivation_index = 0;
    mean_amplitude_emg_LSol      = 0;
    mean_amplitude_emg_LTA      = 0;
    mean_amplitude_emg_RSol      = 0;
    mean_amplitude_emg_RTA      = 0;

%     mean_amplitude_emg_LMG      = 0;
%     mean_amplitude_emg_LST      = 0;
%     %mean_amplitude_emg_LVLat      = 0;
%     mean_amplitude_emg_LRF      = 0;
%     %mean_amplitude_emg_LIl     = 0;

    mean_flexorL_burst_duration = 0;
    mean_extensorL_burst_duration = 0;
    mean_flexorR_burst_duration = 0;
    mean_extensorR_burst_duration = 0;


    if (not(isempty(parameters)))
        for i = 1:(length(parameters))

            mean_gate_duration = mean_gate_duration + parameters(i).gate_duration_sec;

            mean_swing_left = mean_swing_left + parameters(i).swing_duration_left_sec;
            mean_swing_right = mean_swing_right + parameters(i).swing_duration_right_sec;

            mean_step_height_left = mean_step_height_left + parameters(i).step_height_left_mm;
            mean_step_height_right = mean_step_height_right + parameters(i).step_height_right_mm;

            mean_step_length_left = mean_step_length_left + parameters(i).step_length_left_mm;
            mean_step_length_right = mean_step_length_right + parameters(i).step_length_right_mm;

            mean_stridewidth = mean_stridewidth + parameters(i).stridewidth_mm;

            
            % EMG
            mean_coactivation_index = mean_coactivation_index + parameters(i).coactivation_index;
            
            mean_amplitude_emg_LSol = mean_amplitude_emg_LSol + parameters(i).amplitude_emg_LSol;
            mean_amplitude_emg_LTA = mean_amplitude_emg_LTA + parameters(i).amplitude_emg_LTA;
            mean_amplitude_emg_RSol = mean_amplitude_emg_RSol + parameters(i).amplitude_emg_RSol;
            mean_amplitude_emg_RTA = mean_amplitude_emg_RTA + parameters(i).amplitude_emg_RTA;
            
            

%             mean_amplitude_emg_LMG = mean_amplitude_emg_LMG + parameters(i).amplitude_emg_LMG;
%             mean_amplitude_emg_LST = mean_amplitude_emg_LST + parameters(i).amplitude_emg_LST;
%             %mean_amplitude_emg_LVLat = mean_amplitude_emg_LVLat + parameters(i).amplitude_emg_LVLat;
%             mean_amplitude_emg_LRF = mean_amplitude_emg_LRF + parameters(i).amplitude_emg_LRF;
%             %mean_amplitude_emg_LIl = mean_amplitude_emg_LIl + parameters(i).amplitude_emg_LIl;

            mean_flexorL_burst_duration = mean_flexorL_burst_duration + parameters(i).flexorL_burst_duration;
            mean_extensorL_burst_duration = mean_extensorL_burst_duration + parameters(i).extensorL_burst_duration;
            
            
            mean_flexorR_burst_duration = mean_flexorR_burst_duration + parameters(i).flexorR_burst_duration;
            mean_extensorR_burst_duration = mean_extensorR_burst_duration + parameters(i).extensorR_burst_duration;
            
            


        end

        Global.mean_gate_duration = mean_gate_duration/length(parameters);
        Global.mean_swing_left = mean_swing_left/length(parameters);
        Global.mean_swing_right = mean_swing_right/length(parameters);

        Global.mean_step_height_left = mean_step_height_left/length(parameters);
        Global.mean_step_height_right = mean_step_height_right/length(parameters);

        Global.mean_step_length_left = mean_step_length_left/length(parameters);
        Global.mean_step_length_right = mean_step_length_right/length(parameters);

        Global.mean_stridewidth = mean_stridewidth/length(parameters);

        %EMG
        Global.mean_coactivation_index = mean_coactivation_index/length(parameters);
        Global.mean_amplitude_emg_LSol      = mean_amplitude_emg_LSol/length(parameters);
        Global.mean_amplitude_emg_LTA      = mean_amplitude_emg_LTA/length(parameters);
        Global.mean_amplitude_emg_RSol      = mean_amplitude_emg_RSol/length(parameters);
        Global.mean_amplitude_emg_RTA      = mean_amplitude_emg_RTA/length(parameters);

        Global.mean_flexorL_burst_duration = mean_flexorL_burst_duration/length(parameters);
        Global.mean_extensorL_burst_duration = mean_extensorL_burst_duration/length(parameters);
        
        Global.mean_flexorR_burst_duration = mean_flexorR_burst_duration/length(parameters);
        Global.mean_extensorR_burst_duration = mean_extensorR_burst_duration/length(parameters);
        

%         Global.mean_amplitude_emg_LMG      = mean_amplitude_emg_LMG/length(parameters);
%         Global.mean_amplitude_emg_LST      = mean_amplitude_emg_LST/length(parameters);
%         %Global.mean_amplitude_emg_LVLat      = mean_amplitude_emg_LVLat/length(parameters);
%         Global.mean_amplitude_emg_LRF      = mean_amplitude_emg_LRF/length(parameters);
%         %Global.mean_amplitude_emg_LIl      = mean_amplitude_emg_LIl/length(parameters);
%         % LTA_filtered = emgLib.filter_emg(data.LTA,0);
%         % [burst_duration,Global.mean_amplitude_emg,Global.integral_emg,Global.rms_emg] = emgLib.emg_parameters(LTA_filtered,0);
%         % Global.mean_burst_duration = mean(burst_duration);

    end

end

function parameters = gate_global_parameters(parameters,Global)

    %get the global values dependent parameters (for exemple the
    %variability)
    for i = 1:(length(parameters))
         parameters(i).var_step_height_left_mm = abs(Global.mean_step_height_left- parameters(i).step_height_left_mm)/Global.mean_step_height_left*100;
         parameters(i).var_step_height_right_mm = abs(Global.mean_step_height_right- parameters(i).step_height_right_mm)/Global.mean_step_height_right*100;
         parameters(i).var_step_length_left_mm = abs(Global.mean_step_length_left- parameters(i).step_length_left_mm)/Global.mean_step_length_left*100;
         parameters(i).var_step_length_right_mm = abs(Global.mean_step_length_right- parameters(i).step_length_right_mm)/Global.mean_step_length_right*100;
         parameters(i).var_stridewidth = abs(Global.mean_stridewidth- parameters(i).stridewidth_mm)/Global.mean_stridewidth*100;

         
         %EMG 
         parameters(i).var_coactivation_index = abs(Global.mean_coactivation_index-parameters(i).coactivation_index)/Global.mean_coactivation_index*100;
         
         parameters(i).var_amplitude_emg_LSol = abs(Global.mean_amplitude_emg_LSol-parameters(i).amplitude_emg_LSol)/Global.mean_amplitude_emg_LSol*100;
         parameters(i).var_amplitude_emg_LTA = abs(Global.mean_amplitude_emg_LTA-parameters(i).amplitude_emg_LTA)/Global.mean_amplitude_emg_LTA*100;
         parameters(i).var_amplitude_emg_RSol = abs(Global.mean_amplitude_emg_RSol-parameters(i).amplitude_emg_RSol)/Global.mean_amplitude_emg_RSol*100;
         parameters(i).var_amplitude_emg_RTA = abs(Global.mean_amplitude_emg_RTA-parameters(i).amplitude_emg_RTA)/Global.mean_amplitude_emg_RTA*100;
         
         parameters(i).var_flexorL_burst_duration = abs(Global.mean_flexorL_burst_duration-parameters(i).flexorL_burst_duration)/Global.mean_flexorL_burst_duration;
         parameters(i).var_extensorL_burst_duration = abs(Global.mean_extensorL_burst_duration-parameters(i).extensorL_burst_duration)/Global.mean_extensorL_burst_duration;
         
         parameters(i).var_flexorR_burst_duration = abs(Global.mean_flexorR_burst_duration-parameters(i).flexorR_burst_duration)/Global.mean_flexorR_burst_duration;
         parameters(i).var_extensorR_burst_duration = abs(Global.mean_extensorR_burst_duration-parameters(i).extensorR_burst_duration)/Global.mean_extensorR_burst_duration;

         
%          parameters(i).var_amplitude_emg_LMG = abs(Global.mean_amplitude_emg_LMG-parameters(i).amplitude_emg_LMG)/Global.mean_amplitude_emg_LMG*100;
%          parameters(i).var_amplitude_emg_LST = abs(Global.mean_amplitude_emg_LST-parameters(i).amplitude_emg_LST)/Global.mean_amplitude_emg_LST*100;
%          %parameters(i).var_amplitude_emg_LVLat = abs(Global.mean_amplitude_emg_LVLat-parameters(i).amplitude_emg_LVLat)/Global.mean_amplitude_emg_LVLat*100;
%          parameters(i).var_amplitude_emg_LRF = abs(Global.mean_amplitude_emg_LRF-parameters(i).amplitude_emg_LRF)/Global.mean_amplitude_emg_LRF*100;
%          %parameters(i).var_amplitude_emg_LIl = abs(Global.mean_amplitude_emg_LIl-parameters(i).amplitude_emg_LIl)/Global.mean_amplitude_emg_LIl*100;

    end
end

function parameters =  remove_param(parameters)

    %parameters = rmfield(parameters,"gate_duration_sec");

    %parameters = rmfield(parameters,"step_length_right_mm");
    %parameters = rmfield(parameters,"step_length_left_mm");
    %parameters = rmfield(parameters,"step_height_right_mm");
    %parameters = rmfield(parameters,"step_height_left_mm");
    
end
