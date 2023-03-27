
%% gait parameters file

% In this file, we calculate some gait parameters

%% Loading the data

%dataset sain 
data_healthy=load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/4_AML01_3kmh.mat");

% dataset SCI Human
data_SCI=load("SCI Human/DM002_TDM_08_1kmh.mat");

%% gate parameters

gate = cut_gate(data_healthy.data);

function swing_duration()

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
    S_f = filtering(data.LTOE(:,3));

    % calculate the gradient
    G = gradient(S_f);

    %initialisation
    gate = [];

    % calculate the position
    for i = 2:length(G)
        if sign(G(i)) ~= sign(G(i-1)) && sign(G(i)) > 0 
            % if change sign, pass 0 and >0 foot strike
                gate = [gate,i];
        end
    end

end