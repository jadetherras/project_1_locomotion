
%dataset sain 
data=load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/1_AML01_2kmh.mat");

[n,p] = size(data.data.LSol)
t = 1:n;
t = t/1000;

%% EMG data plot

plot(t,[data.data.LSol,data.data.LMG,data.data.LTA,data.data.LST,data.data.LVLat,data.data.LRF,data.data.LIl,data.data.RSol,data.data.RMG,data.data.RTA,data.data.RST,data.data.RVLat,data.data.RRF]),
legend('LSol','LMG','LTA','LST','LVlat','LRF','Lil','RSol','RMG','RTA','RST','RVLat','RRF')

xlabel('Time'), ylabel('mV')
title('plot')

%% Plot of the movement

[n2,t2] = size(data.data.LHIP);
 
t2 = 1:n2;
t2 = t2/120;

%global
plot(t2, [data.data.LHIP,data.data.LKNE,data.data.LANK,data.data.LTOE,data.data.RHIP,data.data.RKNE,data.data.RANK,data.data.RTOE])
legend('LHIP','LKNE','LANK', 'LTOE','RHIP','RKNE','RANK','RTOE')
xlabel('Time'), ylabel('movement')
title('plot')

%left
plot(t2, [data.data.LHIP,data.data.LKNE,data.data.LANK,data.data.LTOE])
legend('LHIP','LKNE','LANK', 'LTOE')
xlabel('Time'), ylabel('movement')
title('plot')

%hip
figure
plot(t2, [data.data.LHIP])
legend('x','y','z')
xlabel('Time'), ylabel('Displacement')
title('plot hip marker')
savefig('figure/H_1_2kmh_hip.fig')

%knee 
figure
plot(t2, [data.data.LKNE])
legend('x','y','z')
xlabel('Time'), ylabel('Displacement')
title('plot knee marker')
savefig('figure/H_1_2kmh_knee.fig')

figure
plot(data.data.LKNE(:,2), data.data.LKNE(:,3))
legend('movement')
xlabel('front/back'), ylabel('up/down')
title('plot knee marker')
savefig('figure/H_1_2kmh_knee_yz.fig')

%ankle
figure
plot(t2, [data.data.LANK])
legend('x','y','z')
xlabel('Time'), ylabel('Displacement')
title('plot ankle marker')
savefig('figure/H_1_2kmh_ankle.fig')

figure
plot(data.data.LANK(:,2), data.data.LANK(:,3)-265)
axis([250 900 -325 325])
legend('movement')
xlabel('front/back'), ylabel('up/down')
title('plot ankle marker')
savefig('figure/H_1_2kmh_ankle_yz.fig')

%toe
figure
plot(t2, [data.data.LTOE])
legend('x','y','z')
xlabel('Time'), ylabel('Displacement')
title('plot toe marker')
savefig('figure/H_1_2kmh_toe.fig')

figure
plot(data.data.LTOE(:,2), data.data.LTOE(:,3))
axis([50 800 -100 650])
legend('movement')
xlabel('front/back'), ylabel('up/down')
title('plot toe marker')
savefig('figure/H_1_2kmh_toe_yz.fig')

%% dataset SCI Human
%data=load("SCI Human/DM002_TDM_08_1kmh.mat");

%[n2,p2] = size(data.data.LHIP)
%
%t2 = 1:n2;
%t2 = t2/100;

%% normal plotting

%hip
figure
plot(t2, [data.data.LHIP])
legend('x','y','z')
xlabel('Time'), ylabel('Displacement')
title('plot hip marker')
savefig('figure/SCI_1kmh_hip.fig')

figure
plot(data.data.LHIP(:,2), data.data.LHIP(:,3))
legend('movement')
xlabel('front/back'), ylabel('up/down')
title('plot hip marker')
savefig('figure/SCI_1kmh_hip_yz.fig')

%knee 
figure
plot(t2, data.data.LKNE)
legend('x','y','z')
xlabel('Time'), ylabel('Displacement')
title('plot knee marker')
savefig('figure/SCI_1kmh_knee.fig')

figure
plot(data.data.LKNE(:,2), data.data.LKNE(:,3))
legend('movement')
xlabel('front/back'), ylabel('up/down')
title('plot knee marker')
savefig('figure/SCI_1kmh_knee_yz.fig')

%ankle
figure
plot(t2, [data.data.LANK])
legend('x','y','z')
xlabel('Time'), ylabel('Displacement')
title('plot ankle marker')
savefig('figure/SCI_1kmh_ankle.fig')

figure
plot(data.data.LANK(:,2), data.data.LANK(:,3)-265)
%axis([250 900 -325 325])
legend('movement')
xlabel('front/back'), ylabel('up/down')
title('plot ankle marker')
savefig('figure/SCI_1kmh_ankle_yz.fig')

%toe
figure
plot(t2, [data.data.LTOE])
legend('x','y','z')
xlabel('Time'), ylabel('Displacement')
title('plot toe marker')
savefig('figure/SCI_1kmh_toe.fig')

figure
plot(data.data.LTOE(:,2), data.data.LTOE(:,3))
legend('movement')
xlabel('front/back'), ylabel('up/down')
title('plot toe marker')
savefig('figure/SCI_1kmh_toe_yz.fig')
%}


%% filter the signal without function

y = highpass(data.data.LTOE(:,2),1e-1,1e2);
y_toe = lowpass(y,0.6,1e2, 'ImpulseResponse','iir');

G = gradient(y_toe);
G2 = gradient(G);

figure
plot(t2,[data.data.LTOE(:,2)/100,y_toe/100,G,G2])
legend('ori','y','grad y','2grad y')
xlabel('Time'), ylabel('Displacement')
title('plot toe marker')
savefig('figure/SCI_1kmh_toe_filter.fig')

% for knee

y = highpass(data.data.LKNE(:,2),1e-1,1e2);
y_knee = lowpass(y,0.6,1e2, 'ImpulseResponse','iir');

figure
plot(t2,[data.data.LKNE(:,2),y_knee])
legend('ori','y')
xlabel('Time'), ylabel('Displacement')
title('plot knee marker')
savefig('figure/SCI_1kmh_knee_filter.fig')

% comparison knee, toe
figure
plot(t2,[y_toe,y_knee])
legend('toe','knee')
xlabel('Time'), ylabel('Displacement')
title('plot toe-knee marker')
savefig('figure/SCI_1kmh_toe_knee_filter.fig')


time = {};
for i = 2:(length(y_toe)-1)
    if sign(y_toe(i+1)-y_toe(i)) ~= sign(y_toe(i)-y_toe(i-1))
        time = [time,i/100];
    end
end

time = {};
col = {};
for i = 2:length(G)
    if sign(G(i)) ~= sign(G(i-1))
        time = [time,i/100];
        if sign(G(i)) > 0
            col = [col,'r'];
        else
            col = [col,'b'];
        end
    end
end

figure
plot(t2,y_toe)
for k = 1:numel(time)
    xline(time{k},col{k})
end
legend('y filtered')
xlabel('Time'), ylabel('Displacement')
title('plot toe marker')
savefig('figure/SCI_1kmh_toe_filter.fig')

%% with function

% toe
[S_f_toe,time_toe,col_toe] = filtering(data.data.LTOE(:,2));

figure
plot(t2,S_f_toe)
for k = 1:numel(time)
    xline(time{k},col{k})
end
legend('y filtered')
xlabel('Time'), ylabel('Displacement')
title('plot toe marker')
savefig('figure/SCI_1kmh_toe_filter.fig')

%knee
[S_f_knee,time_knee,col_knee] = filtering(data.data.LKNE(:,2));

figure
plot(t2,S_f_knee)
for k = 1:numel(time)
    xline(time{k},col{k})
end
legend('y filtered')
xlabel('Time'), ylabel('Displacement')
title('plot knee marker')
savefig('figure/SCI_1kmh_toe_filter.fig')


figure
plot(t2,[S_f_toe,S_f_knee])
legend('toe','knee')
xlabel('Time'), ylabel('Displacement')
title('plot toe-knee marker')
savefig('figure/SCI_1kmh_toe_knee_filter.fig')

%accuracy(time_toe,time_knee)

%% accuracy manager 

function accuracy(toe,knee)
    
    figure
    
    diff = {};
    for i = 1:numel(toe)
        scatter(toe{i},knee{i})
        diff = [diff, abs(toe{i}-knee{i})];
    end
    xlabel('toe'), ylabel('knee')
    title('accuracy')
    
    print(mean(diff))

end

%% filter and calculate pos and calculate gate

function [S_f,time,col] = filtering(S)
    S_f = lowpass(highpass(S,1e-1,1e2),0.6,1e2, 'ImpulseResponse','iir');
    
    G = gradient(S_f);
    time = {};
    col = {};

    for i = 2:length(G)
        if sign(G(i)) ~= sign(G(i-1))
            time = [time,i/100];
            if sign(G(i)) > 0
                col = [col,'r'];
            else
                col = [col,'b'];
            end
        end
    end

end

%% try to plot markers 

% idea : plot the marker in 2d (not x) moving and add the timer on it ! 
% good way to verify but maybe not simple

%to do : 
%function plot 
%fonction more generalised filter, grad, find
%accuracy comparison and do it for all the data + cut and plot