
%% Loading the data

close all;

%dataset 
%data=load("SCI Human/DM002_TDM_1kmh_NoEES.mat");
%name = ' no EES';
%data=load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/3_AML01_1kmh.mat");
%name = ' healthy 1km 01';
data=load("SCI Human/DM002_TDM_08_2kmh.mat");
name = ' ESS 2kmh';

[n,p] = size(data.data.LSol);
t = 1:n;
t = t/data.data.EMG_sr;


%% EMG data plot - exploring dataset

figure

plot(t,data.data.LSol)
hold on
plot(t,data.data.LMG+3*10^-4)
legend('sol','MG')
xlabel('Time'), ylabel('EMG')
title(strcat('EMG data ',name))
saveas(gcf,strcat('figure/EMG_SolvsMG ',name,'.png'))
%% ECG data plot - exploring dataset

[n,p] = size(data.data.LKNE);
t = 1:n;
t = t/data.data.marker_sr;

 figure
 plot(t, [data.data.LTOE])
 legend('x','y','z')
 xlabel('Time'), ylabel('Displacement')
 title(strcat('toe marker vs time ',name))
 saveas(gcf,strcat('figure/toe_time ',name,'.png'))

 figure
 plot(t, [data.data.LANK])
 legend('x','y','z')
 xlabel('Time'), ylabel('Displacement')
 title(strcat('ankle marker vs time ',name))
 saveas(gcf,strcat('figure/ankle_time ',name,'.png'))
 

 figure
 plot(data.data.LTOE(:,2), data.data.LTOE(:,3))
 legend('movement')
 xlabel('front/back'), ylabel('up/down')
 title(strcat('toe marker ',name))
 saveas(gcf,strcat('figure/toe_xy ',name,'.png'))

 figure
 plot(data.data.LANK(:,2), data.data.LANK(:,3))
 legend('movement')
 xlabel('front/back'), ylabel('up/down')
 title(strcat('ankle marker ',name))
 saveas(gcf,strcat('figure/ankle_xy ',name,'.png'))

 figure
 plot3(data.data.LANK(:,1),data.data.LANK(:,2), data.data.LANK(:,3))
 legend('movement')
 ylabel('left/right'),ylabel('front/back'), zlabel('up/down')
 title(strcat('plot ankle marker ',name))
 saveas(gcf,strcat('figure/ankle_xyz ',name,'.png'))

 figure
 plot3(data.data.LTOE(:,1),data.data.LTOE(:,2), data.data.LTOE(:,3))
 legend('movement')
 ylabel('left/right'),ylabel('front/back'), zlabel('up/down')
 title(strcat('plot ankle marker ',name))
 saveas(gcf,strcat('figure/toe_xyz ',name,'.png'))

%% filter the signal

 figure
 plot(t,[data.data.LANK(:,2),base(data.data.LANK(:,2)),filtering(data.data.LANK(:,2))])
 legend('original','first filter','filtfilt')
 xlabel('Time'), ylabel('Displacement')
 title(strcat('plot ankle marker filtered ',name))
 saveas(gcf,strcat('figure/filter ',name,'.png'))


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


function [S_f] = base(S)

   S_f = lowpass(S,0.6,1e2, 'ImpulseResponse','iir');

end

% filtering : take a signal (ex the y value for a marker) and return a
% clean signa
function [S_f] = filtering(S)

    %low, high pass filter
    d1 = designfilt("lowpassiir",FilterOrder=2, HalfPowerFrequency=0.03,DesignMethod="butter");
    S_f = filtfilt(d1,highpass(S,1e-1,1e2));
end

% calculation : take a signal and calculate the gate cycle (using the gradient
function [S_f,time,col] = calculation(S)
    S_f = filtering(S);
    
    G = gradient(S_f);
    time = {};
    col = {};

    for i = 2:length(G)
        if sign(G(i)) ~= sign(G(i-1))
            time = [time,i];
            if sign(G(i)) > 0
                col = [col,'r'];
            else
                col = [col,'b'];
            end
        end
    end

end