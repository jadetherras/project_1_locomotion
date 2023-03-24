
%% Kinematic file

% In this file, we analize the kinematic data
% filter and plot the signal
% cut in gate cycle
% animation of the movement
% plot the gate cycle

%to do :

% calculate the parameters using function

%% Loading the data

%dataset sain 
data_healthy=load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/4_AML01_3kmh.mat");

% dataset SCI Human
data_SCI=load("SCI Human/DM002_TDM_08_1kmh.mat");

%% Plot the movement

% animation : decomment the one you want (healthy, SCI)

%plot_t(data_healthy.data)
%plot_t(data_SCI.data)

% plot a gate : chose the marker for the gate calculation and a gate

%filter the data (here toe) and gate calculation
S_L = filtering(data_healthy.data.LTOE(:,2));
[time_L,pos] = gate(S_L);
N1 = 191; %in pos
N2 = 356;

%plot the gate cycle
plot_gate(data_healthy.data,N1,N2,'B')

%% function

% filtering : take a signal (ex the y value for a marker) and return a
% clean signa
function [S_f] = filtering(S)

    %low, high pass filter
    S_f = lowpass(highpass(S,1e-1,1e2),0.6,1e2, 'ImpulseResponse','iir');

end
    
% gate : take a signal and calculate the gate cycle (using the gradient)
% time take the value foot off and foot strike at the correct timing
% pos is the list of the position of the gate events 
function [time,pos] = gate(S) 

    % calculate the gradient
    G = gradient(S);
    
    %initialisation
    time = {''};
    pos = {};

    % calculate the position
    for i = 2:length(G)
        if sign(G(i)) ~= sign(G(i-1)) % if change sign, pass 0
            if sign(G(i)) > 0
                time = [time,"foot strike"]; 
            else
                time = [time,"foot off"];
            end

            pos = [pos,i];
        else
            time = [time,""];
        end
    end

end

% plot_gate : take the dataset, a starting time and ending time (ex
% start/end of a gate cycle)
% plot the gate as a line for each time 
% data is the dataset
% N1 is the start position
% N2 is the end position
% side is a string giving the limb to plot => L/R/B (left, right,both)
function plot_gate(data,N1,N2,side)
    
    %axis([1 1200 0 1200])
    xlabel('x'), ylabel('y')
    title('kinematic reconstruction')

    
    for n = N1:N2
        %dec = (n-N1)*2;

        if side == 'L' || side == 'B'

            plot([data.LHIP(n,2) data.LKNE(n,2)],[data.LHIP(n,3) data.LKNE(n,3)],'black')
            if n == N1 
                hold on
            end
            plot([data.LKNE(n,2) data.LANK(n,2)],[data.LKNE(n,3) data.LANK(n,3)],'black')
            plot([data.LANK(n,2) data.LTOE(n,2)],[data.LANK(n,3) data.LTOE(n,3)],'black')
        end

        if side == 'R' || side == 'B'

            plot([data.RHIP(n,2) data.RKNE(n,2)],[data.RHIP(n,3) data.RKNE(n,3)],'green')
            if n == N1 && side == 'R'
                hold on
            end
            plot([data.RKNE(n,2) data.RANK(n,2)],[data.RKNE(n,3) data.RANK(n,3)],'green')
            plot([data.RANK(n,2) data.RTOE(n,2)],[data.RANK(n,3) data.RTOE(n,3)],'green')
        end

    end

    hold off
end

% plot_t : plot a animation of a dataset
% the left legs is represented by full circular markers
% the left legs is represented by empty circular markers
% on the corner, a text display "foot strike/foot off" at the correct time
% the toe is use to calculate the gate
% data is the dataset
function plot_t(data)

    %marker = ["LHIP","LKNE", "LANK","LTOE"];
    N = length(data.LHIP(:,1));
    T = 1/120;

    % calculate timing left and right 
    S_L = filtering(data.LTOE(:,2));
    time_L = gate(S_L);

    S_R = filtering(data.RTOE(:,2));
    time_R = gate(S_R);


    %left plot
    hip = scatter(data.LHIP(1,2),data.LHIP(1,3),'o','MarkerFaceColor','red');
    hold on 
    knee = scatter(data.LKNE(1,2),data.LKNE(1,3),'o','MarkerFaceColor','blue');
    ankle = scatter(data.LANK(1,2),data.LANK(1,3),'o','MarkerFaceColor','magenta');
    toe = scatter(data.LTOE(1,2),data.LTOE(1,3),'o','MarkerFaceColor','green');
    

    %right plot
    hip2 = scatter(data.RHIP(1,2),data.RHIP(1,3),'o','red');
    knee2 = scatter(data.RKNE(1,2),data.RKNE(1,3),'o','blue');
    ankle2 = scatter(data.RANK(1,2),data.RANK(1,3),'o','magenta');
    toe2 = scatter(data.RTOE(1,2),data.RTOE(1,3),'o','green');

    hold off 

    axis([1 1200 0 1200])
    xlabel('x'), ylabel('y')
    title('kinematic reconstruction')

    %update the markers and text at each time
    for n = 2:N

        move_L = text(1000,1000,"left : " + time_L(n));
        move_R = text(1000,900,"right : " + time_R(n));

        hip.XData = data.LHIP(n,2);
        hip.YData = data.LHIP(n,3);
        knee.XData = data.LKNE(n,2);
        knee.YData = data.LKNE(n,3);
        ankle.XData = data.LANK(n,2);
        ankle.YData = data.LANK(n,3);
        toe.XData = data.LTOE(n,2);
        toe.YData = data.LTOE(n,3);

        hip2.XData = data.RHIP(n,2);
        hip2.YData = data.RHIP(n,3);
        knee2.XData = data.RKNE(n,2);
        knee2.YData = data.RKNE(n,3);
        ankle2.XData = data.RANK(n,2);
        ankle2.YData = data.RANK(n,3);
        toe2.XData = data.RTOE(n,2);
        toe2.YData = data.RTOE(n,3);
        
        % give the time to see before change
        drawnow
        pause(T)

        delete(move_L)
        delete(move_R)
    end

end

%% try to plot markers 

% idea : plot the marker in 2d (not x) moving and add the timer on it ! 
% good way to verify but maybe not simple

%to do : 
%accuracy comparison and do it for all the data + cut and plot