
%dataset sain 
data_healthy=load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/1_AML01_2kmh.mat");

% dataset SCI Human
data_SCI=load("SCI Human/DM002_TDM_08_1kmh.mat");

%% Plot of the movement

S_L = filtering(data_healthy.data.LKNE(:,2));
time_L = gate(S_L);

plot_t(data_healthy.data)
plot_t(data_SCI.data)


%% function plot and filter

function plot_t(data)

    %marker = ["LHIP","LKNE", "LANK","LTOE"];
    N = length(data.LHIP(:,1));
    T = 1/120;

    S_L = filtering(data.LTOE(:,2));
    time_L = gate(S_L);

    S_R = filtering(data.RTOE(:,2));
    time_R = gate(S_R);


    %left 

    hip = scatter(data.LHIP(1,2),data.LHIP(1,3),'o','MarkerFaceColor','red');
    hold on 
    knee = scatter(data.LKNE(1,2),data.LKNE(1,3),'o','MarkerFaceColor','blue');
    ankle = scatter(data.LANK(1,2),data.LANK(1,3),'o','MarkerFaceColor','magenta');
    toe = scatter(data.LTOE(1,2),data.LTOE(1,3),'o','MarkerFaceColor','green');
    

    %right

    hip2 = scatter(data.RHIP(1,2),data.RHIP(1,3),'o','red');
    knee2 = scatter(data.RKNE(1,2),data.RKNE(1,3),'o','blue');
    ankle2 = scatter(data.RANK(1,2),data.RANK(1,3),'o','magenta');
    toe2 = scatter(data.RTOE(1,2),data.RTOE(1,3),'o','green');

    

    hold off 

    axis([1 1200 0 1200])
    xlabel('x'), ylabel('y')
    title('kinematic reconstruction')

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

        
        drawnow
        pause(T)

        delete(move_L)
        delete(move_R)
    end

end


function [S_f] = filtering(S)
    S_f = lowpass(highpass(S,1e-1,1e2),0.6,1e2, 'ImpulseResponse','iir');
end
    
function time = gate(S)    
    G = gradient(S);
    time = {""};

    plot(G)

    for i = 2:length(G)
        if sign(G(i)) ~= sign(G(i-1))
            if sign(G(i)) > 0
                time = [time,"foot strike"];
            else
                time = [time,"foot off"];
            end
        else
            time = [time,""];
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