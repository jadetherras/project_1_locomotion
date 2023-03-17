
%dataset sain 
data=load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/1_AML01_2kmh.mat");

% dataset SCI Human
%data=load("SCI Human/DM002_TDM_08_1kmh.mat");

% Plot of the movement

[N,t] = size(data.data.LHIP);
 
t = 1:N;
t = t/120;

%% with function

plot_t(data.data)

%% filter and calculate pos and calculate gate

function plot_t(data)

    marker = ["LHIP","LKNE", "LANK","LTOE"];
    N = length(data.LHIP(:,1));
    T = 1/120;

%figure 
%     legend(marker)
%     xlabel('y'), ylabel('z')
%     title('plot marker')
    
    figure
    
    hip = scatter(data.LHIP(1,2),data.LHIP(1,3),'o','MarkerFaceColor','red');
    hold on 
    knee = scatter(data.LKNE(1,2),data.LKNE(1,3),'o','MarkerFaceColor','blue');
    ankle = scatter(data.LANK(1,2),data.LANK(1,3),'o','MarkerFaceColor','yellow');
    toe = scatter(data.LTOE(1,2),data.LTOE(1,3),'o','MarkerFaceColor','green');
    hold off 

    axis([1 1200 0 1200])

    for n = 2:N
        hip.XData = data.LHIP(n,2);
        hip.YData = data.LHIP(n,3);
        knee.XData = data.LKNE(n,2);
        knee.YData = data.LKNE(n,3);
        ankle.XData = data.LANK(n,2);
        ankle.YData = data.LANK(n,3);
        toe.XData = data.LTOE(n,2);
        toe.YData = data.LTOE(n,3);
        drawnow
        pause(T)
    end


end

function [S_f] = filtering(S)
    S_f = lowpass(highpass(S,1e-1,1e2),0.6,1e2, 'ImpulseResponse','iir');
end
    
function [time,col] = gate(S)    
    G = gradient(S);
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