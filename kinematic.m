
%% Kinematic file

% In this file, we visualise the kinematic data and the gate detection 

% filter and plot the signal
% cut in gate cycle
% animation of the movement
% plot the gate cycle

%% Loading the data

close all;

%dataset sain 
%data=load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/3_AML01_1kmh.mat");
%data=load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/4_AML02_3kmh.mat");
%data=load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/3_AML02_1kmh.mat");


% dataset SCI Human
data=load("SCI Human/DM002_TDM_08_1kmh.mat");
%data=load("SCI Human/DM002_TDM_1kmh_NoEES.mat");

%% Plot the movement

% animation : decomment the one you want (healthy, SCI)

%test range

%if you want all the data
N = length(data.data.LHIP(:,1));

%if you want just one gait
start = 1; %this values come from the gate detection, we decided to plot one gait randomly
stop = N;

T = 1/data.data.marker_sr;
%dec = 1000/3600*T*1000; %change here the velocity
dec = 0;

%visualisation 2 limb
double_dec = 1000;
%double_dec = 0;


%filter the data (here toe) and gate calculation
%S_L = filtering(data_healthy.data.LTOE(:,2));
%time_L = gate(S_L);

Gate = cut_gate(data.data,true);

%plot the gate cycle
%plot_gate(data_healthy.data,start,stop,'B',1,dec,double_dec)

%animate the data for comparison with visualisation

animate(data.data, start, stop,dec) 
%animate(data_SCI.data)




%% function

% filtering : take a signal (ex the y value for a marker) and return a
% clean signa
function [S_f] = filtering(S)

    %low, high pass filter
    d1 = designfilt("lowpassiir",FilterOrder=2, HalfPowerFrequency=0.03,DesignMethod="butter");
    S_f = filtfilt(d1,S);
end
    
% gate : take a signal and calculate the gate cycle (using the gradient)
% time take the value foot off and foot strike at the correct timing
function time = gate(S) 

    % calculate the gradient
    G = gradient(S);
    
    %initialisation
    time = {''};

    % calculate the position
    for i = 2:length(G)
        if sign(G(i)) ~= sign(G(i-1)) % if change sign, pass 0
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

% cut_gate : cut the date in gate cycle (using a gradient)
% each value represent a foot strike
% data is the dataset
% gate is a array of position 
function Gate = cut_gate(data,constraint) 
    
    T = 1/data.marker_sr;
    %filtering
     S_L = filtering(data.LANK(:,2));
     S_R = filtering(data.RANK(:,2));

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
                gate.duration = gate.offnext - gate.offL;
                gate.stepL = S_L(gate.offL)-S_L(gate.strikeL);
                gate.stepR = S_R(gate.offR)-S_R(gate.strikeR);
                
                % we don't take the speed in account, bad but anyway it
                % does work because remove only in SCI (slow)
                if not(constraint) || (constraint && gate.duration*T <=4 && gate.stepL>50 && gate.stepR>50)
                    Gate = [Gate,gate];
                elseif constraint  
                    disp(['remove a gate with duration',num2str(gate.duration),'in position',num2str(gate.offL)])
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

% plot_gate : take the dataset, a starting time and ending time (ex
% start/end of a gate cycle)
% plot the gate as a line for each time 
% data is the dataset
% N1 is the start position
% N2 is the end position
% side is a string giving the limb to plot => L/R/B (left, right,both)
% frame say how much to display (if 1 all the timepoint, if 5 every 5
% delta_time for exemple)
% dec permit to introduce a decalage in the y direction. We walk on a
% treadmill so all the gate appear on the same position. If we want to
% introduce a movement in y or visualise more gate cycle we can introduce a
% decalage, like a velocity. if 0 all fixed. 
function plot_gate(data,N1,N2,side,frame, dec, double_dec)

    xlabel('x'), ylabel('y')
    title('kinematic reconstruction')

    Min = ceil(N1/frame);
    Max = ceil(N2/frame);

    color_left = 'black';
    color_right = 'green';
    color_gateoff = 'red';
    color_gatestrike = 'blue';

    Gate = cut_gate(data,true);

    V = 1;
    while V < length(Gate) && Gate(V).offL < Min
        V = V+1;
    end

    for Val = Min:Max
        
        first = true;
        n = Val*frame;

        if n >= Gate(V).offnext && V < length(Gate)
            V = V+1;
        end

        if side == 'L' || side == 'B'
            
            if n == Gate(V).offL 
                plot([data.LHIP(n,2)-dec*(n-Min*frame) data.LKNE(n,2)-dec*(n-Min*frame)],[data.LHIP(n,3)+double_dec data.LKNE(n,3)+double_dec],color_gateoff)
            elseif n == Gate(V).strikeL
                plot([data.LHIP(n,2)-dec*(n-Min*frame) data.LKNE(n,2)-dec*(n-Min*frame)],[data.LHIP(n,3)+double_dec data.LKNE(n,3)+double_dec],color_gatestrike)
            else 
                plot([data.LHIP(n,2)-dec*(n-Min*frame) data.LKNE(n,2)-dec*(n-Min*frame)],[data.LHIP(n,3)+double_dec data.LKNE(n,3)+double_dec],color_left)
            end

            if first
                hold on
                first = false;
            end

            if n == Gate(V).offL 
                plot([data.LKNE(n,2)-dec*(n-Min*frame) data.LANK(n,2)-dec*(n-Min*frame)],[data.LKNE(n,3)+double_dec data.LANK(n,3)+double_dec],color_gateoff)
                plot([data.LANK(n,2)-dec*(n-Min*frame) data.LTOE(n,2)-dec*(n-Min*frame)],[data.LANK(n,3)+double_dec data.LTOE(n,3)+double_dec],color_gateoff)
            elseif n == Gate(V).strikeL
                plot([data.LKNE(n,2)-dec*(n-Min*frame) data.LANK(n,2)-dec*(n-Min*frame)],[data.LKNE(n,3)+double_dec data.LANK(n,3)+double_dec],color_gatestrike)
                plot([data.LANK(n,2)-dec*(n-Min*frame) data.LTOE(n,2)-dec*(n-Min*frame)],[data.LANK(n,3)+double_dec data.LTOE(n,3)+double_dec],color_gatestrike)
            else 
                plot([data.LKNE(n,2)-dec*(n-Min*frame) data.LANK(n,2)-dec*(n-Min*frame)],[data.LKNE(n,3)+double_dec data.LANK(n,3)+double_dec],color_left)
                plot([data.LANK(n,2)-dec*(n-Min*frame) data.LTOE(n,2)-dec*(n-Min*frame)],[data.LANK(n,3)+double_dec data.LTOE(n,3)+double_dec],color_left)
            end
            
        end

        if side == 'R' || side == 'B'
            
            if n == Gate(V).offR
                plot([data.RHIP(n,2)-dec*(n-Min*frame) data.RKNE(n,2)-dec*(n-Min*frame)],[data.RHIP(n,3) data.RKNE(n,3)],color_gateoff)
            elseif n == Gate(V).strikeR
                plot([data.RHIP(n,2)-dec*(n-Min*frame) data.RKNE(n,2)-dec*(n-Min*frame)],[data.RHIP(n,3) data.RKNE(n,3)],color_gatestrike)
            else 
                plot([data.RHIP(n,2)-dec*(n-Min*frame) data.RKNE(n,2)-dec*(n-Min*frame)],[data.RHIP(n,3) data.RKNE(n,3)],color_right)
            end

            if first && side == 'R'
                hold on
                first = false;
            end

            if n == Gate(V).offR 
                plot([data.RKNE(n,2)-dec*(n-Min*frame) data.RANK(n,2)-dec*(n-Min*frame)],[data.RKNE(n,3) data.RANK(n,3)],color_gateoff)
                plot([data.RANK(n,2)-dec*(n-Min*frame) data.RTOE(n,2)-dec*(n-Min*frame)],[data.RANK(n,3) data.RTOE(n,3)],color_gateoff)
            elseif n == Gate(V).strikeR
                plot([data.RKNE(n,2)-dec*(n-Min*frame) data.RANK(n,2)-dec*(n-Min*frame)],[data.RKNE(n,3) data.RANK(n,3)],color_gatestrike)
                plot([data.RANK(n,2)-dec*(n-Min*frame) data.RTOE(n,2)-dec*(n-Min*frame)],[data.RANK(n,3) data.RTOE(n,3)],color_gatestrike)
            else 
                plot([data.RKNE(n,2)-dec*(n-Min*frame) data.RANK(n,2)-dec*(n-Min*frame)],[data.RKNE(n,3) data.RANK(n,3)],color_right)
                plot([data.RANK(n,2)-dec*(n-Min*frame) data.RTOE(n,2)-dec*(n-Min*frame)],[data.RANK(n,3) data.RTOE(n,3)],color_right)
            end

        end


    end

    hold off
end

% animate : plot a animation of a dataset
% the left legs is represented by full circular markers
% the left legs is represented by empty circular markers
% on the corner, a text display "foot strike/foot off" at the correct time
% the toe is use to cal
% culate the gate
% data is the dataset
function animate(data, start, stop,dec)

    %marker = ["LHIP","LKNE", "LANK","LTOE"];
    T = 1/(data.marker_sr*3);

    % calculate timing left and right

    Gate = cut_gate(data,true);
    
    V = 1;
    if start >= Gate(V).offnext && V < length(Gate)
        V = V+1;
    end

    %left plot
    hip = scatter(data.LHIP(start,2),data.LHIP(start,3),'o','MarkerFaceColor','red');
    hold on 
    knee = scatter(data.LKNE(start,2),data.LKNE(start,3),'o','MarkerFaceColor','blue');
    ankle = scatter(data.LANK(start,2),data.LANK(start,3),'o','MarkerFaceColor','magenta');
    toe = scatter(data.LTOE(start,2),data.LTOE(start,3),'o','MarkerFaceColor','green');
    

    %right plot
    hip2 = scatter(data.RHIP(start,2),data.RHIP(start,3),'o','red');
    knee2 = scatter(data.RKNE(start,2),data.RKNE(start,3),'o','blue');
    ankle2 = scatter(data.RANK(start,2),data.RANK(start,3),'o','magenta');
    toe2 = scatter(data.RTOE(start,2),data.RTOE(start,3),'o','green');

    hold off 

    axis([-1000 1000 100 1200])
    xlabel('x'), ylabel('y')
    title('kinematic reconstruction')
    
    change = false;
    %update the markers and text at each time
    for n = (start+1):stop
        
      
        
        hip.XData = data.LHIP(n,2)-dec*(n-start);
        hip.YData = data.LHIP(n,3);
        knee.XData = data.LKNE(n,2)-dec*(n-start);
        knee.YData = data.LKNE(n,3);
        ankle.XData = data.LANK(n,2)-dec*(n-start);
        ankle.YData = data.LANK(n,3);
        toe.XData = data.LTOE(n,2)-dec*(n-start);
        toe.YData = data.LTOE(n,3);

        hip2.XData = data.RHIP(n,2)-dec*(n-start);
        hip2.YData = data.RHIP(n,3);
        knee2.XData = data.RKNE(n,2)-dec*(n-start);
        knee2.YData = data.RKNE(n,3);
        ankle2.XData = data.RANK(n,2)-dec*(n-start);
        ankle2.YData = data.RANK(n,3);
        toe2.XData = data.RTOE(n,2)-dec*(n-start);
        toe2.YData = data.RTOE(n,3);
        
         %normal color

        if change 
            hip.MarkerFaceColor = 'red';
            knee.MarkerFaceColor = 'blue';
            ankle.MarkerFaceColor = 'magenta';
            toe.MarkerFaceColor = 'green';
            hip2.CData = [1,0,0];
            knee2.CData = [0,0,1];
            ankle2.CData = [1,0,1];
            toe2.CData = [0,1,0];
            change = false;
        end


        if n == Gate(V).offL 
            hip.MarkerFaceColor = 'red';
            knee.MarkerFaceColor = 'red';
            ankle.MarkerFaceColor = 'red';
            toe.MarkerFaceColor = 'red';
            change = true;
        elseif n == Gate(V).offR 
            hip2.CData = [1,0,0];
            knee2.CData = [1,0,0];
            ankle2.CData = [1,0,0];
            toe2.CData = [1,0,0];
            change = true;
        elseif n == Gate(V).strikeL
            hip.MarkerFaceColor = 'blue';
            knee.MarkerFaceColor = 'blue';
            ankle.MarkerFaceColor = 'blue';
            toe.MarkerFaceColor = 'blue';
            change = true;
        elseif n == Gate(V).strikeR 
            hip2.CData = [0,0,1];
            knee2.CData = [0,0,1];
            ankle2.CData = [0,0,1];
            toe2.CData = [0,0,1];
            change = true;
            V = V+1;
        end

        % give the time to see before change
        drawnow
        if change
            pause(1)
        else
            pause(T)
        end
    end

end