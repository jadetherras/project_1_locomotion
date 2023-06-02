
%% Kinematic file

% In this file, we visualise the kinematic data and the gate detection 

% filter and plot the signal
% cut in gate cycle
% animation of the movement
% plot the gate cycle

%% Loading the data

close all;

%dataset sain 

%choose a dataset
%data=load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/1_AML01_2kmh.mat"); 
%data=load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/1_AML02_2kmh.mat");
%data=load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/3_AML01_1kmh.mat");
%data=load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/3_AML02_1kmh.mat");
data =load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/4_AML01_3kmh.mat");
%data=load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/4_AML02_3kmh.mat");
%data =load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/2_AML01_3kmh_inclined.mat");
%data =load("Healthy dataset (CHUV recording - 03.03.2023)-20230310/2_AML02_3kmh_inclined.mat");


% dataset SCI Human
%data=load("SCI Human/DM002_TDM_08_1kmh.mat");
%data=load("SCI Human/DM002_TDM_08_2kmh.mat");
%data=load("SCI Human/DM002_TDM_1kmh_NoEES.mat");

%provide a name for the trail
name = ' 1km healthy';

%% Plot the movement


%test range

%if you want all the data
N = length(data.data.LHIP(:,1));

%if you want just one gait
start = 1; %this values come from the gate detection, we decided to plot one gait randomly
stop = 1500;

T = 1/data.data.marker_sr;
dec = 1000/3600*T*1000; %change here the velocity
%dec = 0;

%visualisation 2 limb
%double_dec = 1000;
double_dec = 0;


%filter the data (here toe) and gate calculation
%S_L = filtering(data_healthy.data.LTOE(:,2));
%time_L = gate(S_L);

%if only cut the gate
%Gate = cut_gate(data.data,true,true,name);

%plot the gate cycle
%plot_gate(data.data,start,stop,'R',1,dec,double_dec,name)

%animate the data for comparison with visualisation
animate(data.data, start, stop,dec,name) 



%% function

% filtering : take a signal (ex the y value for a marker) and return a
% clean signa
function [S_f] = filtering(S)

    %low, high pass filter
    d1 = designfilt("lowpassiir",FilterOrder=2, HalfPowerFrequency=0.03,DesignMethod="butter");
    S_f = filtfilt(d1,highpass(S,1e-1,1e2));
end
   

% cut_gate : cut the date in gate cycle (using a gradient)
% each value represent a foot strike
% data is the dataset
% gate is a array of position 
function Gate = cut_gate(data,constraint,plotting,name) 
    
    T = 1/data.marker_sr;

    %filtering
    Sy_L = filtering(data.LANK(:,2));
    Sy_R = filtering(data.RANK(:,2));
    Sz_L = filtering(data.LANK(:,3));
    Sz_R = filtering(data.RANK(:,3));

     % calculate the gradient
    Gyl = normalize(gradient(Sy_L));
    Gyr = normalize(gradient(Sy_R));
    Gzl = normalize(gradient(Sz_L));
    Gzr = normalize(gradient(Sz_R));

    if plotting 
         figure 
         plot(Sy_L)
         hold on 
         plot(Sy_R)
         plot(Sz_L)
         plot(Sz_R)
         hold off
         title(strcat('function',name))
         legend('Sy_L','Sy_R','Sz_L','Sz_R')
         saveas(gcf,strcat('figure/function ',name,'.png'))

        figure 
        plot(Gyl)
        hold on 
        plot(Gyr)
        plot(Gzl)
        plot(Gzr)
        hold off
        title(strcat('derivative ',name))
        legend('Gy_L','Gy_R','Gz_L','Gz_R')
        saveas(gcf,strcat('figure/derivative ',name,'.png'))
     end


    %initialisation
    Gate = [];
    gate = [];
    
    on_gate = false;

    % calculate the position
    i = 1;
    while i <=(length(Gyl)-1)
        i = i+1;
        %start of the forward movement -> foot off
        if sign(Gyl(i)) ~= sign(Gyl(i-1)) && sign(Gyl(i)) < 0 
            
             while abs(Gzl(i)) > 1 && i <=(length(Gyl)-1)
                 i = i+1;
             end

            if isfield(gate,'strikeR')
                
                gate.offnext = i;
                gate.duration = gate.offnext - gate.offL;
                gate.stepL = Sy_L(gate.offL)-Sy_L(gate.strikeL);
                gate.stepR = Sy_R(gate.offR)-Sy_R(gate.strikeR);
                
                % we don't take the speed in account, just arbitrary for
                % step but will remore every oscillation of the legs (not
                % proper step), also remove too long gate (when we probably
                % miss a event)
                if not(constraint) || (constraint && gate.duration*T <=4 && gate.duration*T >0.5 && gate.stepL>80 && gate.stepR>80)
                    Gate = [Gate,gate];
                elseif constraint  
                    disp(['remove a gate with duration',num2str(gate.duration),'in position',num2str(gate.offL)])
                end
            end

            on_gate = true;
            gate = [];
            gate.offL = i;


            while on_gate && i <=(length(Gyl)-1)
                i = i+1;

                %end of forward movement
                if sign(Gyl(i-1)) ~= sign(Gyl(i)) && sign(Gyl(i)) > 0 && isfield(gate,'offL')
                      while abs(Gzl(i)) > 1 && i <=(length(Gyl)-1)
                          i = i+1;
                      end
                     gate.strikeL = i;

                %start other limb forward movement
                elseif sign(Gyr(i-1)) ~= sign(Gyr(i)) && sign(Gyr(i)) < 0 && isfield(gate,'strikeL')
                     while abs(Gzr(i)) > 1 && i <=(length(Gyl)-1)
                         i = i+1;
                     end
                    gate.offR = i;

                %end of other limb forward movement
                elseif sign(Gyr(i-1)) ~= sign(Gyr(i)) && sign(Gyr(i)) > 0 && isfield(gate,'offR')
                    while abs(Gzr(i)) > 1 && i <=(length(Gyl)-1)
                        i = i+1;
                    end
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
function plot_gate(data,N1,N2,side,frame, dec, double_dec,name)
    
    Gate = cut_gate(data,true,false,name);
    Min = ceil(N1/frame);
    Max = ceil(N2/frame);

    color_left = 'black';
    color_right = 'green';
    color_gateoff = 'red';
    color_gatestrike = 'blue';


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

    xlabel('x'), ylabel('y')
    title(strcat('kinematic reconstruction for ',name))
    saveas(gcf,strcat('figure/gate ',name,'.png'))

end

% animate : plot a animation of a dataset
% the left legs is represented by full circular markers
% the left legs is represented by empty circular markers
% on the corner, a text display "foot strike/foot off" at the correct time
% the toe is use to cal
% culate the gate
% data is the dataset
function animate(data, start, stop,dec,name)

    %marker = ["LHIP","LKNE", "LANK","LTOE"];
    T = 1/(data.marker_sr);

    % calculate timing left and right

    Gate = cut_gate(data,true,false,name);
    
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