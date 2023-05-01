classdef emgLib
    methods(Static)
        function smoothed_envelope = filter_emg(emg_signal,plot_signals)
            % INPUTS: 
            %
            %   emg_signal: signal to be filtered
            %   plot_signal: 1 to plot the signals, 0 not to plot
            
            % Set up filter parameters
            cutoff_freq = 50; % Hz
            sampling_freq = 6000; % Hz
            order = 4;
            
            % Design and apply the lowpass filter
            [b, a] = butter(order, cutoff_freq/(sampling_freq/2), 'low');
            filtered_signal = filtfilt(b, a, emg_signal);
            
            % Rectify the filtered signal
            rectified_signal = abs(filtered_signal);
            
            % Compute the envelope using the Hilbert transform
            [u_envelope, ~] = envelope(rectified_signal,500,'peak');
            
            % Smooth the envelope using a moving average filter
            window_size = 1500; % adjust as needed
            smoothed_envelope = movmean(u_envelope, window_size);
            
            % Plot the original signal and the smoothed envelope
            if (plot_signals == 1)
                t = (0:length(emg_signal)-1)/sampling_freq;
                plot(t,emg_signal,t,smoothed_envelope);
                legend('Original Signal','Smoothed Envelope');
                xlabel('Time (s)');
                ylabel('Amplitude');
            end

        end

        function [burst_duration,mean_amplitude,integral,rms] = emg_parameters(emg_signal,plot)
            
            sampling_freq = 6000;

            threshold_on = 0.3 * max(emg_signal);
            threshold_off = 0.1 * max(emg_signal);
            
            % Detect onset and offset of muscle activity
            [onset, offset] = emgLib.detect_bursts(emg_signal, threshold_on, threshold_off);
            
            % Plot the EMG signal, envelope, and detected onset/offset events
            if (plot == 1)
                t = (0:length(emg_signal)-1) / sampling_freq;
                plot(t, emg_signal);
                hold on;
                plot(t, emg_signal, 'r');
                plot(t(onset), emg_signal(onset), 'go');
                plot(t(offset), emg_signal(offset), 'ro');
                xlabel('Time (s)');
                ylabel('Amplitude');
                legend('EMG', 'Envelope', 'Onset', 'Offset');
            end
            
            % Compute burst duration, mean amplitude, integral of muscle activity, and RMS of muscle activity
            burst_duration = (offset - onset) / sampling_freq;
            mean_amplitude = mean(emg_signal(onset:offset));
            integral = trapz(emg_signal(onset:offset)) / sampling_freq;
            rms = sqrt(mean(emg_signal(onset:offset).^2));

        end

        function CI = coactivation_index(antagonist_emg,agonist_emg)
            % Normalize emgs
            antagonist_emg = antagonist_emg/rms(antagonist_emg);
            agonist_emg    = agonist_emg/rms(agonist_emg);
            end_time = length(antagonist_emg);
            CI_vect = zeros(size(antagonist_emg));
            for t = 1:end_time
                lower  = min(abs(antagonist_emg(t)),abs(agonist_emg(t)));
                higher = max(abs(antagonist_emg(t)),abs(agonist_emg(t)));
                CI_vect(t) = lower/higher*(lower+higher);
            end
            CI = trapz(CI_vect)/end_time*100;
        end

        function [onset, offset] = detect_bursts(signal, threshold_on, threshold_off)
            % Detect onset and offset of bursts in an EMG signal
            %
            % INPUTS:
            %   signal: The preprocessed EMG signal
            %   threshold_on: The threshold for detecting the onset of a burst
            %   threshold_off: The threshold for detecting the offset of a burst
            %
            % OUTPUTS:
            %   onset: A vector of onset times (in samples)
            %   offset: A vector of offset times (in samples)
            

            % Find local peaks in the signal
            [peaks, peak_locs] = findpeaks(signal);
            % Sort peaks and peak_locs in ascending order
            [peaks, sort_idx] = sort(peaks);
            peak_locs = peak_locs(sort_idx);
            
            % Initialize onset and offset vectors
            onset = zeros(size(peaks));
            offset = zeros(size(peaks));
            
            % Loop through the peaks to find the corresponding onset and offset
            for i = 1:length(peaks)
                % Find the onset
                j = peak_locs(i);
                while j > 1 && signal(j) > signal(j-1) && signal(j) > threshold_on
                    j = j - 1;
                end
                onset(i) = j;
                
                % Find the offset
                j = peak_locs(i);
                while j < length(signal) && signal(j) > signal(j+1) && signal(j) > threshold_off
                    j = j + 1;
                end
                offset(i) = j;
            end
        
        end

    end
end