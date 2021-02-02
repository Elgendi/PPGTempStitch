function [peak,onset] = ppg_peak_onset_detection_Hilbert(PPG,Fs,view)
    %% 
    %% cite: A Robust PPG Onset and Systolic Peak Detection Algorithm Based On Hilbert Transform (Abhishek Chakraborty)
    %% 
    
    %%sixth-order Butterworth lowpass filter
    Fc = 15;  %cutoff frequency: 15Hz
    [A,B,C,D] = butter(6,Fc/(Fs/2));
    [filter_SOS,g] = ss2sos(A,B,C,D);    
    PPG_filtered = filtfilt(filter_SOS,g,PPG);
    
    APG = diff(diff(PPG_filtered));  % second derivative
    
    if view == 1 || view == 2
        figure;
        ax(1) = subplot(4,1,1);
        plot(PPG);
        title('Original PPG');
        ax(2) = subplot(4,1,2);
        plot(PPG_filtered);
        title('Filtered PPG (6th butter)');
        ax(3) = subplot(4,1,3);
        plot(APG);
        title('APG');
    end
    %% detect onset and systolic peak
    H = hilbert(APG);
    H = imag(H);
    if view == 1 || view == 2
        ax(4) = subplot(4,1,4);
        plot(H);
        title('APG Hilbert transform');
        linkaxes(ax,'x');
    end
    H_index = find(H>0.5*max(H));  % greater than 50% amplitude
    block = {};
    start  = 1;
    for i = 2:length(H_index)
        if H_index(i)-H_index(i-1) >1
            array = {H_index(start:1:i-1)};
            block = [block;array];
            start = i;
        else
            if i == length(H_index)
                array = {H_index(start:1:length(H_index))};
                block = [block;array];
            end
        end
    end
    num = numel(block);
    window = 0.25*Fs;
    H_peak = zeros(num,1);
    H_left = zeros(num,1);
    H_right = zeros(num,1);
    for i = 1:num
        array1 = H(block{i,1});
        [~,Block_index] = max(array1);
         H_peak(i) = Block_index + block{i,1}(1) -1;  % peak
         %% find the left zero-cross point
         for j = 1:H_peak(i)-1
             if H(H_peak(i) - j) * H(H_peak(i) - j+1) <0
                 H_left(i) = H_peak(i) - j; 
                 break;
             end
         end
         
         
         %% find the right zero-cross point
         for j = 1:length(H)-H_peak(i)-1
             if H(H_peak(i) + j) * H(H_peak(i) + j+1) < 0
                 H_right(i) = H_peak(i) + j; 
                 break;
             end
         end
         
    end
    zero_index = find(H_left == 0);
    H_left(zero_index) = [];
    zero_index = find(H_right == 0);
    H_right(zero_index) = [];
    onset = zeros(length(H_left),1);
    peak = zeros(length(H_right),1);
    all_onset_segment = {};
    all_peak_segment = {};
    for i = 1:length(H_left)
        onset_segment_index = max([1,H_left(i)-14]):1:min([H_left(i)+15,length(PPG_filtered)]);
        for j = 2:length(onset_segment_index)-1
            if PPG_filtered(onset_segment_index(j)) - PPG_filtered(onset_segment_index(j-1)) <0 && ...
                PPG_filtered(onset_segment_index(j)) - PPG_filtered(onset_segment_index(j+1)) <0
                onset(i) = onset_segment_index(j);
            end
        end
        all_onset_segment = [all_onset_segment;onset_segment_index];
    end
    for i = 1:length(H_right)
        peak_segment_index = max([1,H_right(i)-floor(window/2)]):1:min([H_right(i)+ceil(window/2),length(PPG_filtered)]);
        for j = 2:length(peak_segment_index)-1
            if PPG_filtered(peak_segment_index(j)) - PPG_filtered(peak_segment_index(j-1)) >0 && ...
                PPG_filtered(peak_segment_index(j)) - PPG_filtered(peak_segment_index(j+1)) >0
                peak(i) = peak_segment_index(j);
            end
        end
        all_peak_segment = [all_peak_segment;peak_segment_index];
    end
    zero_index = find(onset == 0);
    onset(zero_index) = [];
    zero_index = find(peak == 0);
    peak(zero_index) = [];
    if view == 2 || view == 3 
        figure;
        ax(1) = subplot(4,1,1);
        plot(H);
        hold on;
        plot(H_index,H(H_index),'ro');
        title('Hilbert Transform');
        ax(2) = subplot(4,1,2);
        plot(H);
        hold on;
        plot(H_peak,H(H_peak),'ro');
        plot(H_left,H(H_left),'r+');
        plot(H_right,H(H_right),'r*');
        title('Hilbert Transform');
        ax(3) = subplot(4,1,3);
        plot(PPG_filtered);
        hold on
        for i = 1:length(all_onset_segment)
            plot(all_onset_segment{i,1},PPG_filtered(all_onset_segment{i,1}),'r');
        end
        plot(onset,PPG_filtered(onset),'k*');
        title('Filtered PPG (valley segments)');
        ax(4) = subplot(4,1,4);
        plot(PPG_filtered);
        hold on
        for i = 1:length(all_peak_segment)
            plot(all_peak_segment{i,1},PPG_filtered(all_peak_segment{i,1}),'r');
        end
        plot(peak,PPG_filtered(peak),'k+');
        title('Filtered PPG (peak segments)');
        linkaxes(ax,'x');
        suptitle('Hilbert Method');
    end
    %%
end

