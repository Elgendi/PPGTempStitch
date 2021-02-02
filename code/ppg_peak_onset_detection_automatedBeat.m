function [peak,onset,ppg] = ppg_peak_onset_detection_automatedBeat(PPG,Fs,View)
    % cite: automated Beat Onset and Peak Detection Algorithm for Field-Collected Photoplethysmograms (Liangyou Chen)
    %
    
    
    %% find outlier
    missing_point = find(ismissing(PPG) == 1);
    outlier = find(PPG>median(PPG)*20);
    outlier = [missing_point outlier];
    t = 1:length(PPG);
    if ~isempty(outlier)
        w0 = interp1(t,PPG,outlier);
    else
        w0 = PPG;
    end
    %% find main frequency
    waveLength = length(w0);
    y = w0 - mean(w0);

    NFFT = max(256,2^nextpow2(waveLength));
    [pxx,f] = pwelch(y,waveLength,waveLength/2,(NFFT*2)-1,Fs);
    beatFreqRange = f>=0.8 & f<=3;    
    f1 =  f(beatFreqRange);   
    pxx1 = pxx(beatFreqRange);
    [maxPower,MaxIndex] = max(pxx1);
    heartRate = f1(MaxIndex);
    windowSize = round(0.2*(1/heartRate)*Fs);
    w1 = medfilt1(w0,windowSize);
    w1 = smooth(w1,windowSize);
    
    Fc = 1.5*f1(MaxIndex);  %cutoff frequency: 1.5 * heart rate frequency
    [A,B,C,D] = butter(3,Fc/(Fs/2));
    [filter_SOS,g] = ss2sos(A,B,C,D); 
    w2 = filtfilt(filter_SOS,g,w1);
    
    b = smooth(w2,1.5*Fs/heartRate);
    
    s = find(w2>b);
    start = 1;
    blocks = {};
    for i =2:length(s)
        if s(i) - s(i-1) >1
            array = s(start:1:i-1);
            blocks = [blocks;array];
            start = i;
        end
    end
    array = s(start:1:i);
    blocks = [blocks;array];
    
    for i = 1:length(blocks)
        [~,index] = max(w0(blocks{i,1}));
        initialPeaks(i) = index + blocks{i,1}(1) - 1;
    end
    initialPeaksValue = w0(initialPeaks);
    initialPeaksValue_sorted =sort(initialPeaksValue);
    heightValue = initialPeaksValue_sorted(round(length(initialPeaksValue_sorted)*2/3));
    
    falsepeak1 = find(initialPeaksValue<0.5*heightValue);
    t = diff(initialPeaks);
    MAD = median(abs(t-median(t)));
    falsepeak2 = find(abs(t-median(t))>=2*MAD);
    falseIndex = unique([falsepeak1,falsepeak2]);
    
    Peaks_deleteFalse = initialPeaks;
    Peaks_deleteFalse(falseIndex) = [];
    isDeletedPeaks = zeros(numel(Peaks_deleteFalse),1);
    start = 1;
    for i = 1:length(Peaks_deleteFalse)
        num = find(initialPeaks(falseIndex)> start & initialPeaks(falseIndex)<Peaks_deleteFalse(i));
        start = Peaks_deleteFalse(i);
        if ~isempty(num)
            isDeletedPeaks(i) = 1;
        end
        
    end
    
    %%Relocation of missed peaks
    missedPeaks = [];
    for i = 1:length(isDeletedPeaks)
        if i == 233
            ffff = 1;
        end
        if isDeletedPeaks(i) == 1
            if i == 1
                start = 1;
            else
                start = Peaks_deleteFalse(i-1);
            end
            n = 1;
            
            endIndex = start+ ceil(median(t));
            while(endIndex < Peaks_deleteFalse(i) )
                
                segment = w1(max([1,start+ceil(median(t))-ceil(MAD*n/2)]): min([start+ceil(median(t))+ceil(MAD*n/2),length(w1)]));
                
                [~,segmentMaxIndex] = max(segment);
                if ~(segmentMaxIndex == length(segment) || segmentMaxIndex == 1)
                    [~,w0MaxIndex] = max(w0(max([1,start+ceil(median(t))-ceil(MAD*n/2)]): min([start+ceil(median(t))+ceil(MAD*n/2),length(w1)])));
                    missedPeaks = [missedPeaks, w0MaxIndex+start+ceil(median(t))-ceil(MAD*n/2)-1];
                    start = segmentMaxIndex+start+ ceil(median(t)) -ceil(MAD*n/2)-1;
                    n = 0;
                end
                endIndex = start+ ceil(median(t))+ceil(MAD*n/2);
                n = n+1;
            end
            
        end
    end
    truePeaks = sort([missedPeaks Peaks_deleteFalse]);
    %%onset detection
    start = 1;
    onset = [];
    for i = 1:length(truePeaks)
        %step 1: Found ranges where w0 was below both the w2 and b
        w0_1 = w0(start:truePeaks(i));
        w2_1 = w2(start:truePeaks(i));
        b_1 = b(start:truePeaks(i));
        rangeIndex = find(w0_1(:)<w2_1(:) & w0_1(:)<b_1(:));
        %step 2: If there were multiple ranges identified, we ranked them based on their lengths, 
        % and selected the rightmost one from the top two ranges. 
        onset_blocks_index = find(diff(rangeIndex) >1);

        if numel(onset_blocks_index) + 1>1
            onsetBlock = {};
            len = [];
            for j = 1:numel(onset_blocks_index)
                if j == 1
                    onsetBlock = [onsetBlock;rangeIndex(1:onset_blocks_index(1))];
                else
                    onsetBlock = [onsetBlock;rangeIndex(onset_blocks_index(j-1):onset_blocks_index(j))];
                end
                len = [len;numel(onsetBlock{j,1})];
            end
            onsetBlock = [onsetBlock;rangeIndex(onset_blocks_index(j)+1:end)];
            len = [len;numel(onsetBlock{j+1,1})];
            [len_sort,I] = sort(len,'descend');
            I_max = max(I(1:2));
            onsetBlock = onsetBlock{I_max,1};
        else
            onsetBlock = rangeIndex;
        end
        %step 3: Identified the minimum position on the waveform w0 in the selected range as the onset location
        if ~isempty(onsetBlock)
            [~,I_onset] = min(w0(onsetBlock + start - 1)); 
            onset_1 = I_onset + onsetBlock(1) -1 + start - 1;
            onset = [onset,onset_1];
        end
        start = truePeaks(i);
    end
    
    peak = truePeaks;
    ppg = w0;
    if View == 1
        x = (1:length(ppg))/Fs;
        figure;
        subplot(3,2,1);
        plot(x,w0);
        title('w0');
        subplot(3,2,2);
        plot(f1,pxx1);
        title('Power Spectrum');
        subplot(3,2,3);
        plot(x,w0);
        hold on;
        plot(x,w1);
        plot(x,w2);
        plot(x,b);
        legend('w0','w1','w2','b');
        subplot(3,2,4);
        plot(x,w0);
        hold on;
        plot(x,w2);
        plot(x,b);
        plot(x(initialPeaks),w0(initialPeaks),'b+');
        plot(x(initialPeaks(falseIndex)),w0(initialPeaks(falseIndex)),'r*');
        subplot(3,2,5);
        plot(x,w0);
        hold on;
        plot(x,w2);
        plot(x,b);
        plot(x(initialPeaks),w0(initialPeaks),'b+');
        plot(x(initialPeaks(falseIndex)),w0(initialPeaks(falseIndex)),'r*');
        plot(x(truePeaks),w0(truePeaks),'bs');
        subplot(3,2,6);
        plot(x,w0);
        hold on;
        plot(x,w2);
        plot(x,b);
        plot(x(truePeaks),w0(truePeaks),'bs');
        plot(x(onset),w0(onset),'bd');
    end
end