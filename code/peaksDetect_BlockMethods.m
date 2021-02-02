function peaks = peaksDetect_BlockMethods(s,Fs,view,F1,F2,W1,W2,W3,beta)
    if nargin <3
        view = 0;
    end
    if nargin <4
        F1 = 0.5;
        F2 = 10;
    end
    if nargin <6
        W1 = floor(Fs*0.111);
        W2 = floor(Fs*0.667);
        W3 = W1;
        beta = 0.02;
    end

   %% W1 the peak 
   %% W2 the beat
    x = wavefilter(s,Fs,'cheby2',[F1, F2],2);
    index = find(x<0);
    x(index) = 0;
    wave = x.*x;
%     wave = x;
    MApeak = smooth(wave, W1);
    MAbeat = smooth(wave, W2);
    
    z = mean(wave);
    alpha = beta * z;
    THR1 = MAbeat + alpha;
    BlocksOfInterest = zeros(1,length(MApeak));
    for n = 1:length(MApeak)
        if MApeak(n) > THR1(n)
            BlocksOfInterest(n) = 1;
        else
            BlocksOfInterest(n) = 0;
        end
    end

    [BlockOnset,BlockEnd] = find_blockes(BlocksOfInterest);
    THR2 = W3;
    peaks = [];
    for j = 1: length(BlockEnd)
            if BlockEnd(j) - BlockOnset(j) + 1 >= THR2
               [~,locs] = max(wave(BlockOnset(j) : BlockEnd(j)));
               locs = locs + BlockOnset(j) -1;
               peaks = [peaks locs];
            end
    end

    if view == 1
        figure;
        plot(wave,'k');
        hold on;
        plot(MApeak,'g');
        plot(MAbeat,'r');
        plot(THR1,'b')
        %             plot(THR1,'k');
        plot(BlocksOfInterest * max(wave),'k--');
        plot(peaks,wave(peaks),'r*');
        legend('wave','MApeak','MAbeat','THR1');
        title(['Window=' num2str(W1)]);
    end
end

function [BlockOnset,BlockEnd] = find_blockes(BlocksOfInterest)
    BlockOnset = zeros(1,3000);
    BlockEnd = zeros(1,3000);
    count = 1;
    start_locs = 1;
    for i = 2:1:length(BlocksOfInterest)
        if BlocksOfInterest(i-1) == 0 && BlocksOfInterest(i) == 1
            start_locs = i;
            break;
        end
    end
    if start_locs ~= 1  
        for i = start_locs:1:length(BlocksOfInterest)-2
            if BlocksOfInterest(i-1) == 0 && BlocksOfInterest(i) == 1
                BlockOnset(count) = i;
            elseif BlocksOfInterest(i) == 1 && BlocksOfInterest(i+1) == 0
                BlockEnd(count) = i;
                count = count + 1;
            end
        end
    end
    BlockOnset(BlockEnd == 0) = [];
    BlockEnd(BlockEnd == 0) = [];
end