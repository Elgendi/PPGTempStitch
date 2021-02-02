function output = wavefilter(wave,Fs,filterName,FreqRange,order)
        Fn = Fs/2;
        Freq_l = FreqRange(1);
        Freq_h = FreqRange(2);
        
        switch filterName
            case 'butter'
                [A,B,C,D] = butter(order,[Freq_l Freq_h]/Fn);
                [filter_SOS,g] = ss2sos(A,B,C,D);
                output = filtfilt(filter_SOS,g,wave);
%                 freqz(filter_SOS,512,Fs);
            case 'cheby1'
                [A,B,C,D] = cheby1(order,0.1,[Freq_l Freq_h]/Fn);
                [filter_SOS,g] = ss2sos(A,B,C,D);
                output = filtfilt(filter_SOS,g,wave);
            case 'cheby2'
                [A,B,C,D] = cheby2(order,20,[Freq_l Freq_h]/Fn);
                [filter_SOS,g] = ss2sos(A,B,C,D);
%                  freqz(filter_SOS,512,Fs);
                output = filtfilt(filter_SOS,g,wave);
            case 'ellip'
                [A,B,C,D] = ellip(order,0.1,30,[Freq_l Freq_h]/Fn); 
                [filter_SOS,g] = ss2sos(A,B,C,D); 
                output = filtfilt(filter_SOS,g,wave);
            case 'bandpass'
                output = bandpass(wave,FreqRange,Fs);
            case 'smooth'
                output = smooth(wave,order);

            case 'medfilt1'
                output = medfilt1(wave,order);
            case 'wden'
                output= wden(wave,'modwtsqtwolog','s','mln',order,'db2');
        end
end