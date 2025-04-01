function APD90 = APD90(t,v)
    t=t';
    dVdt = diff(v)./diff(t);
    VdotMax = max(dVdt);
    i0 = find(dVdt == VdotMax) + 1;
    t0 = t(i0);

    %essentially gather the index for the closest time point -+5 away
    [~, startIndex] = min(abs(t - (t0 - 5)));
    [~, endIndex] = min(abs(t - (t0 + 5)));

    minDiastolic = min(v(startIndex:endIndex));
    peak = max(v(startIndex:endIndex));
    
    %APD90 magnitude must be calculated with lower bound taken into
    %account
    endPeak = 0.10*(peak+abs(minDiastolic))-abs(minDiastolic);
    vBottom = max(find(abs(v - endPeak) < 0.7));
    
    tf = t(vBottom);
    
    APD90 = tf-t0;
end