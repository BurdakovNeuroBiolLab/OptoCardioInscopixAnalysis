function alpha1 = dfa_alpha1_local(x)
% DFA_ALPHA1_LOCAL: short-term DFA exponent for time series x (IBI).
% Box sizes run from 4 beats up to floor(N/4), capped at 16 beats (the
% standard alpha1 range), so a long series still gives alpha1 rather than
% a mix of short- and long-term scaling. Needs >= 8 values.

    x = x(:);
    N = numel(x);

    % Center and integrate (profile)
    x = x - mean(x);
    y = cumsum(x);

    % Choose box sizes based on N
    % Minimum box size = 4, maximum = floor(N/4), at most 16
    nMin = 4;
    nMax = min(floor(N/4), 16);

    nVals = nMin:nMax;
    F = nan(size(nVals));

    for i = 1:numel(nVals)
        n = nVals(i);
        nSeg = floor(N / n);
        if nSeg < 2
            F(i) = NaN;
            continue;
        end

        rmsSeg = zeros(nSeg,1);
        for s = 1:nSeg
            idx = (s-1)*n + (1:n);
            t = (1:n)';
            p = polyfit(t, y(idx), 1);    % linear trend
            yfit = polyval(p, t);
            rmsSeg(s) = sqrt(mean((y(idx) - yfit).^2));
        end
        F(i) = sqrt(mean(rmsSeg.^2));
    end

    valid = ~isnan(F) & F > 0;
    if sum(valid) < 2
        alpha1 = NaN;
        return;
    end

    logn = log(nVals(valid));
    logF = log(F(valid));

    coeffs = polyfit(logn(:), logF(:), 1);
    alpha1 = coeffs(1);
end
