function factor = normalise_loudness (x, fs, lkfs)
%LOUDNESS_MATCH
%   LOUDNESS_MATCH (X, fs, lkfs) calculates an amplitude multiplication factor
%   which equalizes a sound to match the specified loudness (in LKFS).
%
%   2010-02-23 by MARUI Atsushi
%   http://www.geidai.ac.jp/~marui/matlab/node43.html
%
%   Modified by Nic 04/03/2016

factor = 1.0;
factorHigh = 10^(+60/20);
factorLow = 10^(-60/20);
i = 1;

while i
    s = loudness_itu(factor * x, fs);
    k = (s - lkfs) / lkfs;
    
%     fprintf('%3d:   Factor %7.5f   LKFS %7.3f   (%6.2f%% difference)\n', ...
%         i, factor, s, k*100);
    i = i + 1;
    
    if abs(k) < 0.001
        return;
    end
    
    if k < 0.0
        factorOld = factor;
        factor = (factorLow + factor) / 2.0;
        factorHigh = factorOld;
    elseif k > 0.0
        factorOld = factor;
        factor = (factorHigh + factor) / 2.0;
        factorLow = factorOld;
    end
end

function lkfs = loudness_itu(x, fs)
%LOUDNESS_ITU   compute loudness (LKFS) based on ITU-R BS.1770-2
%   LOUDNESS_ITU(x, fs, mode) compute loudness based on ITU-R BS.1770-2
%   specification.  x is an input signal (mono/stereo/5.1ch) with sampling
%   frequency fs.
%
%   Input signal has to be either in mono (M), stereo (L and R), or
%   5.1ch(L, R, C, LFE, Ls, and Rs) in respective channel order.
%
%   Example: loudness calculation
%   [x,fs] = wavread('soundfile.wav');
%   lkfs = LOUDNESS_ITU(x, fs);
%
%   2012-01-06 MARUI Atsushi
%   http://www.geidai.ac.jp/~marui/matlab/node42.html

% constants
blockSize = 400 / 1000 * fs;   % in samples
overlapSize = 0.75;            % in ratio
channelWeights = [1 1 1 0 sqrt(2) sqrt(2)];   % in L, R, C, LFE, Ls, Rs order
absoluteThreshold = -70;   % in dB
relativeThreshold = -10;   % in dB

% preparation
if size(x,1)~=length(x)
    x = x';
end
numch = size(x,2);

if fs~=48000
    x = resample(x, 48000, fs);
    fs = 48000;
end

switch(numch)
    case 5
        chwat = channelWeights([1 2 3 5 6]);
    otherwise
        chwat = channelWeights(1:numch);
end

% K-filter
B1 = [
    1.53512485958697
    -2.69169618940638
    1.19839281085285
    ];
A1 = [
    1.0
    -1.69065929318241
    0.73248077421585
    ];
B2 = [
    1.0
    -2.0
    1.0
    ];
A2 = [
    1.0
    -1.99004745483398
    0.99007225036621
    ];
y = filter(B2,A2,filter(B1,A1,x));

% Mean square
numBlock = ceil(length(y)/ blockSize) + 1;
yy = zeros(numBlock * blockSize, numch);
yy(1:length(y),:) = y;
j = 0:(length(y) - blockSize)/(blockSize * (1-overlapSize));
z = zeros(length(j), numch);
for n1 = 1:length(j)
    yyy = yy(blockSize*j(n1)*(1-overlapSize)+1:blockSize*(j(n1)*(1-overlapSize)+1)+1,:);
    z(n1,:) = sum(yyy .^ 2) / blockSize;
end
l = -0.691 + 10*log10(sum(repmat(chwat,length(z),1) .* z, 2));

% Gating (absolute)
Jg = l > absoluteThreshold;
zz = repmat(channelWeights(1:numch),length(z),1) .* z;
L_KG = -0.691 + 10*log10(chwat * sum(zz(Jg,:), 1)' / sum(Jg));

% Gating (relative)
Gamma_r = L_KG + relativeThreshold;
Jg = l > Gamma_r;
L_KG = -0.691 + 10*log10(chwat * sum(zz(Jg,:), 1)' / sum(Jg));

% OPTIONAL: draw pretty figure
if false
    t = (blockSize*(j*(1-overlapSize)+1)+1)/fs;
    h = plot(t, l, 'color', [.5 .5 .5]);
    hold on;
    h = line([t(1) t(end)], [L_KG L_KG]);
    set(h, 'Color', 'b');
    set(h, 'LineStyle', '--');
    set(h, 'LineWidth', 1);
    h = line([t(1) t(end)], [Gamma_r Gamma_r]);
    set(h, 'Color', 'b');
    set(h, 'LineStyle', ':');
    set(gca, 'Xlim', [t(1) t(end)]);
    hold off;
    legend({'momentary', 'integrated', 'relative gate'});
end

% output
lkfs = L_KG;