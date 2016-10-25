function SpatializedMat = SpatializeSound(SoundMat,sr,ITD)
%   Spatialize a sound matrix using Inter-temporal Difference. SOUNDMAT is the sound matrix (2 channels only), SR the sampling rate, and ITD the inter-temporal
%   difference (in seconds).

if nargin<2
    ITD =0.0003; %seconds
end

ITDsr = floor(ITD*sr);

SoundMat = reshape(SoundMat,numel(SoundMat)/2,2);

SpatializedMat = [[SoundMat(:,1);zeros(ITDsr,1)]+[zeros(ITDsr,1);SoundMat(:,2)],[zeros(ITDsr,1);SoundMat(:,1)]+[SoundMat(:,2);zeros(ITDsr,1)]];