function [usvFile, F,T,P] = USVdetect(fileName,pathToData,varargin)
%USVdetect ; adapted from USVseg (Tachibana et al) to run as a function
%without opening GUI to be nested within a script

%create input parser
%p=inputParser;

%define default conditions
%defaultSaveCond=true;

%parse optional inputs
%addOptional(p,'SaveCond',defaultSaveCond,@islogical);

%parse
%parse(p,fileName,pathToData,varargin{:});

%read parse results
%savecond=p.Results.SaveCond;

%Set default params
fftsize = 512;
prm.timestep = 0.0005; prm.freqmin = 40000; prm.freqmax = 120000;
prm.threshval = 5.0;   prm.durmin = 0.01;  prm.durmax = 0.3;
prm.gapmin = 0.020;    prm.margin = 0.015;  prm.wavfileoutput = 1;
prm.imageoutput = 1;   prm.imagetype = 0;   prm.traceoutput = 0;
prm.readsize = 300;     prm.mapL = 1000;     prm.mapH = 6000;
vocStarts = prm.readsize;
vocEnds = prm.durmax;
%fetch loading params
timestep = prm.timestep;
readsize = prm.readsize;
% read the sound file
pth = pathToData;
fn = [fileName '.wav'];
fp = [pth fn];
ai = audioinfo(fp);
wavsize = ai.TotalSamples;
fs = ai.SampleRate;
nreadsize = round(readsize*fs);
nread = ceil(wavsize/nreadsize);
%setup files to save
sfn = [fileName '_dat.csv'];
savefp = [pth sfn];
fid = fopen(savefp,'wt');
fprintf(fid,'num,start,end,duration\n');
fclose(fid);
% define values outside the loop
prevn = 0;
prevlast = 0;
med = [];
thresh = [];
threshval = prm.threshval;
freqmin = prm.freqmin;
freqmax = prm.freqmax;
durmin = prm.durmin;
durmax = prm.durmax;
gapmin = prm.gapmin;
margin = prm.margin;
tempEnd = 0;
for r=1:nread
    rng = [prevlast+1 min(r*nreadsize,wavsize)];
    if diff(rng)<fftsize*2, break; end
    [wav, fs] = audioread(fp,rng);
    %process
    %multitaper spectrogram
    step = round(timestep*fs);
    mtsp = multitaperspec(wav,fs,fftsize,timestep);
    fvec = [0; (1:(fftsize/2))'/fftsize*fs];
    tvec = ((0:(size(mtsp,2)-1))*step+fftsize/2)/fs;
    %flattening
    fltnd = flattening(mtsp);
    %save flattened spectrogram in matrix
    tempStart = tempEnd+1;
    tempEnd =tempEnd +length(fltnd);
    specMat(:,tempStart:tempEnd)=fltnd;
    %thresholding
    % threshold calculation with n*sigma (SD) of background noise
    thresh = estimatethresh(fltnd,fs,freqmin,freqmax,threshval);
    thrshd = thresholding(fltnd,fs,freqmin,freqmax,thresh);
    %detection
    %fetch detection params
    [onoffset, ~,~,contflg] = detectonoffset(thrshd,fs,timestep,gapmin,durmin,durmax,margin);
    dur = diff(onoffset,[],2);
    ronoffset = onoffset+(prevlast+1)/fs;
    nsetp = size(mtsp,2);
    %save CSV
    fid = fopen(savefp,'at');
    if ~isempty(ronoffset)
        for n = 1:size(ronoffset,1)
            fprintf(fid,'%d,%f,%f,%f\n',n+prevn,ronoffset(n,1),ronoffset(n,2),dur(n));
        end
    end
    fclose(fid);
    %calc to keep the loop going or close if done
    prevn = size(onoffset,1)+prevn;
    if contflg==1
        if isempty(onoffset)
            prevlast = rng(2) - durmax*fs;
        else prevlast = round(onoffset(end,2)*fs)+prevlast;
        end
    else
        prevlast = rng(2);
    end
end
time = ((0:(size(specMat,2)-1))*step+fftsize/2)/fs;
usvFile = savefp;
P = specMat;
T = time;
F = fvec;
save (fullfile(pathToData,strcat(fileName,'_spect.mat')),'F','T','P');
end