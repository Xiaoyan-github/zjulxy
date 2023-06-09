%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation-free spatial-domain image reconstruction of structured illumination microscopy   
% Copyright (C) 2023 Xiaoyan Li           
%                                                   
% Please cite:                                                                                         
% Xiaoyan Li et al., "Estimation-free spatial-domain image reconstruction of structured illumination microscopy" xx.xx.xx                                                                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------- clear workspace ---------- 
clc;clear;close all;
%% ---------- data information ----------
addpath('.\functions\');
a_num = 3; % number of pattern orientations
p_num = 3; % phase shift times for each pattern orientation
pixelSize = 65; % pixel size. unit: nm
lambdaEmi = 525; % fluorescence emission wavelength (emission maximum). unit: nm (670/525)
NA = 1.49; % numerical aperture of the objective
%% -----------read image file------------
filepath='C:\Users\mk\Desktop\EFSD-SIM for github\testdata\';
filename='1_G_';
fileformat='tif';
for ii=1:a_num * p_num
    Rawimage(:,:,ii)=single(imread([filepath,filename,num2str(ii),'.',fileformat]));
end
[xsize,ysize]=size(Rawimage(:,:,1,1));
% % % interpolation
scaleFactor = 2;
pixelSize = pixelSize/scaleFactor;
xsize = xsize * scaleFactor;
ysize = ysize * scaleFactor;
image = gpuArray((zeros(size(Rawimage))));
for ii=1:a_num * p_num
    image(:,:,ii) = Rawimage(:,:,ii);    
end
sample = imresize(image,scaleFactor,'Method','bilinear');
widefield = mean(sample,3);
widefield = widefield./max(widefield(:));
% % % get notch filter
[Y,X]=meshgrid(1:ysize,1:xsize);
rad = hypot(X-floor(xsize/2+1),Y-floor(ysize/2+1));
cyclesPerMicron = 1/(xsize*pixelSize/1000);
cycl = rad.*cyclesPerMicron;
amp = 1; % ranges from 0 to 1, depends on the background level of raw data, e.g. 1
sgima = 1; % positive number, e.g. 1
attFun = gpuArray(single(1-amp*exp(-cycl.^2/(2*sgima^2))));
[~, OTFde] = generatePSFOTFGPU(xsize, ysize, pixelSize, NA, lambdaEmi);
notchFilter = attFun.*conj(OTFde);
notchFilter = notchFilter./max(notchFilter(:));
% % % get Wiener filter
wnrFactorWF = 0.5;    % wiener factor, e.g. 0.5
wnrFilterWF = conj(OTFde)./(abs(OTFde).^2 + wnrFactorWF);
wnrFactorImg = 0.5;     % wiener factor, e.g. 0.5
wnrFilterImg = conj(OTFde)./(abs(OTFde).^2 + wnrFactorImg);
%% ---------- main reconstruction procedure of EFSD-SIM ---------- 
tic
% processing on raw measurements to get calculated patterns and filtered measurements
for ii=1:a_num * p_num
    Filteredsample(:,:,ii) = real(ifft2(ifftshift((fftshift(fft2(sample(:,:,ii)))).*wnrFilterWF)));
    Filteredwidefield = real(ifft2(ifftshift((fftshift(fft2(widefield))).*wnrFilterImg)));
    pattern(:,:,ii) = ( Filteredsample(:,:,ii))./( Filteredwidefield + 0.1);
    sample(:,:,ii) = real(ifft2(ifftshift((fftshift(fft2(sample(:,:,ii)))).*notchFilter)));    
end
% reconstruction
patternSubMean = pattern - mean(pattern,3);
sampleSubMean = sample - mean(sample,3);
EFSDRec = mean(sampleSubMean.*patternSubMean,3);
toc
%% --------------------------------------------------------------
% % % save EFSD-SIM image
EFSDSIM = EFSDRec.*(EFSDRec>0);
EFSDSIM = EFSDSIM./max(EFSDSIM(:));
EFSDSIM = mat2gray(EFSDSIM);
imwrite(EFSDSIM,strcat(filepath,'EFSD-SIM.tif'));
% % % save widefield image
widefield = mat2gray(widefield);
imwrite(widefield,strcat(filepath,'widefield.tif'));