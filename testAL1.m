%% Code of this paper:
% T. Tang, C. Yang, S. Yan, L. Xu and D. Chen, "A Sparse Bayesian Learning-Based Approach With Indian Buffet 
% Process Prior for Joint Wideband DOA and Frequency Band Estimation", IEEE Trans. Aerosp. Electron. Syst. .
% 
% Copyright 2025 by Tao Tang
% The code is sorted out at the University of Chinese Academy of Sciences
% Date: 01-Sept-2025
% Matlab R2022b
%% demo
clc
clear
close all;

%% signal parameter (two signals)

sig_time = 10;                      % signal duration (s)
fs = 4000;                          % sampling frequency (Hz)
sig_len = fs*sig_time;
theta = [80.73,101.28];             % signal DOAs (degree)
K = length(theta);                  % signal number

f_min1 = 400;  f_max1 = 1000;       % the lower and upper bound of frequency band of signal 1 (Hz)
f_min2 = 1400;  f_max2 = 2000;      % the lower and upper bound of frequency band of signal 2 (Hz)

% option:-5dB and 0dB SNR
% load('receivedArrayData-5dB.mat')    % read array received data -5dB
load('receivedArrayData0dB.mat')       % read array received data 0dB
%% array parameter (ULA) 
M = 16;                      % array sensor number
p = [0:M-1]';    
c = 1500;                    % acoustic speed underwater
f = 2000;        
lamda = c/f;
d = lamda/2;                 % array spacing
J = 100;                     % snapshots number
fft_num = floor(sig_len/J);  % frequency Number 

delta_f = fs/fft_num;        % frequency resolution

% monitored frequency band [400-2000]
fft_num_low = 400/delta_f+1;       
fft_num_up = 2000/delta_f+1;       
fft_num_use = fft_num_up - fft_num_low;

%% array Manifold matrix

resolution = 2;                % angle grid resolution
grid = (0:resolution:180)';

N = length(grid);

% uniform linear array (ULA)
A = zeros(M,N,fft_num);
f_temp = [0:fft_num-1]*fs/fft_num;
for ff = 1:fft_num
    for n = 1:N
        temp = exp(1i*2*pi*d*p*f_temp(ff)/c*cos(grid(n)/180*pi)); 
        A(:,n,ff) = temp; 
    end
end

%% test algorithm

AUsed = A(:,:,fft_num_low+1:fft_num_up);

params.Y = YUsed;       % [M*T*fft_num_use]: array received data in frequency domain
params.A = AUsed;
params.resolution = resolution/180*pi;
params.rho = 1e-2;
% signal variance initialization
temp1 = 0;
for ff = 1:fft_num_use
    temp = AUsed(:,:,ff)'*YUsed(:,:,ff);
    temp1 = mean(abs(temp),2)/50;
    alpha_mean(:,ff) = temp1;
end
params.alpha = alpha_mean;
%
params.maxiter = 2000;                % maximum iteration
params.tolerance = 3e-5;
% noise variance initialization 
sigma2_mean = 0;
for ff = 1:fft_num_use
    temp = mean(var(YUsed(:,:,ff)))/100;
    sigma2_mean = sigma2_mean + temp/fft_num;
end
params.sigma2 = sigma2_mean;     
%

res = SBLIBP(params);

%% feature-to-signal 
alphaEachF = (res.alpha);
Zema = res.Zema;
phi = res.phi;

alpha_all = sum(phi,1);
[peak,locs] = findpeaks(alpha_all);
[~,idx] = sort(peak,'descend');
gridIdx_signals = locs(idx(1:K));

signalIdx = []; 
for kk = 1:K
    idxgrid = gridIdx_signals(kk);
    idx1 = 0; maxpeak = 0;
    for ii = 1:size(res.phi,1)
        phi = res.phi;
        maxValue = phi(ii,idxgrid);
        if maxValue>maxpeak
            idx1 = ii;
            maxpeak = maxValue;
        end
    end  
    signalIdx = [signalIdx, idx1];
end
% indicate frequency bands for each signal
bandSignal = res.Zema(:,signalIdx); 
%% Plot

figure
[X,Y] = meshgrid(f_temp(fft_num_low+1:fft_num_up), grid);
mesh(X,Y,alphaEachF.','EdgeColor','none','FaceColor','interp');view([0,0,1])
xlabel('frequency(Hz)')
ylabel('DOA(degree)')

figure
plot(grid, phi(signalIdx,:)/max(max(phi(signalIdx,:))))
legend(['feature ',num2str(signalIdx(1))],['feature ',num2str(signalIdx(2))])
xlabel('DOA(degree)')
ylabel('amplitude')

figure
temp = res.Zema;
temp(:,size(temp,2)+1) = 0;
idxtemp = setdiff(1:res.featureNum,signalIdx);
temp(:,idxtemp) = 0;
[X,Y] = meshgrid(f_temp(fft_num_low+1:fft_num_up), 0:res.featureNum);
mesh(X, Y, temp.','FaceColor','flat');
view([0,0,1])
xlabel('frequency(Hz)')
ylabel('Index of feature')
%     set(gca,'YTick',1:1:res.featureNum)

aaa = sort(unique([0:5:res.featureNum,signalIdx(1),signalIdx(2)]));
aaa(find(aaa == 0)) = [];
set(gca,'YTick',aaa);
%     ytricks([1:5:res.featureNum,signalIdx(1),signalIdx(2)])
%% DOA refinement

doa_coarse = grid(gridIdx_signals);        % on-grid DOA estimate

resolution_fine = 0.001;
for ii = 1:length(doa_coarse)
    grid_fine(ii,:) = doa_coarse(ii)-resolution : resolution_fine: doa_coarse(ii)+ resolution;
end
N1 = size(grid_fine,2);
% finer grid for DOA refinement
A_fine = zeros(M,N1,fft_num,length(doa_coarse));
for ii = 1:length(doa_coarse)
    for ff = 1:fft_num
        temp = exp(1i*2*pi*d*p*f_temp(ff)/c*cos(grid_fine(ii,:)/180*pi)); % array manifold matrix for DOA refinement
        A_fine(:,:,ff,ii) = temp;
        B_fine(:,:,ff,ii) = -1i * 2 * pi * d * p * f_temp(ff)/c * sin(grid_fine(ii,:)/180*pi) .* temp;
    end
end
A_fine = A_fine(:,:,fft_num_low+1:fft_num_up,:);
B_fine = B_fine(:,:,fft_num_low+1:fft_num_up,:);

doa_spectral = zeros(ii,N1);
A_coarse = AUsed(:,gridIdx_signals,:);
alpha_coarse = alpha_all(gridIdx_signals);

for ii = 1:length(doa_coarse)
    temp = zeros(1,N1);
    fft_numIdx = find(bandSignal(:,ii)>0).';
    for ff = fft_numIdx
        if real(alpha_all(gridIdx_signals(ii)-1)) < real(alpha_all(gridIdx_signals(ii)+1))
            tempalpha = alphaEachF(ff, gridIdx_signals(ii))*AUsed(:,gridIdx_signals(ii),ff)*AUsed(:,gridIdx_signals(ii),ff)' + alphaEachF(ff, gridIdx_signals(ii)+1)*AUsed(:,gridIdx_signals(ii)+1,ff)*AUsed(:,gridIdx_signals(ii)+1,ff)';
        else
            tempalpha = alphaEachF(ff, gridIdx_signals(ii))*AUsed(:,gridIdx_signals(ii),ff)*AUsed(:,gridIdx_signals(ii),ff)' + alphaEachF(ff, gridIdx_signals(ii)-1)*AUsed(:,gridIdx_signals(ii)-1,ff)*AUsed(:,gridIdx_signals(ii)-1,ff)';
        end
        
        Rf = (YUsed(:,:,ff)*YUsed(:,:,ff)')/J;

%             C = res.sigma2 * eye(M) + AUsed(:,:,ff) * diag(alpha) *AUsed(:,:,ff)' - alpha(locs(ii))*AUsed(:,locs(ii),ff)*AUsed(:,locs(ii),ff)';
        C = res.sigma2 * eye(M) + AUsed(:,:,ff) * diag(alphaEachF(ff,:)) *AUsed(:,:,ff)' - tempalpha;
        C_inv = inv(C);
        for n = 1:N1
            afine = A_fine(:,n,ff,ii);
            bfine = B_fine(:,n,ff,ii);
            temp(n) = temp(n) + (afine' * C_inv * ((afine*afine')*C_inv*Rf - Rf*C_inv*(afine*afine')) * C_inv * bfine) ./ (afine'*C_inv*Rf*C_inv*afine); 
        end
     end  
     doa_spectral(ii,:) = 1./abs(real(temp));
end
[~,idx] = max(doa_spectral,[],2);
for ii = 1:length(doa_coarse)
    doa(ii) = grid_fine(ii,idx(ii)); 
end
doa = sort(doa);

% DOA estimation results 
doa_estimation = doa.'
