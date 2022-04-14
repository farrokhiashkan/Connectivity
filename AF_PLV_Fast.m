%% PLV caculation
% filter Type: Morlet
% Date: '14-Jan-2022'
% Date: '14-Apr-2022'
% Ashakan Farrokhi
% IUST

%=== Input:
%      X: [Time * Trial * Channel] 
%      Guard:  sample
%      fs: sampling frequency (Hz);
%      Step: time stride (sample);
%      W_Width: Time segment width (sample);
%      freqs2use: frequency vector (Hz)
%      num_cycles: vector used for morlet (length(num_cycles) = length(freqs2use))


%=== OutPut:
%      PLV: (complex) [channel * channel * Time(Step) * N_trial * length(freqs2use)]
%      TimeIndex: [Time(Step) 2] ; (ms)


%=== Example:
% clear 
% clc
% X = randn(700,10,63);   % X: [Time * Trial * Channel] 
% Guard = 50;
% fs=1000;
% 
% %=== Frequency resolution:
% freqs2use  = logspace(log10(1),log10(200),60)';
% freqs2use    = freqs2use(dsearchn(freqs2use,1):dsearchn(freqs2use,40)); % Hz
% num_cycles    = logspace(log10(3),log10(12),length(freqs2use));
% 
% Step = 1;
% W_Width =size(X,1)-(2*Guard);
% tic
% [PLV,TimeIndex] = AF_PLV_Fast2(X,freqs2use,num_cycles,fs,Guard,W_Width,Step);
% toc
%=========================================


function [PLV,TimeIndex] = AF_PLV_Fast(X,freqs2use,num_cycles,fs,Guard,W_Width,Step)

%=== calculate time index for accross time calculation
% Step = 20;
% W_Width =250;

TimeIndex = size(X,1)-(2*Guard):-1:1;
TimeIndex = buffer(TimeIndex,W_Width,W_Width-Step);
indFirst = find(TimeIndex(1,:)==0,1,'last');
TimeIndex(:,indFirst+1:end);
TimeIndex = flip(TimeIndex,2);
TimeIndex = flip(TimeIndex,1);
indLast = find(TimeIndex(end,:)==0,1,'first');
TimeIndex(:,indLast:end)=[];
clear indFirst indLast 
%=========================================



%%


N_Trial = size(X,2);
N_Channel = size(X,3);
L_Trial = size(X,1);


%=== Cat data: across "Channel"
X_Cat = reshape(X, [L_Trial*N_Trial*N_Channel,1]);
clear X

%=== Frequency resolution:
% freqs2use  = logspace(log10(1),log10(200),60)';
% freqs2use    = freqs2use(dsearchn(freqs2use,1):dsearchn(freqs2use,40)); % Hz

%=== Design wavelet filters:
% num_cycles    = logspace(log10(3),log10(12),length(freqs2use));
wavelet_fft=cell(length(freqs2use),1);
for j=1:length(freqs2use) 

    f = freqs2use(j);
    nc = num_cycles(j);
    fmin = freqs2use(1);
    nc1 = num_cycles(1);
    
    ts=fs*(nc1/(2*pi*fmin));
    time          = -round(4*ts)/fs:1/fs:round(4*ts)/fs;
    half_wavelet  = (length(time)-1)/2;

    n_wavelet     = length(time);
    n_data        = size(X_Cat,1);
    n_convolution = n_wavelet+n_data-1;

    % create wavelet and take FFT
    s = nc/(2*pi*f);
    wavelet_fft{j,1} = fft( exp(2*1i*pi*f.*time) .* exp(-time.^2./(2*(s^2))) ,n_convolution);

end %j
clear f nc fmin nc1 ts time n_wavelet n_data  s j
wavelet_fft = cat(1,wavelet_fft{:});
wavelet_fft = transpose(wavelet_fft);
%=========================================


%=== Data FFTs
data_fft1 = fft(X_Cat,n_convolution);
clear X_Cat

%=== Filtering Data using Wavelet filter bank
convolution_result_fft = ifft(wavelet_fft.* repmat(data_fft1,1,length(freqs2use)) ,n_convolution,1);
convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet,:);
phase = angle(convolution_result_fft);
phase = exp(1i.*phase);
clear wavelet_fft convolution_result_fft data_fft1

phase = reshape(phase,L_Trial,N_Trial,N_Channel,length(freqs2use));
phase = phase(Guard+1:end-Guard,:,:,:);


%=== PLV:
PLV = nan(N_Channel,N_Channel,size(TimeIndex,2),N_Trial,length(freqs2use));
for i=1:N_Channel
    
    for j=i+1:N_Channel
        DifPhase = phase(:,:,i,:) ./ phase(:,:,j,:);
        
        % accross time:
        DifPhase = squeeze(DifPhase(TimeIndex(:),:,:,:));
        DifPhase = reshape(DifPhase,size(TimeIndex,1),size(TimeIndex,2),size(DifPhase,2),size(DifPhase,3));
        PLV(i,j,1:size(TimeIndex,2),:,:) =  squeeze( abs(mean(DifPhase,1)));
    end  %j
    
end %i
TimeIndex = TimeIndex([1 end],:)';
clear DifPhase i j phase


end



























%%