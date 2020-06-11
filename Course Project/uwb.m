%%Subsampling architecture for UWB. (c) 2020, Vivek Jain and Swapnil Sayan Saha
% IEEE CM-1 code has been referenced from: 
% Molisch, Andreas F., et al. "IEEE 802.15. 4a Channel Model-Final Report." IEEE P802 15.04 (2004): 0662.
%% Part 1: Analysis of delay offset on SNR at the output of MF

clear;
clc;
num_trials = 1000; %number of trials
B=10^9; %signal bandwidth (1 GHz)
Ts=1/B; %symbol duration (1 ns)
fc=3.5*10^9; %center frequency (3.5 GHz)
ts= 0.001; %timestep (nS)
Fs=1/(ts*10^(-9)); %sampling frequency
Nf = 10; %number of frames per symbol (10, given in paper)
Tf=Ts/Nf; %frame time Ts = Nf*Tf
t = 0:ts*10^(-9):Ts; %time range
SNR = 20; %ideal SNR
t0=0:10^-14:(ts*10^-9)/2; %offsets
snrs=zeros(3,length(t0)); %store SNR for varying timing offset

for z = 1:num_trials
    %generate pseudorandom sequence
    d = [ones(1,Nf^2/2), -ones(1,(Nf^2/2)+1)];
    d = d(randperm(length(d)));
    
    %generate information signal (BPSK)
    b_symbol = zeros(1, Nf+1);
    for m = 1:Nf+1
        b_set = [-1 1];
        info_choice = randperm(2,1);
        b_symbol(m) = b_set(info_choice);
    end
    
    %generate unit energy transmitted signal (subsampled and at baseband)
    s_tx = zeros(length(t),1);
    for i = 1:length(t)
        for j = -Nf^2/2:1:Nf^2/2
            s_tx(i) = s_tx(i) + sqrt(1/(2*Nf)).*d(j+(Nf^2/2)+1).*b_symbol((Nf/2)+1+floor(j/Nf)).*gauspuls(t(i)-j*Tf,fc,(B/fc));
        end
    end
    
    %generate analytic signal and matched filter response
    s_t = hilbert(s_tx); %analytic signal without timing offset
    h = conj(real(fliplr(s_t)))+(1i*conj(imag(fliplr(s_t)))); %impulse response of matched filter
    f=Fs*(-(length(t)-1)/2:(length(t)-1)/2);
    %induce timing offset and calculate matched filter output
    for q = 1:length(t0) %timing offset
        %add timing offset
        y = fftshift(fft(s_t,length(s_t)));
        y_PS = y.*(exp(-1i*2*pi*f*t0(q)).');
        y_PS(1) = real(y_PS(1));
        y = (ifft(ifftshift(y_PS),length(s_t)));
        
        %generate matched filter output
        mr= sum(real(y).*real(fliplr(h)))+sum(imag(y).*imag(fliplr(h)));
        mi= -sum(real(y).*imag(fliplr(h)))+sum(imag(y).*real(fliplr(h)));
        
        %generate noise variance
        n_var = (sum(real(s_t).^2).*2)./(10^(SNR/10)); %noise variance
        %noise = (randn(1,1)+1i.*randn(1,1)).*sqrt(n_var).*sqrt(0.5); %AWGN
        
        %calculate SNR
        snrs(1,q)=snrs(1,q)+(mr./(n_var./2))./num_trials; %real part
        snrs(2,q)=snrs(2,q)+(mi./(n_var./2))./num_trials; %imaginary part
        snrs(3,q)=snrs(3,q)+(sqrt(mr.^2+mi.^2)./(n_var))./num_trials; %magnitude
    end
end
%Plot results
figure
plot(t0./t0(end), 10.*log10(abs(snrs(1,:))),'LineWidth',2);
hold on;
plot(t0./t0(end), 10.*log10(abs(snrs(2,:))),'LineWidth',2);
plot(t0./t0(end), 10.*log10(abs(snrs(3,:))),'LineWidth',2,'Marker','o');
hold off;
grid minor;
legend('real','imag','mag');
ylim([0 30]);
xlabel('normalized to Ts');
title('Analytic matched filter output for SNR = 20 db');
ylabel('SNR (db)');
saveas(gcf,'Part_1','epsc');

%% Part 2: Analysis of delay offset on BER at the output of MF

clearvars -except B Ts fc Nf Tf
num_trials = 5000; %number of trials
ts= 0.01; %timestep (nS)
Fs=1/(2*ts*10^(-9)); %sampling frequency
t = 0:ts*10^(-9):Ts; %time range
SNR = 0:2:20; %SNR range (0 - 20 db)
bers=zeros(length(SNR),3); %store BER
t0=[(ts/2)*10^-9, 0, (ts)*10^-9]; %offsets

%generate information signal (BPSK)
b_symbol = zeros(1, num_trials);
for m = 1:num_trials
    b_set = [-1 1];
    info_choice = randperm(2,1);
    b_symbol(m) = b_set(info_choice);
end

for k=1:1:length(SNR)
    n_var = (1.*2)./(10^(SNR(k)/10)); %variance of AWGN
    noise = (randn(num_trials,length(t))+1i.*randn(num_trials,length(t))).*sqrt(n_var).*sqrt(0.5); %generate AWGN
    for z = 1:num_trials
        s_tx = b_symbol(z).*(gauspuls(t,fc,(B/fc)).'); %generate transmitted signal
        s_tx=s_tx/rms(s_tx); %normalize the signal to unit energy
        h= ieee_cm_1(1,fc/10^9,B/10^9,fc/10^9); %generate channel impulse response (1 channel, antenna resonant frequency same as fc)
        [los,lind]=max(h(1:length(t))); %find LoS component and index
        s_ref=gauspuls(t(lind),fc,(B/fc)).'; %generate reference pulse for LoS component and index
        s_rx=conv(h(1:length(t),1),s_tx); %convolve transmitted signal with channel impulse response
        f=Fs*(-(length(t)-1)/2:(length(t)-1)/2);
        
        for q = 1:length(t0)
            %Add noise and timing offset
            y = fftshift(fft(s_rx(1:length(t)),length(t)))/length(t);
            y1 = y.*(exp(-1i*2*pi*f*t0(q)).');
            y1(1) = real(y1(1)); % band aid fix
            y = length(t)*(ifft(ifftshift(y1),length(t)));
            y=y+noise(z,:).';
            
            %Detect the signal at output of matched filter
            y_hat=sign((real(y(lind).*s_ref)));
            %calculate BER
            if y_hat~=b_symbol(z)
                bers(k,q)=bers(k,q)+1/num_trials; 
            end
        end
    end
end
%plot results
figure
for i=1:1:3
    semilogy(SNR,bers(:,i),'LineWidth',2,'Marker','o');
    hold on;
end
hold off;
grid minor;
legend('Offset=0', 'Offset=0.5Ts',  'Offset=Ts');
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
title('BER vs SNR for various timing offsets')
saveas(gcf,'Part_2','epsc');











%% Function for IEEE CM-1 Channel Model


function[h] = ieee_cm_1(num_channels,fc,fs,f0)
%IEEE CM-1 Channel Model
% input:
% num_channels: number of channel impulse responses to generate
% fc: center frequency (GHz)
% fs: Bandwidth (GHz) (e.g. if fc=6, fs=8 if signal spans from 2-10 Ghz),
% f0: antenna frequency dependence (GHz) (e.g. 5)

%output: h (each column denotes for each channel in num_channels)
%clear;
randn('state',12); % initialize state of function for repeatability, seed = 12
rand('state',12); % initialize state of function for repeatability, seed = 12
%% Declare channel model parameters for IEEE CM-1 (Residential LOS)
% MPC arrival
Lam = 0.047; Lmean = 3;
lambda_1 = 1.54; lambda_2 = 0.15; beta = 0.095;
%MPC decay
Gam = 22.61; gamma_0 = 12.53; Kgamma = 0; sigma_cluster = 2.75;
% Small-scale Fading
m0 = 0.67; Km = 0; sigma_m0 = 0.28; sigma_Km = 0;
% Frequency Dependence
kappa = 1.12;
ts = 1/fs; % sampling frequency / time resolution

%% Compute continuous time impulse responses (Output: h_ct, t_ct, t0, np)
std_L = 1/sqrt(2*Lam); % std dev (nsec) of cluster arrival spacing
% std dev (nsec) of ray arrival spacing
std_lam_1 = 1/sqrt(2*lambda_1);
std_lam_2 = 1/sqrt(2*lambda_2);
h_len = 1000; % there must be a better estimate of # of paths than this
ngrow = 1000; % amount to grow data structure if more paths are needed
h = zeros(h_len,num_channels); %h is returned as a matrix with num_channels columns, each column
%holding a random realization of the channel model (an impulse response)
t = zeros(h_len,num_channels); %t is organized as h, but holds the time instances (in nsec) of the paths whose
%signed amplitudes are stored in h
t0 = zeros(1,num_channels); %arrival time of the first cluster for each realization
np = zeros(1,num_channels); % np is the number of paths for each realization.
%the k'th realization of the channel impulse response is the sequence
%of (time,value) pairs given by (t(1:np(k),k), h(1:np(k),k))
for k = 1:num_channels % loop over number of channels
    tmp_h = zeros(size(h,1),1);
    tmp_t = zeros(size(h,1),1);
    Tc = 0; % First cluster arrival occurs at time 0
    t0(k) = Tc;
    L = max(1, poissrnd(Lmean)); % number of clusters
    cluster_index = zeros(1,L);
    path_ix = 0;
    nak_m = [];
    for ncluster = 1:L
        % Determine Ray arrivals for each cluster
        Tr = 0; % first ray arrival defined to be time 0 relative to cluster
        cluster_index(ncluster) = path_ix+1; % remember the cluster location
        gamma = Kgamma*Tc + gamma_0; % delay dependent cluster decay time
        Mcluster = sigma_cluster*randn;
        Pcluster = 10*log10(exp(-1*Tc/Gam))+Mcluster; % total cluster power
        Pcluster = 10^(Pcluster*0.1);
        Tr_len = 10*gamma;
        while (Tr < Tr_len)
            t_val = (Tc+Tr); % time of arrival of this ray
            h_val = Pcluster/gamma*exp(-Tr/gamma)/(beta*lambda_1+(1-beta)*lambda_2+1);
            path_ix = path_ix + 1; % row index of this ray
            if path_ix > h_len
                % grow the output structures to handle more paths as needed
                tmp_h = [tmp_h; zeros(ngrow,1)];
                tmp_t = [tmp_t; zeros(ngrow,1)];
                h = [h; zeros(ngrow,num_channels)];
                t = [t; zeros(ngrow,num_channels)];
                h_len = h_len + ngrow;
            end
            tmp_h(path_ix) = h_val;
            tmp_t(path_ix) = t_val;
            if rand < beta
                Tr = Tr + (std_lam_1*randn)^2 + (std_lam_1*randn)^2;
            else
                Tr = Tr + (std_lam_2*randn)^2 + (std_lam_2*randn)^2;
            end
            % generate log-normal distributed nakagami m-factor
            m_mu = m0 - Km*t_val;
            m_std = sigma_m0 - sigma_Km*t_val;
            nak_m = [nak_m, lognrnd(m_mu, m_std)];
        end
        Tc = Tc + (std_L*randn)^2 + (std_L*randn)^2;
    end
    for path = 1:path_ix
        h_val = (gamrnd(nak_m(path), tmp_h(path)/nak_m(path))).^(1/2);
        tmp_h(path) = h_val;
    end
    np(k) = path_ix; % number of rays (or paths) for this realization
    [sort_tmp_t,sort_ix] = sort(tmp_t(1:np(k))); % sort in ascending time order
    t(1:np(k),k) = sort_tmp_t;
    h(1:np(k),k) = tmp_h(sort_ix(1:np(k)));
end
h_ct = h;
t_ct = t;
%clearvars -except htemp beta f0 fc fs Gam gamma_0 kappa Kgamma km Lam lambda_1 lambda_2 lambda_mode Lmean m0 nlos num_channels sigma_cluster sigma_Km sigma_m0 ts h_ct t_ct np t0

%% Change to Complex Baseband Channel
h_ct_len = size(h_ct, 1);
phi = zeros(h_ct_len, 1);
for k = 1:num_channels
    phi = rand(h_ct_len, 1).*(2*pi);
    h_ct(:,k) = h_ct(:,k) .* exp(phi .* 1i);
end

%% Reduce continuous-time result to a discrete-time result
min_Nfs = 100; % GHz
N = max( 1, ceil(min_Nfs*ts) ); % N*fs = N/ts is the intermediate sampling frequency before decimation
N = 2^nextpow2(N); % make N a power of 2 to facilitate efficient multi-stage decimation
Nfs = N / ts;
t_max = max(t_ct(:)); % maximum time value across all channels
h_len = 1 + floor(t_max * Nfs); % number of time samples at resolution ts / N
hN = zeros(h_len,num_channels); %discrete time representation of channel (time resolution ts / N)
for k = 1:num_channels
    np_k = np(k); % number of paths in this channel
    t_Nfs = 1 + floor(t_ct(1:np_k,k) * Nfs); % vector of quantized time indices for this channel
    for n = 1:np_k
        hN(t_Nfs(n),k) = hN(t_Nfs(n),k) + h_ct(n,k);
    end
end
if N > 1
    h = resample(hN, 1, N); % decimate the columns of hN by factor N
else
    h = hN;
end

%% Add Frequency Dependency

h_len = length(h(:,1));
f = [fc-fs/2 : fs/h_len/2 : fc+fs/2]./f0;
f = f.^(-2*(kappa));
f = [f(h_len : 2*h_len), f(1 : h_len-1)]';
for c = 1:num_channels
    % add the frequency dependency
    h2 = zeros(2*h_len, 1);
    h2(1 : h_len) = h(:,c); % zero padding
    fh2 = fft(h2);
    fh2 = fh2 .* f;
    h2 = ifft(fh2);
    h(:,c) = h2(1:h_len);
    % Normalize the channel energy to 1
    h(:,c) = h(:,c)/sqrt(h(:,c)' * h(:,c) );
end
end