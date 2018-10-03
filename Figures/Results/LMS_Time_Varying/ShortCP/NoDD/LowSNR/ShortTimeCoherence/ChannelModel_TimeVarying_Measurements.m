% clear all;
close all;
clc;
addpath('functions');
%% Defining some useful variables
%Constants
C = physconst('light'); %ms

% Modem variables
M = 4; % M-QAM;
FFTLength = 2^7;
BitsPerSymbol = log2(M);
SymbolsPerFrame = 1000;
% SymbolsPerFrame = round(1500*8/BitsPerSymbol/FFTLength); % Ethernet v2 MTU?
CPLength = 10;
frameSize = FFTLength*SymbolsPerFrame*BitsPerSymbol;
NoOfFrames = 1;

% Tranmission and Reception Variables
f0 = 2.45e9; % Carrier frequency [Hz]
v = 10*0.013; % UE velocity [m/s] (speed of a garden snail?)
% v = 0;
B = 0.5e6; % OFDM Symbol Bandwidth [Hz]
T = 1/B; % Sample period [s]
Ts = T * (FFTLength + CPLength); % OFDM Symbol Period [s]

% Equalisation variables
mu = 0.1;
% mu = 0.5;
decision_directed = false;
training_symbols = SymbolsPerFrame;
% training_symbols = 1;
NLMS_Numerical_Conditioning_Constant = 0.005;

% Wireless Channel Variables
EbNoVec = [12,18,20,22,24,26,28];

% Fading parameters
A = -20; % difference between maximum and negligible path power. dB
A_linear = 10^(A/10);
tau_d = 0.75*T; % RMS delay spread
T_m = -tau_d*log(A_linear); % Maximum delay spread. s
f_0 = 1/T_m; % coherence bandwidth. Hz

fd = v/(C/f0); % Doppler frequency [Hz]
T_0 = 9/(16.*pi.*fd); % 0.5 coherence time.[s]

pathDelays = [0,1,2].*T;
p = (1./tau_d).*exp(-1.*pathDelays./tau_d);
g = sqrt(T.^2.*p);
pathGains = 10.*log10(g);

%Simulation Measurement variables
NoOfRealisations = 100;
MeanSquareError = zeros(FFTLength, SymbolsPerFrame, NoOfRealisations);
berTmpVec = zeros(length(EbNoVec), NoOfRealisations, NoOfFrames);

% Wiener Solution Variables
WienerRx = zeros(FFTLength, SymbolsPerFrame);
WienerFilterCoeff = zeros(FFTLength,SymbolsPerFrame);
WienerError = zeros(FFTLength, SymbolsPerFrame, NoOfRealisations);
WienerMSE = zeros(FFTLength);
WienerBERTmpVec = zeros(length(EbNoVec), NoOfRealisations, NoOfFrames);
%% Defining System Objects
% QPSK modulation
QPSKmod = comm.QPSKModulator('BitInput', true);
qpskDemod = comm.QPSKDemodulator('BitOutput',true);

% OFDM Modulation
ofdmQpskMod = comm.OFDMModulator( ...
    'FFTLength',            FFTLength, ...
    'NumGuardBandCarriers', [0;0], ...
    'InsertDCNull',         false, ...
    'PilotInputPort',       false, ...
    'CyclicPrefixLength',   CPLength, ...
    'NumSymbols',           SymbolsPerFrame, ...
    'NumTransmitAntennas',  1);
ofdm4QAMDemod = comm.OFDMDemodulator(ofdmQpskMod);

rayChan = comm.RayleighChannel( ...
        'PathDelays', pathDelays, ...
        'AveragePathGains', pathGains, ...
        'NormalizePathGains', true, ...
        'PathGainsOutputPort', true, ...
        'MaximumDopplerShift', fd, ...
        'SampleRate', B, ...
        'DopplerSpectrum', doppler('Jakes'));
    
% rayChan.Visualization = 'Frequency response';
% AWGN Channel
awgnChannel = comm.AWGNChannel( ...
    'NoiseMethod', 'Variance', ...
    'VarianceSource', 'Input port');
    
ricChan=comm.RicianChannel( ...
    'PathDelays', pathDelays, ...
    'AveragePathGains', pathGains, ...
    'NormalizePathGains', true, ...
    'PathGainsOutputPort', true, ...
    'MaximumDopplerShift', fd, ...
    'KFactor', 3, ...
    'DirectPathDopplerShift', 0, ...
    'DirectPathInitialPhase', 0, ...
    'SampleRate', B, ...
    'DopplerSpectrum', doppler('Jakes'));   
% ricChan.Visualization = 'Frequency response';

% multipathChan = rayChan;
multipathChan = ricChan;

% gainScope = dsp.TimeScope( ...
%     'SampleRate', multipathChan.SampleRate, ...
%     'TimeSpan', SymbolsPerFrame/multipathChan.SampleRate, ...
%     'Name', 'Multipath Gain', ...
%     'ShowGrid', true, ... 
%     'YLimits', [-40, 10], ...
%     'YLabel', 'Gain (dB)');

%% Simulation
% Simulation
% MSE_flag for toggling between mean square error and bit error rate
% plots
MSE_flag = 1;
if (MSE_flag == 1)
    M = 1;
else
    M = length(EbNoVec);
end


for m = 1:M
    for k = 1:NoOfRealisations
        multipathChan.reset();
        for i = 1:NoOfFrames
            decision_directed = false;
            x1=newRandomBinaryFrame(frameSize);
            Tx_QPSK = reshape(QPSKmod(x1), [FFTLength,SymbolsPerFrame]);
            Tx = ofdmQpskMod(Tx_QPSK);

            % AWGN channel system variables       
            snr = EbNoVec(m) + 10*log10(BitsPerSymbol);
            powerDB = 10*log10(var(Tx));
            noiseVar = 10.^(0.1*(powerDB-snr));

            [TxMultipath, multipathTaps] = multipathChan(Tx);
            Noiseless_Rx = ofdm4QAMDemod(TxMultipath);
            TxMultipath = awgnChannel(TxMultipath, noiseVar);
            Rx = ofdm4QAMDemod(TxMultipath);
            
            T0TsRatioOrderOfMagnitude = floor(log10(T_0/Ts));
            FilterLength = 1*10^(T0TsRatioOrderOfMagnitude-1);
            for l = 1:(SymbolsPerFrame/FilterLength)
                for j = 1:size(Rx,1)
                    [WienerRx(j,(l-1)*FilterLength+1:l*FilterLength), ...
                        WienerFilterCoeff(j,...
                        (l-1)*FilterLength+1:l*FilterLength)] = ...
                        One_Tap_Wiener_Filter(Rx(j,...
                        (l-1)*FilterLength+1:l*FilterLength), ...
                        Tx_QPSK(j,...
                        (l-1)*FilterLength+1:l*FilterLength), ...
                        Noiseless_Rx(j,...
                        (l-1)*FilterLength+1:l*FilterLength), noiseVar);
                    WienerError(j,...
                        (l-1)*FilterLength+1:l*FilterLength,k) = ...
                        WienerRx(j,...
                        (l-1)*FilterLength+1:l*FilterLength)-...
                        Tx_QPSK(j,(l-1)*FilterLength+1:l*FilterLength);
    %                     WienerMSE(j) = var(WienerError(j,:));
                end
            end
            
            % Equaliser
            e = zeros(FFTLength, size(Rx, 2)); % Initialising e (not necessary but a small optimisation)
            H_hat = zeros(FFTLength,FilterLength);
%             for j = 1:size(Rx,1)
%                 [~, H_hat(j,:)] = ...
%                     One_Tap_Wiener_Filter(Rx(j,1:FilterLength), ...
%                     Tx_QPSK(j,1:FilterLength), ...
%                     Noiseless_Rx(j,1:FilterLength), noiseVar);
%             end
            H_hat = mean(H_hat,2);
            for j = 1:size(Rx,2)
                if(~decision_directed)
                    % One tap LMS
                    e(:,j) = Tx_QPSK(:,j) - conj(H_hat).*Rx(:,j);
%                     H_hat = H_hat + (mu./(NLMS_Numerical_Conditioning_Constant+abs(Rx(:,j)).^2)).*Rx(:,j).*conj(e(:,j));
                    H_hat = H_hat + mu.*Rx(:,j).*conj(e(:,j));
                    Rx(:,j) = conj(H_hat).*Rx(:,j); 
                    if(j == training_symbols)
                        decision_directed = ~decision_directed;
                    end
                else
                    equalisedRx = conj(H_hat).*Rx(:,j);
        %             scatterplot(Rx(:,j));
                    dd = qpskDemod(equalisedRx);
                    dd = QPSKmod(dd);
        %                 figure();
        %                 scatterplot(dd);
                    e(:, j) = dd - equalisedRx;
                    H_hat = H_hat + mu.*Rx(:,j).*conj(e(:,j));
                    % NLMS
%                     H_hat = H_hat + (mu./(NLMS_Numerical_Conditioning_Constant+abs(Rx(:,j)).^2)).*Rx(:,j).*conj(e(:,j));
                    Rx(:,j) = equalisedRx;
                end
                if (MSE_flag == 1)
                    MeanSquareError(:,j,k) = Tx_QPSK(:,j) - Rx(:,j);
                end
            end
        %         figure();
    %         scatterplot(reshape(Rx, [FFTLength*SymbolsPerFrame 1]));
            if (MSE_flag == 0)
                Rx = qpskDemod(reshape(Rx, [FFTLength*SymbolsPerFrame 1]));
                berTmpVec(m,k,i) = sum(Rx ~= x1)/frameSize;
                WienerRx1 = qpskDemod(reshape(WienerRx, [FFTLength*SymbolsPerFrame 1]));
                WienerBERTmpVec(m,k,i) = sum(WienerRx1 ~= x1)/frameSize;
            end
        end
    %     if ( k < 15 )
    %         scatterplot(reshape(Rx,[FFTLength*SymbolsPerFrame 1]));  
    %         pause;
    %     end
    end
end
if (MSE_flag == 0)
    berVec = mean(mean(berTmpVec,3),2);
    WienerBERVec = mean(mean(WienerBERTmpVec,3),2);
    semilogy(EbNoVec, berVec);
    hold on;
    semilogy(EbNoVec,WienerBERVec);
    xlabel("E_b/N_0");
    ylabel("BER");
    title("BER vs E_b/N_0 for a Rayleigh fading channel");
else
    MeanSquareError = mean((MeanSquareError.*conj(MeanSquareError)),3);
%     WienerMSE = mean(mean(WienerError.*conj(WienerError),2),3);
    WienerMSE = mean(WienerError.*conj(WienerError),3);
    for i = 1:1
        figure();
        hold on;
        plot(MeanSquareError(i,:));
        plot(1:size(MeanSquareError,2),WienerMSE(i,:),'r');
        xlabel("OFDM Symbols")
        ylabel("Mean Square Error")
        title("Mean Square Error vs No. of Symbols");
    end
    legend('LMS', 'Wiener Solution');
end

