function [throughput_out,BER_out,PER_out,RSS_out,SNR_out] = heMUSimulateScenario(cfg,tgax,cfgSim)
%heMUSimulateScenario Simulate the HE-MU scenario
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   THROUGHPUT = heMUSimulateScenario(CFG,TGAX,CFGSIM) simulates a MU
%   scenario and returns the throughput.
%
%   CFG is the configuration of the HE-MU transmission and is a format
%   configuration object of type
%   <a href="matlab:help('wlanHEMUConfig')">wlanHEMUConfig</a>.
%
%   TGAX is a cell array containing channels models between the transmitter
%   and each user.
%
%   CFGSIM is a structure containing the simulation parameters.

%   Copyright 2018-2019 The MathWorks, Inc.

% Simulation options
c = 3e8;
f = 6.015e9;
numPackets = cfgSim.NumPackets;
noiseFloor = cfgSim.NoiseFloor;
idleTime = cfgSim.IdleTime; % us
transmitPower = cfgSim.TransmitPower; % dBm
RSS = 10*log10(((10^(transmitPower/10))*((c/(4*pi*f))^2))) - 20*log10(cfgSim.Distance);
pathloss = transmitPower - RSS;
numPL = numel(pathloss);

jammers = 4;
BER_out = zeros(numPL,jammers);
PER_out = zeros(numPL,jammers);
RSS_out = zeros(numPL,jammers);
SNR_out = zeros(numPL,jammers);
throughput_out = zeros(numPL,jammers);

for jammerType = 1:jammers
    
    standardPacket = randi([0 1], 2000, 1, 'int8'); % packet preamble
    
    % Perform synchronization with 11ac components
    timingSync = true;
    cfoCorrection = true;
    pilotPhaseTracking = true;
    
    % Get the field indices for extract fields from the PPDU
    ind = wlanFieldIndices(cfg);
    
    % Get allocation info
    allocInfo = ruInfo(cfg);
    
    s = RandStream.getGlobalStream(); % Save random stream to restore at end
    numRUs = numel(cfg.RU);
    packetErrorRate = zeros(numPL,allocInfo.NumUsers);
    throughput = zeros(numPL,allocInfo.NumUsers);
    % Simulate all path losses
    for ipl = 1:numPL
        
        % Clone the channel before running to ensure the same realization is
        % used for all path losses and schemes
        %tgaxClone = cloneChannels(tgax);
        
        % Set random substream index per iteration to ensure that each
        % iteration uses a repeatable set of random numbers
        stream = RandStream('combRecursive','Seed',0);
        stream.Substream = ipl;
        RandStream.setGlobalStream(stream);
        
        A = 10.^((transmitPower-30)/20); % Voltage gain (attenuation)
        
        % Create an instance of the AWGN channel per SNR point simulated
        fs = wlanSampleRate(cfg);
        awgnChannel = comm.AWGNChannel('NoiseMethod','Variance', ...
            'Variance',10^((noiseFloor - 30 + pathloss(ipl))/10)); % 10^((-89.9 - 30 + pathloss)/10)
        
        numPacketErrors = zeros(allocInfo.NumUsers,1);
        
        % Simulate multiple packets
        BER_PL = 0;
        RSS_PL = 0;
        SNR_PL = 0;
        
        for pktIdx = 1:numPackets
            % Generate random data for all users
            psduLength = getPSDULength(cfg) - length(standardPacket)/8;
            txPSDU = cell(allocInfo.NumUsers,1);
            for userIdx = 1:allocInfo.NumUsers
                txPSDU{userIdx} = [standardPacket; randi([0 1],psduLength(userIdx)*8,1,'int8')];
            end
            
            % Generate waveform with idle period
            tx = wlanWaveformGenerator(txPSDU,cfg,'IdleTime',idleTime*1e-6);
            
            % Scale to achieve desired transmit power
            tx = tx*A;
            mag_tx = abs(tx);
            noisevector = wgn(length(mag_tx),1,noiseFloor, 'dBm');
            %TxSNR = snr(mag_tx, noisevector);
            % disp(["SNR: ", TxSNR]);
            chanBW = cfg.ChannelBandwidth;
            
            % Per user processing
            userIdx = 1;
            
            % Loop over RUs
            BER_Packet = 0;
            RSS_Packet = 0;
            SNR_Packet = 0;
            for ruIdx = 1:numRUs
                % Loop over users within an RU
                BER_RU = 0;
                RSS_RU = 0;
                SNR_RU = 0;
                for ruUserIdx = 1:allocInfo.NumUsersPerRU(ruIdx)
                    % Pass waveform through a channel for each user
                    %rx = tgaxClone{userIdx}(tx);
                    
                    % Add noise at receiver (assuming pathloss and noise floor)
                    if jammerType == 1
                        jam = 0;
                    elseif jammerType == 2
                        jam = A*constantJammer(psduLength(userIdx)*8, cfg, idleTime);
                    elseif jammerType == 3
                        jam = A*deceptiveJammer(psduLength(userIdx)*8, cfg, idleTime, standardPacket);
                    else
                        jam = A*randomJammer(psduLength(userIdx)*8, cfg, idleTime, standardPacket, 0.5);
                    end
                    
                    if sum(jam) == 0
                        rx = awgnChannel(tx);
                    else
                        rx = awgnChannel(tx) + awgnChannel(jam);
                    end
                    
                    noise = awgnChannel(noisevector);
                    mag_rx = abs(rx);
                    RxSNR = snr(mag_rx, noise);
                    RSS = RxSNR + noiseFloor;
                    %disp(["SNR: ", RxSNR]);
                    
                    SNR_RU = SNR_RU + RxSNR;
                    RSS_RU = RSS_RU + RSS;
                    
                    if timingSync
                        % Packet detect and determine coarse packet offset
                        coarsePktOffset = wlanPacketDetect(rx,chanBW);
                        if isempty(coarsePktOffset) % If empty no L-STF detected; packet error
                            numPacketErrors(userIdx) = numPacketErrors(userIdx)+1;
                            BER_RU = 1;
                            userIdx = userIdx+1;
                            continue; % Go to next loop iteration
                        end
                    else
                        coarsePktOffset = 0;  %#ok<UNRCH>
                    end
                    
                    if cfoCorrection
                        % Extract L-STF and perform coarse frequency offset correction
                        lstf = rx(coarsePktOffset+(ind.LSTF(1):ind.LSTF(2)),:);
                        coarseFreqOff = wlanCoarseCFOEstimate(lstf,chanBW);
                        rx = helperFrequencyOffset(rx,fs,-coarseFreqOff);
                    end
                    
                    if timingSync
                        % Extract the non-HT fields and determine fine packet offset
                        nonhtfields = rx(coarsePktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
                        finePktOffset = wlanSymbolTimingEstimate(nonhtfields,chanBW);
                        
                        % Determine final packet offset
                        pktOffset = coarsePktOffset+finePktOffset;
                    else
                        pktOffset = 4; %#ok<UNRCH>
                    end
                    
                    % If packet detected outwith the range of expected delays from
                    % the channel modeling; packet error
                    if pktOffset>50
                        numPacketErrors(userIdx) = numPacketErrors(userIdx)+1;
                        userIdx = userIdx+1;
                        continue; % Go to next loop iteration
                    end
                    
                    if cfoCorrection
                        % Extract L-LTF and perform fine frequency offset correction
                        rxLLTF = rx(pktOffset+(ind.LLTF(1):ind.LLTF(2)),:);
                        fineFreqOff = wlanFineCFOEstimate(rxLLTF,chanBW);
                        rx = helperFrequencyOffset(rx,fs,-fineFreqOff);
                    end
                    
                    % HE-LTF demodulation and channel estimation
                    rxHELTF = rx(pktOffset+(ind.HELTF(1):ind.HELTF(2)),:);
                    heltfDemod = wlanHEDemodulate(rxHELTF,'HE-LTF',cfg,ruIdx);
                    [chanEst,pilotEst] = heLTFChannelEstimate(heltfDemod,cfg,ruIdx);
                    
                    % HE data demodulate
                    rxData = rx(pktOffset+(ind.HEData(1):ind.HEData(2)),:);
                    demodSym = wlanHEDemodulate(rxData,'HE-Data',cfg,ruIdx);
                    
                    % Pilot phase tracking
                    if pilotPhaseTracking
                        % Average single-stream pilot estimates over symbols (2nd dimension)
                        pilotEstTrack = mean(pilotEst,2);
                        demodSym = heCommonPhaseErrorTracking(demodSym,pilotEstTrack,cfg,ruIdx);
                    end
                    
                    % Estimate noise power in HE fields
                    ruOFDMInfo = wlanHEOFDMInfo('HE-Data',cfg,ruIdx);
                    nVarEst = heNoiseEstimate(demodSym(ruOFDMInfo.PilotIndices,:,:),pilotEst,cfg,ruIdx);
                    
                    % Equalize
                    [eqSym,csi] = heEqualizeCombine(demodSym,chanEst,nVarEst,cfg,userIdx);
                    
                    % Discard pilot subcarriers
                    rxDataUser = eqSym(ruOFDMInfo.DataIndices,:,:);
                    csiData = csi(ruOFDMInfo.DataIndices,:);
                    
                    % Demap and decode bits
                    rxPSDU = wlanHEDataBitRecover(rxDataUser,nVarEst,csiData,cfg,userIdx,'LDPCDecodingMethod','layered-bp');
                    
                    % Compare bit error
                    BER = sum(rxPSDU~=txPSDU{userIdx})/length(rxPSDU);
                    packetError = 0;
                    if BER > 0.2 % 20% BER = Packet Error
                        packetError = packetError + 1;
                        numPacketErrors(userIdx) = numPacketErrors(userIdx)+packetError;
                    end
                    userIdx = userIdx+1;
                    BER_RU = BER_RU + BER;
                end
                BER_Packet = BER_Packet + (BER_RU/allocInfo.NumUsersPerRU(ruIdx));
                RSS_Packet = RSS_Packet + (RSS_RU/allocInfo.NumUsersPerRU(ruIdx));
                SNR_Packet = SNR_Packet + (SNR_RU/allocInfo.NumUsersPerRU(ruIdx));
            end
            BER_PL = BER_PL + (BER_Packet/numRUs);
            RSS_PL = RSS_PL + (RSS_Packet/numRUs);
            SNR_PL = SNR_PL + (SNR_Packet/numRUs);
        end
        Total_BER = BER_PL/numPackets;
        Total_RSS = RSS_PL/numPackets;
        Total_SNR = SNR_PL/numPackets;
        % Calculate packet error rate (PER) at SNR point
        packetErrorRate(ipl,:) = numPacketErrors/numPackets;
        
        % Calculate the raw throughput
        txtime = (length(tx)/fs)*1e6-idleTime; % Packet duration
        for userIdx = 1:allocInfo.NumUsers
            throughput(ipl,userIdx) = 1e-6*((numPackets-numPacketErrors(userIdx))*(8*cfg.User{userIdx}.APEPLength)/(numPackets*(txtime+idleTime)*1e-6));
        end
        BER_out(ipl,jammerType) = 100*Total_BER;
        PER_out(ipl,jammerType) = 100*sum(packetErrorRate(ipl,:))/2;
        RSS_out(ipl,jammerType) = Total_RSS;
        SNR_out(ipl,jammerType) = Total_SNR;
        throughput_out(ipl,jammerType) = sum(throughput(ipl,:));
        
        disp([' JAMMER: ' num2str(jammerType)])
        
        disp([' Pathloss ' num2str(pathloss(ipl),'%2.1f') ' dB,'...
            ' -- AP throughput ' num2str(throughput_out(ipl,jammerType),'%2.1f') ' Mbps'...
            ' -- Packet Error Rate ' num2str(PER_out(ipl,jammerType),'%2.1f') '%'...
            ' -- Bit Error Rate ' num2str(BER_out(ipl,jammerType),'%2.1f') '%'...
            ' -- RSS ' num2str(RSS_out(ipl,jammerType),'%2.1f') ' dBm'...
            ' -- SNR ' num2str(SNR_out(ipl,jammerType),'%2.1f') ' dB']);
    end
    RandStream.setGlobalStream(s); % Restore random stream
end
end

% function tgaxClone = cloneChannels(tgax)
% % Clone the channels before running to ensure the same realization can be
% % used in the next simulation
% tgaxClone = cell(size(tgax));
% numUsers = numel(tgax);
% for i=1:numUsers
%     tgaxClone{i} = clone(tgax{i});
% end
% end
