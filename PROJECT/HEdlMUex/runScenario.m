function [jammerCorrectlyIdentified] = runScenario(jammerType)

clear all; close all; clc;

%% TODO: 
%   IMPLEMENT JAMMING ON SIGNAL
%       - Use Switch/Case
%       - Define adversay signal (make this a function or just pass in
%       adversary signal parameters/object to heMUSimulateScenario()
%       - Add adversary signal to recevied signal @ line 92-95 in func
%   ADD IN DECISION TREE BY ANALYZING PARAMETERS
%   ADD IN BOOLEAN YES/NO FOR IF ADVERSARY IS IDENTIFIED

%%
% Define what jamming technique is used:
%   - 0 = No Jammer
%   - 1 = Constant
%   - 2 = Deceptive
%   - 3 = Random
jammer = jammerType;

cfgOFDMA = wlanHEMUConfig(192);

numTx = 6; % Number of transmit antennas
guardInterval = 0.8; % Guard interval in Microseconds

% Configure common parameters for all users
cfgOFDMA.NumTransmitAntennas = numTx;
cfgOFDMA.GuardInterval = guardInterval;

MCS = 1; % QPSK = 1
APEPLength = 1000;

% Configure per user parameters
% STA #1 (RU #1)
cfgOFDMA.User{1}.NumSpaceTimeStreams = 2;
cfgOFDMA.User{1}.MCS = MCS;
cfgOFDMA.User{1}.APEPLength = APEPLength;

%% Channel Model Configuration
% A TGax indoor channel model is used in this example. An individual
% channel is used to simulate the link between the AP and each user. A TGax
% channel object, |tgaxBase| is created with properties relevant for all
% users. In this example, the delay profile (Model-D) and number of receive
% antennas are common for all users. Model-D is considered non-line of
% sight when the distance between transmitter and receiver is greater than
% or equal to 10 meters. This is described further in
% <docid:wlan_ref#mw_43b5900e-69e1-4636-b084-1e72dbd46293 wlanTGaxChannel>.
% A fixed seed is used for the channel to allow repeatability.

% Create channel configuration common for all users
tgaxBase = wlanTGaxChannel;
tgaxBase.DelayProfile = 'Model-D';     % Delay profile
tgaxBase.NumTransmitAntennas = numTx;  % Number of transmit antennas
tgaxBase.NumReceiveAntennas = 2;       % Each user has two receive antennas
tgaxBase.TransmitReceiveDistance = 10; % Non-line of sight distance
tgaxBase.ChannelBandwidth = cfgOFDMA.ChannelBandwidth;
tgaxBase.SampleRate = wlanSampleRate(cfgOFDMA);
% Set a fixed seed for the channel
tgaxBase.RandomStream = 'mt19937ar with seed';
tgaxBase.Seed = 5;

%%
% Next a channel is created for each user. The channel for each user is a
% clone of the |tgaxBase|, but with a unique |UserIndex| property, and is
% stored in a cell array |tgax|. The |UserIndex| property of each
% individual channel is set to provide a unique channel for each user. The
% resultant channels are used in the simulation as shown below.
%
% <<../heDownlinkMUChannelDiagram.png>>
%

% A cell array stores the channel objects, one per user
numUsers = numel(cfgOFDMA.User); % Number of users simulated in this example
tgax = cell(1,numUsers);

% Generate per-user channels
for userIdx = 1:numUsers
    tgax{userIdx} = clone(tgaxBase);
    tgax{userIdx}.UserIndex = userIdx; % Set unique user index
end

%% Beamforming Feedback
% Transmit beamforming for both OFDMA and MU-MIMO relies on knowledge of
% the channel state between transmitter and receiver at the beamformer.
% Feedback of the per-subcarrier channel state is provided by each STA by
% channel sounding. A null data packet (NDP) is transmitted by the AP, and
% each STA uses this packet to determine the channel state. The channel
% state is then fed-back to the AP. The same process is used for
% 802.11ac(TM) in the <docid:wlan_ug#example-VHTBeamformingExample
% 802.11ac Transmit Beamforming> and
% <docid:wlan_ug#example-VHTMUMIMOPrecodingExample 802.11ac
% Multi-User MIMO Precoding> examples, but an HE single user NDP packet is
% used instead of a VHT packet. In this example, the feedback is considered
% perfect; there is no noise present for channel sounding and the feedback
% is uncompressed. The |heUserBeamformingFeedback| helper function detects
% the NDP and uses channel estimation to determine the channel state
% information. Singular value decomposition (SVD) is then used to calculate
% the beamforming feedback.

% Create an NDP with the correct number of space-time streams to generate
% enough LTF symbols
cfgNDP = wlanHESUConfig('APEPLength',0,'GuardInterval',0.8); % No data in an NDP
cfgNDP.ChannelBandwidth = tgaxBase.ChannelBandwidth;
cfgNDP.NumTransmitAntennas = cfgOFDMA.NumTransmitAntennas;
cfgNDP.NumSpaceTimeStreams = cfgOFDMA.NumTransmitAntennas;

% Generate NDP packet - with an empty PSDU as no data
txNDP = wlanWaveformGenerator([],cfgNDP);

% For each user STA, pass the NDP packet through the channel and calculate
% the feedback channel state matrix by SVD.
staFeedback = cell(1,numUsers);
for userIdx = 1:numel(tgax)
    % Received waveform at user STA with 50 sample padding. No noise.
    rx = tgax{userIdx}([txNDP; zeros(50,size(txNDP,2))]);

    % Get the full-band beamforming feedback for a user
    staFeedback{userIdx} = heUserBeamformingFeedback(rx,cfgNDP);
end

%% Simulation Parameters
% Different path losses are simulated in this example. The same path loss
% and noise floor is applied to all users. For each path loss simulated, 10
% packets are passed through the channel. Packets are separated by 20
% microseconds.

cfgSim = struct;
cfgSim.NumPackets = 10;       % Number of packets to simulate for each path loss
cfgSim.Pathloss = (100);      % Path losses to simulate in dB
cfgSim.TransmitPower = 30;    % AP transmit power in dBm
cfgSim.NoiseFloor = -89.9;    % STA noise floor in dBm
cfgSim.IdleTime = 20;         % Idle time between packets in us

%% Simulation with OFDMA
% The scenario is first simulated with the OFDMA configuration and
% transmit beamforming.
%
% The steering matrix for each RU is calculated using the feedback from the
% STAs. The |heMUCalculateSteeringMatrix| helper function calculates the
% beamforming matrix for an RU given the CSI feedback.

% For each RU, calculate the steering matrix to apply
for ruIdx = 1:numel(cfgOFDMA.RU)
    % Calculate the steering matrix to apply to the RU given the feedback
    steeringMatrix = heMUCalculateSteeringMatrix(staFeedback,cfgOFDMA,cfgNDP,ruIdx);

    % Apply the steering matrix to each RU
    cfgOFDMA.RU{ruIdx}.SpatialMapping = 'Custom';
    cfgOFDMA.RU{ruIdx}.SpatialMappingMatrix = steeringMatrix;
end

%%
% The |heMUSimulateScenario| helper function performs the simulation. The
% pre-HE preamble of 802.11ax is backwards compatible with 802.11ac,
% therefore in this example the front-end synchronization components for a
% VHT waveform are used to synchronize the HE waveform at each STA. For
% each packet and path loss simulated the following processing steps occur:
%
% # A PSDU is created and encoded to create a single packet waveform.
% # The waveform is passed through an evolving TGax channel model and AWGN
% is added to the received waveform. The channel state is maintained
% between packets.
% # The packet is detected.
% # Coarse carrier frequency offset is estimated and corrected.
% # Fine timing synchronization is established.
% # Fine carrier frequency offset is estimated and corrected.
% # The HE-LTF is extracted from the synchronized received waveform. The
% HE-LTF is OFDM demodulated and channel estimation is performed.
% # The HE Data field is extracted from the synchronized received
% waveform and OFDM demodulated.
% # Common pilot phase tracking is performed to track any residual carrier
% frequency offset.
% # The phase corrected OFDM symbols are equalized with the channel
% estimate.
% # Noise estimation is performed using the demodulated data field pilots
% and single-stream channel estimate at pilot subcarriers.
% # The equalized symbols are demodulated and decoded to recover the PSDU.
% # The recovered PSDU is compared to the transmitted PSDU to determine if
% the packet has been recovered successfully.
%
% The simulation is run for the OFDMA configuration.

disp('Simulating OFDMA...');
throughputOFDMA = heMUSimulateScenario(cfgOFDMA,tgax,cfgSim);

% Sum throughput for all STAs and plot for all configurations
figure;
plot(cfgSim.Pathloss,sum(throughputOFDMA,2),'-x');
grid on;
xlabel('Pathloss (dB)');
ylabel('Throughput (Mbps)');
legend('OFDMA');
title('Raw AP Throughput');


end

