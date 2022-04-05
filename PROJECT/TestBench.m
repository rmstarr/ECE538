%% Test Bench
% Area to perform tests and analyses of results of our designed 
% communication architecture and jammers

%% Creating the 802.11ax Network

% Assuming that I am using an OFDMA-created signal as done in the
% heMUDownlinkMUExample.m file 
numTx = 6;                          % Number of TX antennas
guardInterval = 0.8;                % Specified in microseconds
axOFDMA = wlanHEMUConfig(112);

% Configure per user parameters -- pretty much replicating example
axOFDMA.NumTransmitAntennas = numTx;
axOFDMA.GuardInterval = guardInterval;

% Note that I am using QPSK with coding rate = 1/2 (MCS=1) for each of the
% stations in this scenario

% STA#1
axOFDMA.User{1}.NumSpaceTimeStreams = 2;
axOFDMA.User{1}.MCS = 1;
axOFDMA.User{1}.APEPLength = 1000;

% STA#2
axOFDMA.User{2}.NumSpaceTimeStreams = 2;
axOFDMA.User{2}.MCS = 1;
axOFDMA.User{2}.APEPLength = 1000;

% STA#3
axOFDMA.User{3}.NumSpaceTimeStreams = 2;
axOFDMA.User{3}.MCS = 1;
axOFDMA.User{3}.APEPLength = 1000;

% STA#4
axOFDMA.User{4}.NumSpaceTimeStreams = 2;
axOFDMA.User{4}.MCS = 1;
axOFDMA.User{4}.APEPLength = 1000;

% Simulate scenario uses wLanWaveformGenerator function (L77), which has a
% parameter to specify the delay time. However, if you don't modify it, its
% default value is 0 -- seems like that is the best option to create a
% jammer is to make the idleTime equal to 0 seconds

% They create their own spatial mapping -- we can use a SpatialMapping that
% is 'Direct' or 'Fourier' if we prefer/cannot learn the SpatialMapping
% concept quickly enough for our implementation

% Not sure how to "inject" in the jammer. At best, I think this would be a
% case where we have to run a sample environment a bunch of times and
% randomly add in the jammer -- we can't add it in once the loop gets
% running. Can use this idea to our advantage however -- at the end of
% every single loop, the algorithm can take a snapshot of the metrics and
% use the results to make an assessment on the health of the network. This
% assessment can lead it to go through its decision tree and determine
% whether we were jammed. Therefore, the results won't be in super-real
% time. The results would only occur after the simulation loop is completed

% Looks like we'll pretty much need to use most of the heMUSimulateScenario
% example, but return more metrics and add in a way to inject a jammer. We
% can inject a jammer at random using the switch case decisions that I have
% below. Might want to add in a case 0 where we simply skip the jamming
% (that way we, and our algorithm aren't biased into always expecting a
% jamming event to occur). We can make an array of numbers 0-5 and set a
% certain percentage of the array to be filled with 0's. Then, in the
% switch case, a random number is selected. If 0, no jamming takes place,
% else, select from one of the mapped jammers.

% Note that when they're calculating their amplitude that results from the
% expected path loss, they're assuming an alpha value of 2 (as they're
% doing 10^((transmitPower - 30) / 20)

% Pulling metrics:
% PHY:
%   - SNR - should be able to take the output of the tx variable (line 77,
%   where the wlanWaveformGenerator function is used. Can then call
%   function to get the SNR of the signal and see if significant dips
%   occur, which could suggest that the signal is getting jamming
%   - RSS - derived from the path loss? Therefore, it seems like we can't
%   use this metric because the function sets values to use, it doesn't
%   derive them based on distances / environment loss conditions
% Link:
%   - PER - calculated by simulate scenario function, just need to set as
%   output
%   - Throughput - currently an output of the simulate scenario function
% Network:
%   - Don't see any real ways that the network layer is used in the example
%   scenario. Do we need to add our own "field" for an IP address (which we
%   can specify as a string). Then, when the network is created, we have
%   some good nodes and some jammer nodes that come and go (that way we can
%   introduce new IP addresses into the network, and we're not 100% sure if
%   they're good nodes or if they're jammers)


%% Selecting Between Jammers
% Planning to do 5 different jammers (assume we assign them the following):
% Constant Jammer = 1
% Deceptive Jammer = 2
% Random Jammer = 3 
% RTS/CTS Jammer = 4
% Data/Ack Jammer = 5
% "Good" node = 6 (serves as a control node)

% Create a control block to randomly select between the various jammers
possibleScenarioConfigs = 6;

% Way to select between each jammer randomly, one time
jammerNumUniqueCall = randperm(possibleScenarioConfigs); % Generates number 1-5 one time only

% Iterate over the array and see what the order of the jammers is
disp('Doing the unique test run of calling all jammers once')
for i=1:length(jammerNumUniqueCall)
    % Pull the index out of jammerNumUniqueCall
    jammerNum = jammerNumUniqueCall(i);
    
    % Select the jammer to use
    switch jammerNum
        case 1
            disp('Selected Constant Jammer');
            % Create a barrage jammer object and add it to disrupt network
            % Will show the increase in the noise, which decreases SNR and
            % throughput, increases PER -- should be the easiest to detect
        case 2
            disp('Selected Deceptive Jammer');
        case 3
            disp('Selected Random Jammer');
            % Inside the method, give it a 50-50 shot to model itself after
            % either the barrage jammer or the deceptive jammer
            % Also, specify a sleep period. Can then do one of two things:
            % Say we simulate for 10 packets. If the sleep window is 0.2,
            % select 2 indices to skip (i.e. drop out on iterations 3 and
            % 4, for example)
        case 4
            disp('Selected RTS/CTS Jammer');
        case 5
            disp('Selected Data/ACK Jammer');
        case 6
            % 
            disp('Add just a normal node (ensure we do not make false positives');
        otherwise
            disp('An error occurred');
            
    end
end

% Select between jammers in a random list (can call jammers in any order,
% can call jammers more times than others)
jammerNumRandomCall = randi([1, possibleScenarioConfigs], 1, 10);

% Iterate over this array and see what the order of the jammers is
disp('Doing the random selection of the jammers')
for j=1:length(jammerNumRandomCall)
    % Pull the index out of jammerNumUniqueCall
    jammerNum = jammerNumRandomCall(j);
    
    % Select the jammer to use
    switch jammerNum
        case 1
            disp('Selected Constant Jammer')
        case 2
            disp('Selected Deceptive Jammer')
        case 3
            disp('Selected Random Jammer')
        case 4
            disp('Selected RTS/CTS Jammer')
        case 5
            disp('Selected Data/ACK Jammer') 
        otherwise
            disp('An error occurred')
    end
end



%% Perform Simulations


% With our network configured, we'll want to do one initial run of our
% network. We will collect some initial metrics for the algorithm to know
% what the "normal" performance metrics of the channel looks like

[initRxSNR, initPacketErr, initThroughput] = ...
    simulateNetworkComms(ieeeCfg, chnMdl, ntwkSim);

initialResults = [initRxSNR, initPacketErr, initThroughput];


% We can determine how many different iterations (or scenarios) we want to
% loop through. For each scenario, we'll select from one of our jammers (or
% we'll select just a normal node). Then, we'll perform our simulation with
% the given configuration. We will also likely want to pass the jammer
% number in as an input

% Result will end up looking like the following (assume we run 20 different
% scenarios)
numTestScenarios = 20;
possibleJammerConfigs = 6;

% Generate the array of jammers / good points (use code above)
scenarioCfgs = randi([1, possibleJammerConfigs], 1, numTestScenarios);

for i=1:length(scenarioCfgs)
    % Get our jammer configuration
    jammer = scenarioCfgs(i);
    % Select the jammer to use
    switch jammer
        % Configure to use a constant jammer
        case 1
            disp('Selected Constant Jammer')
            constantJammer();
        % Configure to use a deceptive jammer
        case 2
            disp('Selected Deceptive Jammer')
        % Configure to use a random jammer
        case 3
            disp('Selected Random Jammer')
        % Configure to use an RTS/CTS jammer
        case 4
            disp('Selected RTS/CTS Jammer')
        % Configure to use a Data/ACK jammer
        case 5
            disp('Selected Data/ACK Jammer') 
        % Add a normal node to the configuration
        case 6
            disp('Adding a normal node as a control experiment. Not all nodes are bad!')
        otherwise
            disp('An error occurred. Please debug the code.')
    end
    
    % With the configuration selected, run the simulation
    [rxSNR, packetErr, throughput] = ...
        simulateNetworkComms(ieeeCfg, chnMdl, ntwkSim);
    
    % Now, do we want to save the results to an array, or do we want to
    % call our "algorithm" to assess if we're getting jammed right here?
    % Likely the latter event. In which case, we'll make a call to the
    % algorithm here -- do we want it to just detect that we're jamming, or
    % do we also want to assess what type of jamming is being performed?
    % NOTE: There should be two inputs here:
    % - 1) The initial data results that we got from just our normal nodes
    % - 2) The results that we got from the simulation that was just
    % completed.
    [jamminngStatus, jammingType] = ...
        jammingDetectionDecision(initialResults, [rxSNR, packetErr, throughput]);
    
end

% Maybe we could do a fun little "timing diagram" that shows the results,
% where a "red line" shows the jamming scenarios that we used. A blue line
% can then be used to show the detection decision that the algorithm
% reached

% Could do two plots: One plot that is just "Did we detect the jamming
% event?" and one plot for "did we detect which jamming event occurred?"

