% Add Comments Later


function [jammerOut] = randomJammer(dataLength, axCfg, idleTime, header, numPackets)

packets = [1:numPackets]; 
packetsToJam = sort(randsample(packets, floor(0.3*10)));
if packets
    jammerSig = [header; randi([0 1], dataLength, 1, 'int8')];
    jammerOut = wlanWaveformGenerator(jammerSig, axCfg, 'IdleTime', idleTime*1e-6);
else
    jammerOut = zeros(dataLength); 
end