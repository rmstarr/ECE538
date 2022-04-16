% Add Comments Later


function [jammerOut] = deceptiveJammer(dataLength, axCfg, idleTime, header)

jammerSig = [header; randi([0 1], dataLength, 1, 'int8')];
jammerOut = wlanWaveformGenerator(jammerSig, axCfg, 'IdleTime', idleTime*1e-6);

end
