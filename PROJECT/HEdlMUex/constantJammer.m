% Add Comments Later


function [jammerOut] = constantJammer(dataLength, axCfg, idleTime)

jammerSig = randi([0 1], dataLength, 1, 'int8');
jammerOut = wlanWaveformGenerator(jammerSig, axCfg, 'IdleTime',idleTime*1e-6);

end
