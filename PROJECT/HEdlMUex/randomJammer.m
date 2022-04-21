function [jammerOut] = randomJammer(dataLength, axCfg, idleTime, header, dutyCycle)
jamVals = 0:2;
probs = [(1-dutyCycle) (dutyCycle/2) (dutyCycle/2)];
selectOne = 1;
rng = randsample(jamVals, selectOne, true, probs);

switch rng
    case 0
        jammerOut = zeros(1, dataLength);
    case 1
        jammerOut = constantJammer(dataLength, axCfg, idleTime);
    case 2
        jammerOut = deceptiveJammer(dataLength, axCfg, idleTime, header);
    otherwise
        disp('Error in return value')
end 

end
