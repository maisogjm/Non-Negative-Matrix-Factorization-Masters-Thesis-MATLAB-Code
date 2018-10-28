function [ simdat ] = SimulateData_KD_N(trueK)

switch(trueK)
    case (3)
        simdat = SimulateData_KD3();
    case (4)
        simdat = SimulateData_KD4();
    case (5)
        simdat = SimulateData_KD5();
    case (6)
        simdat = SimulateData_KD6();
    case (8)
        simdat = SimulateData_KD8();
    case (10)
        simdat = SimulateData_KD10();
    case (12)
        simdat = SimulateData_KD12();
    case (14)
        simdat = SimulateData_KD14();
    case (16)
        simdat = SimulateData_KD16();
    otherwise
        error('SimulateData_KD_N: %d not implemented',trueK)
end
