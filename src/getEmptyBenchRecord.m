function benchRecord = getEmptyBenchRecord(numberOfPoints, numberOfSimulations)

benchRecord = [];
benchRecord.activation                  = zeros(numberOfPoints,numberOfSimulations);
benchRecord.time                        = zeros(numberOfPoints,numberOfSimulations);
benchRecord.cpuTime                     = zeros(numberOfPoints,numberOfSimulations);
benchRecord.normFiberForceAlongTendon   = zeros(numberOfPoints,numberOfSimulations);
benchRecord.normFiberLength             = zeros(numberOfPoints,numberOfSimulations);    
benchRecord.pennationAngle              = zeros(numberOfPoints,numberOfSimulations);
benchRecord.normFiberVelocity           = zeros(numberOfPoints,numberOfSimulations);
benchRecord.pennationAngVelocity        = zeros(numberOfPoints,numberOfSimulations);
benchRecord.fiberStiffnessAlongTendon   = zeros(numberOfPoints,numberOfSimulations);
benchRecord.tendonStiffnessAlongTendon  = zeros(numberOfPoints,numberOfSimulations);
benchRecord.muscleStiffness             = zeros(numberOfPoints,numberOfSimulations);
benchRecord.fiberVelocity               = zeros(numberOfPoints,numberOfSimulations);
benchRecord.fiberVelocityAlongTendon    = zeros(numberOfPoints,numberOfSimulations);
benchRecord.tendonVelocity              = zeros(numberOfPoints,numberOfSimulations);
benchRecord.pathLength                  = zeros(numberOfPoints,numberOfSimulations);
benchRecord.pathVelocity                = zeros(numberOfPoints,numberOfSimulations);

benchRecord.dSystemEnergyLessWork       = zeros(numberOfPoints,numberOfSimulations);
benchRecord.systemEnergyLessWork        = zeros(numberOfPoints,numberOfSimulations);
benchRecord.tendonPotentialEnergy       = zeros(numberOfPoints,numberOfSimulations);
benchRecord.fiberPotentialEnergy        = zeros(numberOfPoints,numberOfSimulations);
benchRecord.fiberActiveWork             = zeros(numberOfPoints,numberOfSimulations);
benchRecord.dampingWork                 = zeros(numberOfPoints,numberOfSimulations);
benchRecord.boundaryWork                = zeros(numberOfPoints,numberOfSimulations);

benchRecord.tendonPower                 = zeros(numberOfPoints,numberOfSimulations);
benchRecord.fiberParallelElementPower   = zeros(numberOfPoints,numberOfSimulations);
benchRecord.fiberActivePower            = zeros(numberOfPoints,numberOfSimulations);
benchRecord.dampingPower                = zeros(numberOfPoints,numberOfSimulations);
benchRecord.boundaryPower               = zeros(numberOfPoints,numberOfSimulations);  

benchRecord.activeFiberForce             = zeros(numberOfPoints,numberOfSimulations);  
benchRecord.activeFiberForceAlongTendon  = zeros(numberOfPoints,numberOfSimulations);  
benchRecord.passiveFiberForce            = zeros(numberOfPoints,numberOfSimulations); 
benchRecord.passiveFiberForceAlongTendon = zeros(numberOfPoints,numberOfSimulations); 

benchRecord.fiberActiveForceLengthMultiplier  = zeros(numberOfPoints,numberOfSimulations);  
benchRecord.fiberPassiveForceLengthMultiplier = zeros(numberOfPoints,numberOfSimulations);  
benchRecord.fiberForceVelocityMultiplier      = zeros(numberOfPoints,numberOfSimulations);  
benchRecord.normDamping                  = zeros(numberOfPoints,numberOfSimulations);  