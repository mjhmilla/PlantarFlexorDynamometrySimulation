function benchRecordSet = insertResultsIntoSet(singleBenchRecord, ...
                                benchRecordSet, simulationNumber)


simNo=simulationNumber;

benchRecordSet.activation(:,simNo)                  = singleBenchRecord.activation(:,1)                 ;
benchRecordSet.cpuTime(:,simNo)                     = singleBenchRecord.cpuTime(:,1)                    ;
benchRecordSet.normFiberForceAlongTendon(:,simNo)   = singleBenchRecord.normFiberForceAlongTendon(:,1)  ;
benchRecordSet.normFiberLength(:,simNo)             = singleBenchRecord.normFiberLength(:,1)            ;    
benchRecordSet.pennationAngle(:,simNo)              = singleBenchRecord.pennationAngle(:,1)             ;
benchRecordSet.normFiberVelocity(:,simNo)           = singleBenchRecord.normFiberVelocity(:,1)          ;
benchRecordSet.pennationAngVelocity(:,simNo)        = singleBenchRecord.pennationAngVelocity(:,1)       ;
benchRecordSet.fiberStiffnessAlongTendon(:,simNo)   = singleBenchRecord.fiberStiffnessAlongTendon(:,1)  ;
benchRecordSet.tendonStiffnessAlongTendon(:,simNo)  = singleBenchRecord.tendonStiffnessAlongTendon(:,1) ;
benchRecordSet.muscleStiffness(:,simNo)             = singleBenchRecord.muscleStiffness(:,1)            ;
benchRecordSet.fiberVelocity(:,simNo)               = singleBenchRecord.fiberVelocity(:,1)              ;
benchRecordSet.fiberVelocityAlongTendon(:,simNo)    = singleBenchRecord.fiberVelocityAlongTendon(:,1)   ;
benchRecordSet.tendonVelocity(:,simNo)              = singleBenchRecord.tendonVelocity(:,1)             ;
benchRecordSet.pathLength(:,simNo)                  = singleBenchRecord.pathLength(:,1)                 ;
benchRecordSet.pathVelocity(:,simNo)                = singleBenchRecord.pathVelocity(:,1)               ;
benchRecordSet.dSystemEnergyLessWork(:,simNo)       = singleBenchRecord.dSystemEnergyLessWork(:,1)      ;
benchRecordSet.systemEnergyLessWork(:,simNo)        = singleBenchRecord.systemEnergyLessWork(:,1)       ;
benchRecordSet.tendonPotentialEnergy(:,simNo)       = singleBenchRecord.tendonPotentialEnergy(:,1)      ;
benchRecordSet.fiberPotentialEnergy(:,simNo)        = singleBenchRecord.fiberPotentialEnergy(:,1)       ;
benchRecordSet.fiberActiveWork(:,simNo)             = singleBenchRecord.fiberActiveWork(:,1)            ;
benchRecordSet.dampingWork(:,simNo)                 = singleBenchRecord.dampingWork(:,1)                ;
benchRecordSet.boundaryWork(:,simNo)                = singleBenchRecord.boundaryWork(:,1)               ;
benchRecordSet.tendonPower(:,simNo)                 = singleBenchRecord.tendonPower(:,1)                ;
benchRecordSet.fiberParallelElementPower(:,simNo)   = singleBenchRecord.fiberParallelElementPower(:,1)  ;
benchRecordSet.fiberActivePower(:,simNo)            = singleBenchRecord.fiberActivePower(:,1)           ;
benchRecordSet.dampingPower(:,simNo)                = singleBenchRecord.dampingPower(:,1)               ;
benchRecordSet.boundaryPower(:,simNo)               = singleBenchRecord.boundaryPower(:,1)              ;  

