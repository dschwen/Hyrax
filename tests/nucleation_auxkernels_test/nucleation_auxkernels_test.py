from options import *

AuxNucleationRate_refine0_test = { INPUT : 'AuxNucleationRate_refine0.i',
         EXODIFF : ['AuxNucleationRate_refine0.e'],
         PETSC_VERSION : ['>=3.1.0']}


AuxNucleationRate_refine1_test = { INPUT : 'AuxNucleationRate_refine1.i',
         EXODIFF : ['AuxNucleationRate_refine1.e'],
         PETSC_VERSION : ['>=3.1.0']}

AuxSupersaturation_test = { INPUT : 'AuxSupersaturation.i',
         EXODIFF : ['AuxSupersaturation.e'] }

AuxNucleationProbability_test = { INPUT : 'AuxNucleationProbability.i',
         EXODIFF : ['AuxNucleationProbability.e'] }
