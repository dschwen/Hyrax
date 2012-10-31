from options import *

Value_test = { INPUT : 'Value_test.i',
         EXODIFF : ['testValue.e'],
         PETSC_VERSION : ['>=3.1.0']}

ValueForcing_test = { INPUT : 'ValueForcing_test.i',
         EXODIFF : ['testValueForcing.e'],
         PETSC_VERSION : ['>=3.1.0'],
         ABS_ZERO: 1e-8}

#ValueNucleation_test = { INPUT : 'ValueNucleation_test.i',
#         EXODIFF : ['testValueNucleation.e'],
#         PETSC_VERSION : ['>=3.1.0'],
#         ABS_ZERO: 1e-8}

#ValueNucleationForcing_test = { INPUT : 'ValueNucleationForcing_test.i',
#         EXODIFF : ['testValueNucleationForcing.e'],
#         PETSC_VERSION : ['>=3.1.0']}

#ACNucleationForcing_test = { INPUT : 'ACNucleationForcing_test.i',
#         EXODIFF : ['testACNucleationForcing.e'],
#         PETSC_VERSION : ['>=3.1.0']}

