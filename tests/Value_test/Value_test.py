from options import *

Value_test = { INPUT : 'Value_test.i',
         EXODIFF : ['testValue.e']}

ValueForcing_test = { INPUT : 'ValueForcing_test.i',
         EXODIFF : ['testValueForcing.e'],
         ABS_ZERO: 1e-8}

ValueNucleation_test = { INPUT : 'ValueNucleation_test.i',
         EXODIFF : ['testValueNucleation.e'],
         ABS_ZERO: 1e-8}

ValueNucleationForcing_test = { INPUT : 'ValueNucleationForcing_test.i',
         EXODIFF : ['testValueNucleationForcing.e']}

ACNucleationForcing_test = { INPUT : 'ACNucleationForcing_test.i',
         EXODIFF : ['testACNucleationForcing.e']}
