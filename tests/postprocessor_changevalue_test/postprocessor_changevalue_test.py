from options import *

#Postprocessor_changeValue_test = { INPUT : 'Postprocessor_changeValue_test.i',
#         EXODIFF : ['Postprocessor_changeValue_out.e']}

NucleationPostprocessor_test = { INPUT : 'NucleationPostprocessor_multi_test.i',
         EXODIFF : ['NucleationPostprocessor_multi_out.e']}

NucleationPostprocessor_parallel_test = { INPUT : 'NucleationPostprocessor_multi_test.i',
         EXODIFF : ['NucleationPostprocessor_multi_out.e'],
         PREREQ : ['NucleationPostprocessor_test'],
         MIN_PARALLEL : 2}

OneSeed_test = { INPUT : 'OneSeed_test.i',
         EXODIFF : ['OneSeed_out.e']}

OneSeed_restart_test = { INPUT : 'OneSeed_restart_test.i',
         EXODIFF : ['OneSeed_restart_out.e']}
