### Define which function to use
set l1vtxAlgo "sliding_window_vtx_algo_full"

# vertexing implementation
set l1vtxTopFunc bhv_find_pv

# reference implementation
set l1vtxRefFunc bhv_find_pv_ref

# set to zero to turn off C validation (runs but does not check against the reference implementation)
set l1vtxValidate 1

# set to zero to turn off track quality cuts
set l1tkQualityCuts 0

## version of the IP Core output
set l1vtxIPVersion 1.0