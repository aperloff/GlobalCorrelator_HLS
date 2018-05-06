# get the configuration
source config_hls_sliding_window_vtx.tcl

# open the project, don't forget to reset
open_project -reset proj_bswv
#set_top bhv_add_track
set_top ${l1vtxTopFunc}
add_files firmware/bram_sliding_window_vtx.cpp -cflags "-DTESTBOARD -std=c++0x -DQUALITY=${l1tkQualityCuts}"
# -DCMSSW
add_files -tb bram_sliding_window_vtx_test.cpp  -cflags "-DTESTBOARD -std=c++0x -DTOP_FUNC=${l1vtxTopFunc} -DREF_FUNC=${l1vtxRefFunc} -DVALIDATE=${l1vtxValidate} -DQUALITY=${l1tkQualityCuts}"
add_files -tb bram_sliding_window_vtx_ref.cpp -cflags "-DTESTBOARD -std=c++0x -DQUALITY=${l1tkQualityCuts}"
add_files -tb ../utils/pattern_serializer.cpp -cflags "-DTESTBOARD"
add_files -tb ../utils/DiscretePFInputs_IO.h -cflags "-DTESTBOARD -std=c++0x"
add_files -tb ../regionizer/data/barrel_sectors_1x12_TTbar_PU140.dump
add_files -tb ../regionizer/data/barrel_alltracks_sectors_1x12_TTbar_PU140.dump
#Can use -DCPP0X in cflags and then check for existence of CPP0X in C++ code


# reset the solution
#open_solution -reset "solution_sliding_parallel_find_max"
open_solution -reset "solution_sliding_debug"
#set_part {xc7k160tfbg484-2} -tool vivado
#set_part {xcku9p-ffve900-2-i-EVAL}
#set_part {xc7vx690tffg1927-2}
#set_part {xcku115-flvf1924-2-i}
set_part {xcvu9p-flga2104-2-i-EVAL}
#200MHz
#create_clock -period 5 -name default
#240MHz
#create_clock -period 4.16667 -name default
#set_clock_uncertainty 1.5
#320MHz
#320MHz might be a problem. Need to check the RTL implementation
create_clock -period 3.125 -name default

config_interface -trim_dangling_port
config_core DSP48 -latency 2

# do stuff
csim_design
csynth_design
#cosim_design -trace_level all
#export_design -format ip_catalog -vendor "cern-cms" -version ${l1vtxIPVersion} -description "${l1vtxTopFunc}"

# exit Vivado HLS
exit
