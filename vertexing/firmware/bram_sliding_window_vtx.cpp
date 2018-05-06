#include "bram_sliding_window_vtx.h"
#if (defined(__GXX_EXPERIMENTAL_CXX0X__) and defined(CMSSW))
    #include "utility.h"
#endif

///////////////////////////
// -- BHV_FIND_PV FUNCTION
void bhv_find_pv(ptsum_t hist[BHV_NBINS], zbin_t *pvbin_hw, z0_t *pv_hw, pt_t *sumpt) {
//void bhv_find_pv(TkObjExtended tracks[BHV_NSECTORS][BHV_NTRACKS], zbin_t *pvbin_hw, z0_t *pv_hw, pt_t *sumpt) {
    #pragma HLS pipeline II=1

    /*
    #pragma HLS array_partition variable=tracks complete
    ptsum_t hist[BHV_NBINS];
    for (unsigned int b = 0; b < BHV_NBINS; ++b) {
        #pragma HLS unroll
        hist[b] = 0;
    }
    bhv_merge_sectors_and_tracks(tracks, hist, QUALITY);
    */

    //Without the following pragma the interface pragma will work
    //Only really makes sense if unrolling or using multidimensional array
    #pragma HLS array_partition variable=hist complete dim=1
    //These increase LUT resource usage slightly (not sure why needed then)
    //ug902-vivado-high-level-synthesis.pdf page 459
    //#pragma HLS interface bram port=hist //ap_memory
    //ug902-vivado-high-level-synthesis.pdf page 184
    //#pragma HLS resource  variable=hist core=ram_2p//RAM_T2P_BRAM

    //
    // Setup the intermediate and final storage objects
    //
    zbin_t ibest = 0, b_max = 0;
    slidingsum_t sigma_max = 0;//, pt_weighted_position = 0;
    zsliding_t zvtx_sliding;
    ptsum_t binpt_max[BHV_WINDOWSIZE];
    #pragma HLS array_partition variable=binpt_max complete dim=1
    slidingsum_t hist_window_sums[BHV_NSUMS];
    #pragma HLS array_partition variable=hist_window_sums complete dim=1
    INITIALIZELOOP: for (unsigned int b = 0; b < BHV_NBINS-BHV_WINDOWSIZE+1; ++b) {
        hist_window_sums[b] = 0;
    }

    #if (defined(__GXX_EXPERIMENTAL_CXX0X__) and defined(CMSSW))
        show<ptsum_t,BHV_NBINS>(hist,80,0,-1,std::cout,"","\e[92m");
    #endif
    bhv_compute_sums(hist,hist_window_sums);

    //#if (defined(__GXX_EXPERIMENTAL_CXX0X__) or defined(CMSSW))
    //    show<slidingsum_t,BHV_NBINS-BHV_WINDOWSIZE+1>(hist_window_sums);
    //#endif

    // Search through all of the bins to find the one with the highest sum pT
    //  within the window
    //bhv_find_max(hist_window_sums,b_max,sigma_max);
    bhv_parallel_find_max<slidingsum_t,zbin_t,BHV_NSUMS>(hist_window_sums,b_max,sigma_max);

    //
    // Save the pT information for just the highest pT window
    //
    WINDOWSELECTLOOP: for (unsigned int w = 0; w < BHV_WINDOWSIZE;  ++w) {
        binpt_max[w] = hist[b_max+w];
        //binpt_max[w] = bhv_access_array<zbin_t,ptsum_t,BHV_NBINS>(hist,b_max,zbin_t(w));
    }

    *sumpt = sigma_max;
    zvtx_sliding = weighted_position_ref(b_max,binpt_max,sigma_max);
    *pv_hw = zvtx_sliding;
    //vivado_hls is smart enough to inline this function right here
    //That is why we can split off this simple line of code
    bin_plus_half_window_ref<zbin_t,zbin_t>(b_max,ibest);
    *pvbin_hw = ibest;

    return;
}

//bool dummy(z0_t z0) { return (z0 > 0); }
