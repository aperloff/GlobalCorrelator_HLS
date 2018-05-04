#include "firmware/bram_sliding_window_vtx.h"
#if (defined(__GXX_EXPERIMENTAL_CXX0X__) or defined(CMSSW))
    #include <algorithm>
    #include <iterator>
    #if DEBUG == 3
        #include "firmware/utility.h"
    #endif
#endif

void bhv_find_pv_ref(TkObjExtended tracks[BHV_NSECTORS][BHV_NTRACKS], zbin_t & pvbin, z0_t & pv, int &pvsum) {
    ptsum_t histos[BHV_NSECTORS][BHV_NBINS];
    for (unsigned int is = 0; is < BHV_NSECTORS; ++is) {
        for (unsigned int b = 0; b < BHV_NBINS; ++b) {
            histos[is][b] = 0;
        }
        for (unsigned int it = 0; it < BHV_NTRACKS; ++it) {
            pt_t tkpt = (tracks[is][it].hwPt >> 1);
            //if (track_quality_check_ref(tracks[is][it])) {
                //
                // Saturation or truncation
                //
                if (tkpt > BHV_MAXPT) {
                    tkpt = (BHV_TRUNCSAT) ? BHV_MAXPT : 0;
                } 

                //
                // Check bin validity of bin found for the current track
                //
                zbin_vt bin = fetch_bin_ref(tracks[is][it].hwZ0);
                assert(bin.bin >= 0 && bin.bin < BHV_NBINS);

                //
                // If the bin is valid then sum the tracks
                //
                if (bin.valid) {
                    int newsum = int(histos[is][bin.bin]) + tkpt;
                    histos[is][bin.bin] = ptsum_t(newsum > BHV_MAXBIN ? BHV_MAXBIN : newsum);
                    #if(DEBUG==3)
                        if (bin.bin==37) {
                            std::cout<<"bhv_find_pv_ref::is="<<is<<"\tcurrent track pT=" << tkpt << "\tcurrent sum=" <<histos[is][bin.bin] << std::endl;
                        }
                    #endif
                }
            //}
        }
    }

    //
    // Setup the intermediate and final storage objects
    //
    unsigned int b_max = 0, sigma_max = 0;
    ptsum_t binpt[BHV_WINDOWSIZE], binpt_max[BHV_WINDOWSIZE];

    //
    // Loop through all bins, taking into account the fact that the last bin is
    //  nbins-window_width+1
    // Save the small amount of intermediate information if the sum pT (sigma) of
    //  the current bin is the highest so far
    //
    for (unsigned int b = 0; b < BHV_NBINS-BHV_WINDOWSIZE+1; ++b) {
        int sigma = 0;
        for (unsigned int w = 0; w < BHV_WINDOWSIZE;  ++w) {
            binpt[w] = 0;
            for (unsigned int is = 0; is < BHV_NSECTORS; ++is) {    
                binpt[w] += histos[is][b+w];
            }
            sigma += binpt[w];
        }
        #if (defined(__GXX_EXPERIMENTAL_CXX0X__) or defined(CMSSW)) && (DEBUG == 3)
            std::cout << "bin = " << b << std::endl;
            show<ptsum_t,BHV_WINDOWSIZE>(binpt);
        #endif
        if (sigma > sigma_max) {
            sigma_max = sigma;
            b_max = b;
            std::copy(std::begin(binpt), std::end(binpt), std::begin(binpt_max));
        }
    }
    //
    //Find the weighted position only for the highest sum pT window
    // 
    //vivado_hls is smart enough to inline this function right here
    //That is why we can split off this simple line of code
    bin_plus_half_window_ref<unsigned int,zbin_t>(b_max,pvbin);
    pvsum = sigma_max;
    pv = weighted_position_ref(b_max,binpt_max,sigma_max);
}