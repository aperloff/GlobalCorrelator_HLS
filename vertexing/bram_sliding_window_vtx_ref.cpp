#include "firmware/bram_sliding_window_vtx.h"
#include <iostream>

void bhv_find_pv_ref(TkObjExtended tracks[BHV_NSECTORS][BHV_NTRACKS], zbin_t & pvbin, z0_t & pv, int &pvsum) {
    ptsum_t histos[BHV_NSECTORS][BHV_NBINS];
    for (unsigned int is = 0; is < BHV_NSECTORS; ++is) {
        for (unsigned int b = 0; b < BHV_NBINS; ++b) {
            histos[is][b] = 0;
        }
        for (unsigned int it = 0; it < BHV_NTRACKS; ++it) {
            pt_t tkpt = (tracks[is][it].hwPt >> 1);
            if (track_quality_check_ref(tracks[is][it])) {
                // saturation or truncation
                if (tkpt > BHV_MAXPT) {
                    tkpt = (BHV_TRUNCSAT) ? BHV_MAXPT : 0;
                } 
                zbin_vt bin = fetch_bin_ref(tracks[is][it].hwZ0);
                assert(bin.bin >= 0 && bin.bin < BHV_NBINS);
                if (bin.valid) {
                    int newsum = int(histos[is][bin.bin]) + tkpt;
                    histos[is][bin.bin] = ptsum_t(newsum > BHV_MAXBIN ? BHV_MAXBIN : newsum);
                    //if(bin.bin==37) std::cout<<"bhv_find_pv_ref::is="<<is<<"\tcurrent track pT=" << tkpt << "\tcurrent sum=" <<histos[is][bin.bin] << std::endl;
                }
            }
        }
    }

    int ibest = 0, sbest = 0;
    int sigma_max = -999, zvtx_sliding = -999;
    ptsum_t binpt[BHV_WINDOWSIZE];
    for (unsigned int b = 0; b < BHV_NBINS-BHV_WINDOWSIZE; ++b) {
        int sigma = 0;
        for (unsigned int w = 0; w < BHV_WINDOWSIZE;  ++w) {
            binpt[w] = 0;
            for (unsigned int is = 0; is < BHV_NSECTORS; ++is) {    
                binpt[w] += histos[is][b+w];
            }
            sigma += binpt[w];
        }
        //std::cout << "bin = " << b << std::endl;
        //show<ptsum_t,BHV_WINDOWSIZE>(binpt);
        if ( sigma > sigma_max) {
            sigma_max = sigma;
            zvtx_sliding = weighted_position_ref(b,binpt,sigma);
            //vivado_hls is smart enough to inline this function right here
            //That is why we can split off this simple line of code
            bin_plus_half_window_ref<int,int>(b,ibest);
            sbest = sigma;
        }
    }
    pvbin = ibest;
    pvsum = sbest;
    pv = zvtx_sliding;
}