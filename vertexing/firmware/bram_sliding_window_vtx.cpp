#include "bram_sliding_window_vtx.h"

void bhv_add_track(zbin_vt zbin, pt_t tkpt, ptsum_t hist[BHV_NBINS]) {
    #pragma HLS pipeline II=2
    #pragma HLS interface ap_memory port=hist
    #pragma HLS resource  variable=hist core=ram_1p
    if (zbin.valid) {
        pt_t pt = (tkpt >> 1);
        if (pt > BHV_MAXPT) pt = BHV_MAXPT;
        int sum = int(hist[zbin.bin])+pt;
        hist[zbin.bin] = (sum & (BHV_MAXBIN+1)) ? BHV_MAXBIN : sum;
    }
}

zbin_t bhv_find_pv(ptsum_t hist[BHV_NBINS], pt_t *sumpt) {
    #pragma HLS pipeline II=72
    #pragma HLS array_partition variable=hist complete dim=1
    #pragma HLS interface ap_memory port=hist
    #pragma HLS resource  variable=hist core=ram_1p
    zbin_t ibest = 0; int sbest = 0;
    int sigma_max = -999, zvtx_sliding = -999;
    for (unsigned int b = 0; b < BHV_NBINS-BHV_WINDOWSIZE; ++b) {
        int sigma = 0;
        for (unsigned int w = 0; w < BHV_WINDOWSIZE;  ++w) {
            sigma += hist[b+w];
        }
        if ( sigma > sigma_max) {
            sigma_max = sigma;
            for (unsigned int w = 0; w < BHV_WINDOWSIZE;  ++w) {
                zvtx_sliding += (hist[b+w]*bin_center_ref(b+w));
            }
            if(sigma!=0) zvtx_sliding /= sigma;
            else         zvtx_sliding = bin_center_ref(b+(BHV_WINDOWSIZE >> 1));
            ibest = b+(BHV_WINDOWSIZE >> 1);
            sbest = sigma;
        }
    }
    *sumpt = sbest;
    return ibest;
}

bool dummy(z0_t z0) { return (z0 > 0); }
