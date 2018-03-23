#include "bram_sliding_window_vtx.h"

///////////////////////////
// -- BHV_ADD_TRACK FUNCTION
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

///////////////////////////
// -- BHV_FIND_PV FUNCTION
zbin_t bhv_find_pv(ptsum_t hist[BHV_NBINS], pt_t *sumpt, z0_t *pv_hw) {
    #pragma HLS pipeline II=36
    //Without the following pragma the interface pragma will work
    //Only really makes sense if unrolling or using multidimensional array
    //#pragma HLS array_partition variable=hist complete dim=1
    //#pragma HLS array_partition variable=hist cyclic factor=4 dim=1
    //These decrease resource usage, but incease interval
    //#pragma HLS interface ap_memory port=hist
    #pragma HLS interface bram port=hist
    #pragma HLS resource  variable=hist core=ram_2p
    zbin_t ibest = 0; int sbest = 0;
    slidingsum_t sigma_max = 0;//, pt_weighted_position = 0;
    zsliding_t zvtx_sliding;
    ptsum_t binpt[BHV_WINDOWSIZE];
    BINLOOP: for (unsigned int b = 0; b < BHV_NBINS-BHV_WINDOWSIZE+1; ++b) {
        slidingsum_t sigma = 0;
        SIGMASUMLOOP: for (unsigned int w = 0; w < BHV_WINDOWSIZE;  ++w) {
            sigma += hist[b+w];
            binpt[w] = hist[b+w];
        }
        if ( sigma > sigma_max) {
            sigma_max = sigma;
            //for (unsigned int w = 0; w < BHV_WINDOWSIZE;  ++w) {
            //    pt_weighted_position += (hist[b+w]*bin_center_ref(zbin_t(b+w)));
            //}
            //if(sigma!=0) {
            //    inverse_t inv;
            //    inversion<slidingsum_t, inverse_t>(sigma, inv);
            //    zvtx_sliding = pt_weighted_position*inv; //LUT replacement for zvtx_sliding /= sigma;
            //}
            //else {
            //    zvtx_sliding = (BHV_WINDOWSIZE >> 1) + ((BHV_WINDOWSIZE%2!=0) ? 0.5 : 0.0);//bin_center_ref(zbin_t(b+(BHV_WINDOWSIZE >> 1)));
            //}
            //zvtx_sliding += b;
            //zvtx_sliding = bin_center_ref(zvtx_sliding);
            zvtx_sliding = weighted_position_ref(b,binpt,sigma);
            //vivado_hls is smart enough to inline this function right here
            //That is why we can split off this simple line of code
            bin_plus_half_window_ref<zbin_t,zbin_t>(zbin_t(b),ibest);
            sbest = sigma;
        }
    }
    *sumpt = sbest;
    *pv_hw = zvtx_sliding;
    return ibest;
}

//bool dummy(z0_t z0) { return (z0 > 0); }
