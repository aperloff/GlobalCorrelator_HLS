#include "../../firmware/data.h"
#include "../../regionizer/firmware/regionizer.h"

typedef ap_uint<9> ptsum_t;
typedef ap_uint<18> twoptsums_t;
typedef ap_uint<10> zbin_t;
struct zbin_vt {
    zbin_t bin;
    bool   valid;
};
// the current z0 scale is 20 counts per CM (1 bit = 0.5 mm)
// we divide down by 8 to have 1 bit = 0.4 cm. 
// with 71 bins, we cover +/- 14 cm, i.e. about +/- 3 sigma(Z)

// the standard pt scale is 0.25 GeV per unit
// if we max the track pt to 50 GeV, it would be 200 units
// we divide down by two
#define BHV_MAXPT 100
//#define BHV_MAXBIN 511
#define BHV_MAXBIN 511
#define BHV_CHI2MAX 100
#define BHV_PTMINTRA 2
#define BHV_NSTUBSMIN 4

#define BNV_SHIFT 3
#define BHV_NBINS 72
#define BHV_NHALFBINS (BHV_NBINS/2)
#define BHV_NSECTORS 2*N_IN_SECTORS
#define BHV_NTRACKS 18
#define BHV_WINDOWSIZE 3

//0 : truncation. Tracks with PT above PTMAX are ignored
//1 : saturation . Tracks with PT above PTMAX are set to PT=PTMAX.
#define BHV_TRUNCSAT 1


inline zbin_vt fetch_bin_ref(z0_t z0) {
    int zbin = (z0 >> BNV_SHIFT) + BHV_NHALFBINS; 
    bool valid = true;
    if (zbin < 0) { zbin = 0; valid = false; }
    if (zbin > BHV_NBINS-1) { zbin = 0; valid = false; }
    zbin_vt ret;
    ret.bin = zbin;
    ret.valid = valid;
    return ret;
}
inline z0_t bin_center_ref(zbin_t iz) {
    //std::cout<<"bin_center_ref::iz = " << iz << std::endl;
    //std::cout<<"bin_center_ref::BHV_NHALFBINS = " << BHV_NHALFBINS << std::endl;
    int z = int(iz) - BHV_NHALFBINS;
    //std::cout<<"bin_center_ref::(z << BNV_SHIFT) = " << (z << BNV_SHIFT) << std::endl;
    //std::cout<<"bin_center_ref::( 1 << (BNV_SHIFT-1) ) = " << ( 1 << (BNV_SHIFT-1) ) << std::endl;
    return z0_t((z << BNV_SHIFT) + ( 1 << (BNV_SHIFT-1) ));
}
inline ptsum_t window_sumpt(z0_t b, ptsum_t hist[BHV_NSECTORS][BHV_NBINS], ptsum_t (&binpt)[BHV_WINDOWSIZE]) {
    ptsum_t sigma = 0;
    for (unsigned int w = 0; w < BHV_WINDOWSIZE;  ++w) {
        for (unsigned int is = 0; is < BHV_NSECTORS; ++is) {    
            binpt[w] += hist[is][b+w];
        }
        sigma += binpt[w];
    }
    return sigma;
}
inline z0_t weighted_position_ref(zbin_t iz, ptsum_t binpt[BHV_WINDOWSIZE], pt_t sigma) {
    z0_t zvtx_sliding = 0;
    ap_int<21> zvtx_sliding_tmp = 0;
    //std::cout << "zvtx_sliding (start) = " << zvtx_sliding_tmp << std::endl;
    //std::cout << "weighted_position_ref::iz = " << iz << std::endl;
    //std::cout << "weighted_position_ref::zvtx_sliding (weighted sum) = " << std::flush;
    for(unsigned int w = 0; w < BHV_WINDOWSIZE; ++w) {
        //std::cout << "(" << binpt[w] << " * " << bin_center_ref(iz+w) << ")" << std::flush;
        //if(w<BHV_WINDOWSIZE-1) std::cout<< " + " << std::flush;
        zvtx_sliding_tmp += (binpt[w]*bin_center_ref(iz+w));
    }
    //std::cout << " = " << zvtx_sliding_tmp << std::endl;
    //std::cout << "weighted_position_ref:: sigma = " << sigma << std::endl;
    if(sigma!=0) zvtx_sliding_tmp /= sigma;
    else         zvtx_sliding_tmp = bin_center_ref(iz+(BHV_WINDOWSIZE >> 1));
    //std::cout << "weighted_position_ref::zvtx_sliding (final) = " << zvtx_sliding_tmp << std::endl;
    zvtx_sliding = zvtx_sliding_tmp;
    return zvtx_sliding;
}
inline bool track_quality_check_ref(TkObjExtended track) {
    // track quality cuts
    //if (fabs(z) > ZMAX ) continue; //ZMAX=25.
    if (track.hwChi2 > BHV_CHI2MAX) return false;
    if (track.hwPt < BHV_PTMINTRA) return false;
    if (track.hwStubs < BHV_NSTUBSMIN) return false;
    //if (nPS < nStubsPSmin) return false; //nStubsPSmin=3

    // quality cuts from Louise S, based on the pt-stub compatibility (June 20, 2014)
    //int trk_nstub  = (int) trackIter ->getStubRefs().size();
    //float chi2dof = chi2 / ((trk_nstub<<1)-4);
    //if(doPtComp)
    //float trk_consistency = trackIter ->getStubPtConsistency();
    //if (trk_nstub == 4) {
    //   if (fabs(eta)<2.2 && trk_consistency>10) return false;
    //   else if (fabs(eta)>2.2 && chi2dof>5.0) return false;
    //}
    //if(doTightChi2)
    //if (pt>10.0 && chi2dof>5.0) return false;

    return true;
}
void bhv_add_track(zbin_vt zbin, pt_t tkpt, ptsum_t hist[BHV_NBINS]) ;
zbin_t bhv_find_pv(ptsum_t hist[BHV_NBINS], pt_t *sumpt) ;
void bhv_find_pv_ref(TkObjExtended tracks[BHV_NSECTORS][BHV_NTRACKS], zbin_t & pvbin, z0_t & pv, int &pvsum) ;
bool dummy(z0_t z0) ;
