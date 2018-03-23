#include "../../firmware/data.h"
#include "../../regionizer/firmware/regionizer.h"
#include "utility.h"

#define DEBUG 0

typedef ap_uint<9> ptsum_t;
typedef ap_uint<18> twoptsums_t;
typedef ap_uint<10> zbin_t;
struct zbin_vt {
    zbin_t bin;
    bool   valid;
};
// Types used for LUT
// ap_fixed<X,Y>
//https://web.csl.cornell.edu/courses/ece5775/pdf/lecture04.pdf
//Quantization behavior:
//Will increase the latency slightly
//Default mode: AP_TRN (truncation)
//Other rounding modes: AP_RND, AP_RND_ZERO, AP_RND_INF, ...
//Overflow behavior:
//AP_WRAP, AP_SAT
typedef ap_uint<11> slidingsum_t;
// Type used for LUT output
typedef ap_ufixed<15,1,AP_RND,AP_SAT> inverse_t;
//typedef ap_ufixed<15,1> inverse_t;
// Type used to store the result of the division
#define AP_FIXED_SIZE 20 //17
#define AP_FIXED_DEC 10 //7
typedef ap_fixed<AP_FIXED_SIZE,AP_FIXED_DEC,AP_RND,AP_SAT> zsliding_t;
//typedef ap_fixed<AP_FIXED_SIZE,AP_FIXED_DEC> zsliding_t;
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


///////////////////////////
// -- Inversion FUNCTION
// size of the LUT
#define N_TABLE_SIZE 1533 //Maximum number is 2045 for some reason (SIGSEGV otherwise)
template<class data_T, int N_TABLE>
void init_inversion_table(data_T table_out[N_TABLE]) {
    // Inversion function:
    //   result = 1/x
    for (int ii = 0; ii < N_TABLE; ii++) {
        // First, convert from table index to X-value (unsigned 8-bit, range 0 to +1533)
        float in_val = 1533.0*ii/float(N_TABLE);
        // Next, compute lookup table function
        table_out[ii] = (in_val>0) ? 1.0/in_val : 0.0;
        if (DEBUG>=2) std::cout << "table_out[" << ii << "] = " << table_out[ii] << std::endl;
    }
    return;
}

template<class data_T, class res_T, int TABLE_SIZE>
void inversion(data_T &data_den, res_T &res) {
    // Initialize the lookup table
    res_T inversion_table[TABLE_SIZE];
    init_inversion_table<res_T, TABLE_SIZE>(inversion_table);

    // Index into the lookup table based on data
    int index;

    //#pragma HLS PIPELINE

    if (data_den < 0) data_den = 0;
    if (data_den > TABLE_SIZE-1) data_den = TABLE_SIZE-1;
    index = data_den;
    res = inversion_table[index];
    if (DEBUG>=2) std::cout << "res = " << res << std::endl;
    return;
}

// Default table size provided here:
template<class data_T, class res_T>
void inversion(data_T &data_den, res_T &res) { 
    /* Get the inversion value from the LUT */
    if(data_den==0) {
        std::cout << "WARNING::inversion::data_den==0" << std::endl;
        return;
    }
    inversion<data_T, res_T, N_TABLE_SIZE>(data_den, res); 
    return;
}


///////////////////////////
// -- TRACK_QUALITY_CHECK_REF FUNCTION
//    Returns true if a track passes the quality requirements
//     and false otherwise
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


///////////////////////////
// -- FETCH_BIN_REF FUNCTION
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


///////////////////////////
// -- BIN_CENTER_REF FUNCTIONS
//    Input: Unsigned integer bin index [0,71]
//    Output: Signed integer bin center [-284,284] 
inline z0_t bin_center_ref(zbin_t iz) {
    if (DEBUG) std::cout<<"bin_center_ref::iz = " << iz << std::endl;
    if (DEBUG) std::cout<<"bin_center_ref::BHV_NHALFBINS = " << BHV_NHALFBINS << std::endl;
    int z = int(iz) - BHV_NHALFBINS;
    if (DEBUG) std::cout<<"bin_center_ref::(z << BNV_SHIFT) = " << (z << BNV_SHIFT) << std::endl;
    if (DEBUG) std::cout<<"bin_center_ref::( 1 << (BNV_SHIFT-1) ) = " << ( 1 << (BNV_SHIFT-1) ) << std::endl;
    return z0_t((z << BNV_SHIFT) + ( 1 << (BNV_SHIFT-1) ));
}
//    Input: Unsigned fixed [0,71]
//    Output: Signed integer bin center [-284,284] 
inline z0_t bin_center_ref(zsliding_t iz) {
    if (DEBUG) std::cout<<"bin_center_ref::iz = " << iz << std::endl;
    if (DEBUG) std::cout<<"bin_center_ref::BHV_NHALFBINS = " << BHV_NHALFBINS << std::endl;
    zsliding_t z = iz - BHV_NHALFBINS;
    if (DEBUG) std::cout<<"bin_center_ref::(z << BNV_SHIFT) = " << (z << BNV_SHIFT) << std::endl;
    if (DEBUG) std::cout<<"bin_center_ref::( 1 << (BNV_SHIFT-1) ) = " << ( 1 << (BNV_SHIFT-1) ) << std::endl;
    if (DEBUG) std::cout<<"bin_center_ref::(z << BNV_SHIFT)+( 1 << (BNV_SHIFT-1) ) = " <<(z << BNV_SHIFT) + ( 1 << (BNV_SHIFT-1) ) << std::endl;
    if (DEBUG) std::cout<<"bin_center_ref::z0_t((z << BNV_SHIFT)+( 1 << (BNV_SHIFT-1) )) = " <<z0_t((z << BNV_SHIFT) + ( 1 << (BNV_SHIFT-1) )) << std::endl;
    return z0_t((z << BNV_SHIFT) + ( 1 << (BNV_SHIFT-1) ));
}

///////////////////////////
// -- WINDOW_SUMPT FUNCTION
//    Sums the pT within a window given bin b 
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


///////////////////////////
// -- WEIGHTED_POSITION_REF FUNCTION
//    Returns the position of the vertex from [-288,288] 
inline z0_t weighted_position_ref(zbin_t iz, ptsum_t binpt[BHV_WINDOWSIZE], slidingsum_t sigma) {
    #pragma HLS array_partition variable=binpt complete dim=1
    zsliding_t zvtx_sliding = 0;
    slidingsum_t zvtx_sliding_sum = 0;
    if (DEBUG) std::cout << "zvtx_sliding (start) = " << zvtx_sliding_sum << std::endl;
    if (DEBUG) std::cout << "weighted_position_ref::iz = " << iz << std::endl;
    if (DEBUG) std::cout << "weighted_position_ref::zvtx_sliding_sum (weighted sum) = " << std::flush;
    WINDOWLOOP: for(unsigned int w = 0; w < BHV_WINDOWSIZE; ++w) {
        if (DEBUG) std::cout << "(" << binpt[w] << " * " << w << ")" << std::flush;
        if (DEBUG && w<BHV_WINDOWSIZE-1) std::cout<< " + " << std::flush;
        zvtx_sliding_sum += (binpt[w]*w);
    }
    if (DEBUG) std::cout << " = " << zvtx_sliding_sum << std::endl;
    if (DEBUG) std::cout << "weighted_position_ref::sigma = " << sigma << std::endl;
    if(sigma!=0) {
        inverse_t inv;
        //vivado_hls is smart enough to inline this function right here
        inversion<slidingsum_t, inverse_t>(sigma, inv);
        zvtx_sliding = zvtx_sliding_sum*inv;
    }
    else {
        zvtx_sliding = (BHV_WINDOWSIZE >> 1) + ((BHV_WINDOWSIZE%2!=0) ? 0.5 : 0.0);//bin_center_ref(iz+(BHV_WINDOWSIZE >> 1));
    }
    if (DEBUG) std::cout << "weighted_position_ref::zvtx_sliding (final) = " << zvtx_sliding << std::endl;
    zvtx_sliding += iz;
    //vivado_hls is smart enough to inline this function right here
    zvtx_sliding = bin_center_ref(zvtx_sliding);
    return zvtx_sliding;
}


///////////////////////////
// -- BIN_PLUS_HALF_WINDOW_REF FUNCTION
//    Returns the bin based position of the vertex as the
//     starting bin of the window plus half of the window size.
//    Rounding occurs when determining half of the window size.
//    The return value is an unsigned int.
template<class bin_t, class ret_t>
inline void bin_plus_half_window_ref(bin_t bin, ret_t &ret) {
    ret=bin+(BHV_WINDOWSIZE >> 1);
}

///////////////////////////
// -- FORWARD DECLARED FUNCTIONS
void bhv_add_track(zbin_vt zbin, pt_t tkpt, ptsum_t hist[BHV_NBINS]) ;
zbin_t bhv_find_pv(ptsum_t hist[BHV_NBINS], pt_t *sumpt, z0_t *pv_hw);
void bhv_find_pv_ref(TkObjExtended tracks[BHV_NSECTORS][BHV_NTRACKS], zbin_t & pvbin, z0_t & pv, int &pvsum) ;
//bool dummy(z0_t z0) ;
