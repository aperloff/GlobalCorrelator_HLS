#include "../../firmware/data.h"
#include "../../regionizer/firmware/regionizer.h"
#if (defined(__GXX_EXPERIMENTAL_CXX0X__) and defined(CMSSW))
  #define DEBUG 1
#else
  #define DEBUG 0
#endif
#if DEBUG > 0
    #include <iostream>
#endif

//////////////
// -- TYPEDEFS
typedef ap_uint<11> ptsum_t;
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

///////////////////////////
// -- PRECOMPILER FUNCTIONS
//From: https://stackoverflow.com/questions/27581671/how-to-compute-log-with-the-preprocessor
#define IS_REPRESENTIBLE_IN_D_BITS(D, N)                \
   (((unsigned long) N >= (1UL << (D - 1)) && (unsigned long) N < (1UL << D)) ? D : -1)
#define BITS_TO_REPRESENT(N)                             \
   (N == 0 ? 1 : (31                                     \
                  + IS_REPRESENTIBLE_IN_D_BITS( 1, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS( 2, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS( 3, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS( 4, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS( 5, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS( 6, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS( 7, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS( 8, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS( 9, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS(10, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS(11, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS(12, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS(13, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS(14, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS(15, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS(16, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS(17, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS(18, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS(19, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS(20, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS(21, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS(22, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS(23, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS(24, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS(25, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS(26, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS(27, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS(28, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS(29, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS(30, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS(31, N)    \
                  + IS_REPRESENTIBLE_IN_D_BITS(32, N)    \
                  )                                      \
      )
#define PWRTWO(x) (1 << (x))

///////////////////////
// -- DEFINE STATEMENTS
#define BHV_MAXPT 100
//#define BHV_MAXBIN 511
#define BHV_MAXBIN 511
#define BHV_CHI2MAX 100
#define BHV_PTMINTRA 2
#define BHV_NSTUBSMIN 4
//0 : truncation. Tracks with PT above PTMAX are ignored
//1 : saturation . Tracks with PT above PTMAX are set to PT=PTMAX.
#define BHV_TRUNCSAT 1

#define BNV_SHIFT 3
#define BHV_NBINS 72
#define BHV_NHALFBINS (BHV_NBINS/2)
#define BHV_NSECTORS 2*N_IN_SECTORS
#define BHV_NTRACKS 18
#define BHV_WINDOWSIZE 1
#define BHV_NSUMS BHV_NBINS-BHV_WINDOWSIZE+1
#define BITS_TO_STORE_NSUMS BITS_TO_REPRESENT(BHV_NSUMS)+1
#define MAXIMUM_SEARCH_SIZE PWRTWO(BITS_TO_STORE_NSUMS-1)
#define MAXIMUM_SEARCH_SIZE_ITERATIONS BITS_TO_REPRESENT(MAXIMUM_SEARCH_SIZE)-1

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
        #if (DEBUG==2)
            std::cout << "table_out[" << ii << "] = " << table_out[ii] << std::endl;
        #endif
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
    #if (DEBUG==2)
        std::cout << "res = " << res << std::endl;
    #endif
    return;
}

// Default table size provided here:
template<class data_T, class res_T>
void inversion(data_T &data_den, res_T &res) { 
    /* Get the inversion value from the LUT */
    if(data_den==0) {
        printf("WARNING::inversion::data_den==0\n");
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
    int z = int(iz) - BHV_NHALFBINS;
    #if (DEBUG==1)
        std::cout<<"bin_center_ref::iz = " << iz << std::endl;
        std::cout<<"bin_center_ref::BHV_NHALFBINS = " << BHV_NHALFBINS << std::endl;
        std::cout<<"bin_center_ref::(z << BNV_SHIFT) = " << (z << BNV_SHIFT) << std::endl;
        std::cout<<"bin_center_ref::( 1 << (BNV_SHIFT-1) ) = " << ( 1 << (BNV_SHIFT-1) ) << std::endl;
    #endif
    return z0_t((z << BNV_SHIFT) + ( 1 << (BNV_SHIFT-1) ));
}
//    Input: Unsigned fixed [0,71]
//    Output: Signed integer bin center [-284,284] 
inline z0_t bin_center_ref(zsliding_t iz) {
    zsliding_t z = iz - BHV_NHALFBINS;
    #if (DEBUG==1)
        std::cout<<"bin_center_ref::iz = " << iz << std::endl;
        std::cout<<"bin_center_ref::BHV_NHALFBINS = " << BHV_NHALFBINS << std::endl;
        std::cout<<"bin_center_ref::(z << BNV_SHIFT) = " << (z << BNV_SHIFT) << std::endl;
        std::cout<<"bin_center_ref::( 1 << (BNV_SHIFT-1) ) = " << ( 1 << (BNV_SHIFT-1) ) << std::endl;
        std::cout<<"bin_center_ref::(z << BNV_SHIFT)+( 1 << (BNV_SHIFT-1) ) = " <<(z << BNV_SHIFT) + ( 1 << (BNV_SHIFT-1) ) << std::endl;
        std::cout<<"bin_center_ref::z0_t((z << BNV_SHIFT)+( 1 << (BNV_SHIFT-1) )) = " <<z0_t((z << BNV_SHIFT) + ( 1 << (BNV_SHIFT-1) )) << std::endl;
    #endif
    return z0_t((z << BNV_SHIFT) + ( 1 << (BNV_SHIFT-1) ));
}

///////////////////////////
// -- WINDOW_SUMPT FUNCTION
//    Sums the pT within a window given bin b
/*
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
*/

///////////////////////////
// -- WEIGHTED_POSITION_REF FUNCTION
//    Returns the position of the vertex from [-288,288] 
inline z0_t weighted_position_ref(zbin_t iz, ptsum_t binpt[BHV_WINDOWSIZE], slidingsum_t sigma) {
    #pragma HLS array_partition variable=binpt complete dim=1
    zsliding_t zvtx_sliding = 0;
    slidingsum_t zvtx_sliding_sum = 0;
    #if (DEBUG==1) 
        std::cout << "zvtx_sliding (start) = " << zvtx_sliding_sum << std::endl;
        std::cout << "weighted_position_ref::iz = " << iz << std::endl;
        std::cout << "weighted_position_ref::zvtx_sliding_sum (weighted sum) = " << std::flush;
    #endif
    WINDOWLOOP: for(ap_uint<BITS_TO_REPRESENT(BHV_NBINS)> w = 0; w < BHV_WINDOWSIZE; ++w) {
        #if (DEBUG==1)
            std::cout << "(" << binpt[w] << " * " << w << ")" << std::flush;
            if (w<BHV_WINDOWSIZE-1) std::cout<< " + " << std::flush;
        #endif
        //Might be a problem in future due to timing
        //Check if the RTL doesn't meet the timing requirement
        zvtx_sliding_sum += (binpt[w]*w);
    }
    #if (DEBUG==1)
        std::cout << " = " << zvtx_sliding_sum << std::endl;
        std::cout << "weighted_position_ref::sigma = " << sigma << std::endl;
    #endif
    if(sigma!=0) {
        inverse_t inv;
        //vivado_hls is smart enough to inline this function right here
        inversion<slidingsum_t, inverse_t>(sigma, inv);
        zvtx_sliding = zvtx_sliding_sum*inv;
    }
    else {
        zvtx_sliding = (BHV_WINDOWSIZE >> 1) + ((BHV_WINDOWSIZE%2!=0) ? 0.5 : 0.0);//bin_center_ref(iz+(BHV_WINDOWSIZE >> 1));
    }
    #if (DEBUG==1) 
        std::cout << "weighted_position_ref::zvtx_sliding (final) = " << zvtx_sliding << std::endl;
    #endif
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
// -- BHV_ACCESS_ARRAY FUNCTION
//    Returns the bin value in an array corresponding to the
//     index created by summing i and j;
template<class bin_t, class ret_t, int ARRAY_SIZE>
inline ret_t bhv_access_array(ret_t array[], bin_t i, bin_t j) {
    return array[i+j];
}

///////////////////////////////
// -- BIN_COMPUTE_SUMS FUNCTION
//    Computes the sum of the bins within a window and places those
//     sums in a second array.
//    The array of computed sums is returned by reference.
inline void bhv_compute_sums(ptsum_t hist[BHV_NBINS], slidingsum_t hist_window_sums[BHV_NSUMS]) {
    // Loop through all bins, taking into account the fact that the last bin is
    //  nbins-window_width+1
    BINLOOP: for (unsigned int b = 0; b < BHV_NSUMS; ++b) {
        SUMLOOP: for (unsigned int w = 0; w < BHV_WINDOWSIZE;  ++w) {
            hist_window_sums[b] += hist[b+w];
        }
    }
}

///////////////////////////
// -- BIN_FIND_MAX FUNCTION
//    A simple linear search algorithm for finding the index
//     of an array with the maximum value.
//    Both the maximum value and the bin index are returned
//     by reference.
inline void bhv_find_max(slidingsum_t hist_window_sums[BHV_NSUMS], zbin_t &b_max, slidingsum_t &sigma_max) {
    SUMLOOP: for (unsigned int b = 0; b < BHV_NSUMS; ++b) {
        if(hist_window_sums[b] > sigma_max) {
            sigma_max = hist_window_sums[b];
            b_max = b;
        }
    }
}

////////////////////////////////////
// -- BIN_PARALLEL_FIND_MAX FUNCTION
//    A parallel search algorithm for finding the index
//     of an array with the maximum value.
//    Both the maximum value and the bin index are returned
//     by reference.
//    Also in this section are the associated helper functions
template<class array_t, int SIZE>
inline void initialize_array(array_t arr[]) {
  INITIALIZELOOP: for (unsigned int b = 0; b < SIZE; ++b) {
        #pragma HLS UNROLL
        arr[b] = 0;
    }
}

template<class array_t, class array2_t, int SIZE_SRC, int SIZE_DST>
inline void copy_to_p2_array(array_t input_array[SIZE_SRC], array_t output_array[SIZE_DST], array2_t output_array_index[SIZE_DST]) {
  COPYP2LOOP1: for(unsigned int i=0; i<SIZE_SRC; ++i) {
    #pragma HLS UNROLL
    output_array[i]       = input_array[i];
    output_array_index[i] = i;
  }
  COPYP2LOOP2: for(unsigned int i=SIZE_SRC; i<SIZE_DST; ++i) {
    #pragma HLS UNROLL
    output_array[i]       = 0;
    output_array_index[i] = MAXIMUM_SEARCH_SIZE;
  }
}

template<class array_t, int SIZE>
inline void copy_array(array_t input_array[SIZE], array_t output_array[SIZE]) {
  COPYLOOP: for(unsigned int i=0; i<SIZE; ++i) {
    #pragma HLS UNROLL
    output_array[i] = input_array[i];
  }
}

template<class type1_t, class type2_t>
inline void comparator(type1_t bin1, type1_t bin2, type2_t binindex1, type2_t binindex2, type1_t &res, type2_t &resindex) {
  if (bin1 >= bin2) {
    res = bin1;
    resindex = binindex1; 
  }
  else {
    res = bin2;
    resindex = binindex2;
  }
}

template<class value_t, class index_t, int SIZE>
inline void bhv_parallel_find_max(value_t input_array[SIZE], index_t &b_max, value_t &sigma_max) {
  value_t values_array[MAXIMUM_SEARCH_SIZE];
  index_t index_array[MAXIMUM_SEARCH_SIZE];
  value_t larger_of_pair_array[MAXIMUM_SEARCH_SIZE];
  index_t larger_of_pair_index_array[MAXIMUM_SEARCH_SIZE];
  #pragma HLS array_partition variable=input_array complete dim=1
  #pragma HLS array_partition variable=values_array complete dim=1
  #pragma HLS array_partition variable=index_array complete dim=1
  #pragma HLS array_partition variable=larger_of_pair_array complete dim=1
  #pragma HLS array_partition variable=larger_of_pair_index_array complete dim=1

  copy_to_p2_array<value_t,index_t,SIZE,MAXIMUM_SEARCH_SIZE>(input_array,values_array,index_array);

  ITERATIONLOOP: for (unsigned int iteration=0; iteration<MAXIMUM_SEARCH_SIZE_ITERATIONS; ++iteration) {
        #pragma HLS UNROLL
        initialize_array<value_t,MAXIMUM_SEARCH_SIZE>(larger_of_pair_array);
        initialize_array<index_t,MAXIMUM_SEARCH_SIZE>(larger_of_pair_index_array);
        COMPARATORLOOP: for (int pair=0; pair<PWRTWO(MAXIMUM_SEARCH_SIZE_ITERATIONS-iteration); pair+=2) {
            #pragma HLS UNROLL
            comparator<value_t,index_t>(values_array[pair], values_array[pair+1], index_array[pair], index_array[pair+1], larger_of_pair_array[pair>>1],larger_of_pair_index_array[pair>>1]);
        }
        copy_array<value_t,MAXIMUM_SEARCH_SIZE>(larger_of_pair_array,values_array);
        copy_array<index_t,MAXIMUM_SEARCH_SIZE>(larger_of_pair_index_array,index_array);
    }

    b_max     = index_array[0];
    sigma_max = values_array[0];
}

///////////////////////////
// -- FORWARD DECLARED FUNCTIONS
void bhv_add_track(zbin_vt zbin, pt_t tkpt, ptsum_t hist[BHV_NBINS]) ;
void bhv_find_pv(ptsum_t hist[BHV_NBINS], zbin_t *pvbin_hw, z0_t *pv_hw, pt_t *sumpt);
void bhv_find_pv_ref(TkObjExtended tracks[BHV_NSECTORS][BHV_NTRACKS], zbin_t & pvbin, z0_t & pv, int &pvsum) ;
//bool dummy(z0_t z0) ;
