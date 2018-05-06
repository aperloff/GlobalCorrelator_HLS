#include <cstdio>
#include "firmware/bram_sliding_window_vtx.h"
#include "../utils/DiscretePFInputs_IO.h"
#include "../utils/pattern_serializer.h"
#define NTEST 250

int main() {
    DiscretePFInputs inputs("barrel_alltracks_sectors_1x12_TTbar_PU140.dump");
    //DiscretePFInputs inputs("barrel_sectors_1x12_TTbar_PU140.dump");

    TkObjExtended track_in[BHV_NSECTORS][BHV_NTRACKS], track_tmp[2*BHV_NTRACKS];
    //TkObj track_in[BHV_NSECTORS][BHV_NTRACKS], track_tmp[2*BHV_NTRACKS];

    unsigned int ngood = 0, ntot = 0, ncmssw_good = 0, nagree = 0, nerr = 0;
    double resol = 0;
    double cmssw_resol = 0;
    // run multiple tests

    FILE *test_in  = fopen("bhv_dump_in.txt","w");
    FILE *test_out = fopen("bhv_dump_out.txt","w");
    unsigned int test_in_frame = 0, test_out_frame = 0;

    for (int test = 1; test <= NTEST; ++test) {
        // read the event
        if (!inputs.nextEvent()) break;
        if (inputs.event().regions.size() != N_IN_SECTORS) { printf("ERROR: Mismatching number of input regions: %lu\n", inputs.event().regions.size()); return 2; }
        for (int is = 0; is < N_IN_SECTORS; ++is) {
            const Region & r = inputs.event().regions[is];
            dpf2fw::convert<2*BHV_NTRACKS>(r.track, track_tmp);
            for (int i = 0, t = 0; i < BHV_NTRACKS; ++i) {
                track_in[2*is+0][i] = track_tmp[t++];
                track_in[2*is+1][i] = track_tmp[t++];
            }
        }

        z0_t   pv_gen = round(inputs.event().genZ0*l1tpf_int::InputTrack::Z0_SCALE);
        z0_t   pv_cmssw = round(inputs.event().z0*l1tpf_int::InputTrack::Z0_SCALE);
        zbin_t pvbin_ref = 0;
        z0_t   pv_ref = 0;
        int    ptsum_ref = 0;

        REF_FUNC(track_in, pvbin_ref, pv_ref, ptsum_ref);


        //
        // Make MP7 input file
        //
        for (unsigned int it = 0; it < 18; ++it) {
            for (int go = 1; go >= 0; --go) {
                fprintf(test_in,"Frame %4d : %1d  %3d  %6d\n", test_in_frame++, go, int(fetch_bin_ref(track_in[0][it].hwZ0).bin), int(track_in[0][it].hwPt>>1));
            }
        }

        //
        // Make MP7 output file
        //
        ptsum_t hist[BHV_NBINS];
        for (unsigned int b = 0; b < BHV_NBINS; ++b) hist[b] = 0;
        for (int is = 0; is < BHV_NSECTORS; ++is) {
            for (unsigned int it = 0; it < BHV_NTRACKS; ++it) {
                if ((QUALITY && track_quality_check_ref(track_in[is][it])) || (!QUALITY)) {
                    bhv_add_track(fetch_bin_ref(track_in[is][it].hwZ0), track_in[is][it].hwPt, hist);
                }
                if (it == 17 && is == 0) {
                    ptsum_t max = 0; zbin_t bmax = BHV_NBINS;
                    for (unsigned int b = 0; b < BHV_NBINS; b += 2) {
                       fprintf(test_out,"Frame %4d : %3d   %6d   %6d\n", test_out_frame++, int(b), int(hist[b]), int(hist[b+1]));
                    }
                }
            }
        }

        pt_t ptsum_hw;
        z0_t pv_hw;
        zbin_t pvbin_hw;
        //TOP_FUNC(track_in, &pvbin_hw, &pv_hw, &ptsum_hw);
        TOP_FUNC(hist, &pvbin_hw, &pv_hw, &ptsum_hw);
        
        #ifdef TESTBOARD
            if (!VALIDATE) continue;
        #endif

        bool match_gen = abs(int(pv_ref-pv_gen)) <= 10;
        if (abs(int(pv_ref-pv_gen)) <= 10) { ngood++; resol += std::pow(double(pv_gen - pv_ref), 2); }
        ntot++;
        if (abs(int(pv_gen-pv_cmssw)) <= 10) { ncmssw_good++; cmssw_resol += std::pow(double(pv_gen - pv_cmssw), 2); }

        if (abs(int(pv_ref-pv_cmssw)) <= 10) nagree++;
        printf("%sGEN PV %+4d    CMSSW PV %+4d  bin %3d :  REF %+4d  bin %3d, ptsum %8d, diff %+4d :  HW %+4d  bin %3d, ptsum %8d diff %+4d  :  MATCH %+1d%s\n",
               (!match_gen)?"\e[1;31m":"", int(pv_gen), int(pv_cmssw), int(fetch_bin_ref(pv_cmssw).bin), int(pv_ref), int(pvbin_ref), ptsum_ref, int(pv_ref-pv_gen),
               int(pv_hw), int(pvbin_hw), int(ptsum_hw), int(pv_hw-pv_gen), int(pv_hw-pv_ref),(!match_gen)?"\e[0m":"");
        if (pv_ref != pv_hw || ptsum_hw != ptsum_ref) { nerr++; printf("ERROR!!!\n"); }
    }

    fclose(test_in);
    fclose(test_out);
    printf("Good matches: CMSSW %4d/%4d = %.3f  REF %4d/%4d = %.3f  (REF vs CMSSW: %4d/%4d = %.3f). ERRORS: %d\n", 
        ncmssw_good, ntot, float(ncmssw_good)/(ntot), ngood, ntot, float(ngood)/(ntot), nagree, ntot, float(nagree)/(ntot), nerr);
    printf("Resolution: CMSSW: %5.3f mm   REF %5.3f mm\n", 0.5*sqrt(cmssw_resol/(ncmssw_good)), 0.5*sqrt(resol/(ngood)));
    return 0;
}
