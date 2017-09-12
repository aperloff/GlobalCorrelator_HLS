
#include "DiscretePFInputs.h"
#include "src/data.h"
#include <vector>

namespace dpf2fw {

    // convert inputs from discrete to firmware
    void convert(const l1tpf_int::PropagatedTrack & in, TkObj &out) {
        out.hwPt = in.hwPt;
        out.hwPtErr = in.hwCaloPtErr;
        out.hwEta = in.hwEta; // @calo
        out.hwPhi = in.hwPhi; // @calo
        out.hwZ0 = in.hwZ0;
    }
    void convert(const l1tpf_int::CaloCluster & in, HadCaloObj & out) {
        out.hwPt = in.hwPt;
        out.hwEmPt = in.hwEmPt;
        out.hwEta = in.hwEta;
        out.hwPhi = in.hwPhi;
        out.hwIsEM = in.isEM;
    }
    void convert(const l1tpf_int::CaloCluster & in, EmCaloObj & out) {
        out.hwPt = in.hwPt;
        out.hwPtErr = in.hwPtErr;
        out.hwEta = in.hwEta;
        out.hwPhi = in.hwPhi;
    }
    void convert(const l1tpf_int::Muon & in, MuObj & out) {
        out.hwPt = in.hwPt;
        out.hwPtErr = 0; // does not exist in input
        out.hwEta = in.hwEta; // @calo
        out.hwPhi = in.hwPhi; // @calo
    }

    template<unsigned int NMAX, typename In, typename Out>
    void convert(const std::vector<In> & in, Out out[NMAX]) {
        for (unsigned int i = 0, n = std::min<unsigned int>(NMAX, in.size()); i < n; ++i) {
            convert(in[i], out[i]);
        }
    }


};
