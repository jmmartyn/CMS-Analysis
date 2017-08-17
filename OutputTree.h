#ifndef EmergingJetGenAnalysis_EMJGenAnalyzer_OutputTree_h
#define EmergingJetGenAnalysis_EMJGenAnalyzer_OutputTree_h

#include <vector>

#include "TTree.h"

using std::vector;

namespace EMJGen
{
  class OutputTree {
  public:
    OutputTree() { Init(); }
    void Init();

    void Branch(TTree* tree);

    int                     run                 ;
    int                     lumi                ;
    int                     event               ;
    vector<int>             genpartall_pid      ;//save pid's of all the partciles in the event
    
    vector<int>             genpart_index               ;
    vector<float>           genpart_pt                  ;
    vector<float>           genpart_mass                  ;
    vector<float>           genpart_eta;
    vector<float>           genpart_phi;
    vector<int>           genpart_pid;
    vector<int>           genpart_ndau;
    vector<int>           genpart_ndauch;
    vector<float>           genpart_xdecay;
    vector<float>           genpart_ydecay;
    vector<float>           genpart_zdecay;

    vector<float>           genpart_KKPhoton_pt                  ;
    vector<float>           genpart_KKPhoton_mass                  ;
    vector<float>           genpart_KKPhoton_eta                  ;
    vector<float>           genpart_KKPhoton_phi                  ;

    vector<float>           genpart_Photon_pt                  ;
    vector<float>           genpart_Photon_mass                  ;
    vector<float>           genpart_Photon_eta                  ;
    vector<float>           genpart_Photon_phi                  ;

    vector<float>           genpart_Radion_pt                  ;
    vector<float>           genpart_Radion_mass                  ;
    vector<float>           genpart_Radion_eta                  ;
    vector<float>           genpart_Radion_phi                  ;

    vector<int>             genjet_index               ;
    vector<float>           genjet_pt                  ;
    vector<float>           genjet_eta;
    vector<float>           genjet_phi;

    vector<float>           genpart_Gluon_pt                  ;
    vector<float>           genpart_Gluon_mass                  ;
    vector<float>           genpart_Gluon_eta                  ;
    vector<float>           genpart_Gluon_phi                  ;

    vector<float>           genpart_MET                  ;
  };
}

void
EMJGen::OutputTree::Init() {
  run                 = -1;
  lumi                = -1;
  event               = -1;
  genpartall_pid              .clear();
  genpart_index               .clear();
  genpart_pt                  .clear();
  genpart_mass                  .clear();
  genpart_eta                  .clear();
  genpart_phi                  .clear();
  genpart_pid                  .clear();
  genpart_ndau                  .clear();
  genpart_ndauch                  .clear();
  genpart_xdecay                  .clear();
  genpart_ydecay                  .clear();
  genpart_zdecay                  .clear();

  genpart_KKPhoton_pt                  .clear();
  genpart_KKPhoton_mass                  .clear();
  genpart_KKPhoton_eta                  .clear();
  genpart_KKPhoton_phi                  .clear();

  genpart_Photon_pt                  .clear();
  genpart_Photon_mass                  .clear();
  genpart_Photon_eta                  .clear();
  genpart_Photon_phi                  .clear();

  genpart_Radion_pt                  .clear();
  genpart_Radion_mass                  .clear();
  genpart_Radion_eta                  .clear();
  genpart_Radion_phi                  .clear();

  genpart_Gluon_pt                  .clear();
  genpart_Gluon_mass                  .clear();
  genpart_Gluon_eta                  .clear();
  genpart_Gluon_phi                  .clear();

  genpart_MET                  .clear();

  genjet_index               .clear();
  genjet_pt                  .clear();
  genjet_eta                  .clear();
  genjet_phi                  .clear();
}

void
EMJGen::OutputTree::Branch(TTree* tree) {
#define BRANCH(tree, branch) (tree)->Branch(#branch, &branch);
  BRANCH(tree, run                 );
  BRANCH(tree, lumi                );
  BRANCH(tree, event               );

  BRANCH(tree, genpartall_pid               );
  BRANCH(tree, genpart_index               );
  BRANCH(tree, genpart_pt                  );
  BRANCH(tree, genpart_mass                  );
  BRANCH(tree, genpart_eta                  );
  BRANCH(tree, genpart_phi                  );
  BRANCH(tree, genpart_pid                  );
  BRANCH(tree, genpart_ndau                  );
  BRANCH(tree, genpart_ndauch                  );
  BRANCH(tree, genpart_xdecay                  );
  BRANCH(tree, genpart_ydecay                  );
  BRANCH(tree, genpart_zdecay                  );

  BRANCH(tree, genpart_KKPhoton_pt                  );
  BRANCH(tree, genpart_KKPhoton_mass                  );
  BRANCH(tree, genpart_KKPhoton_eta                  );
  BRANCH(tree, genpart_KKPhoton_phi                  );

  BRANCH(tree, genpart_Photon_pt                  );
  BRANCH(tree, genpart_Photon_mass                  );
  BRANCH(tree, genpart_Photon_eta                  );
  BRANCH(tree, genpart_Photon_phi                  );

  BRANCH(tree, genpart_Radion_pt                  );
  BRANCH(tree, genpart_Radion_mass                  );
  BRANCH(tree, genpart_Radion_eta                  );
  BRANCH(tree, genpart_Radion_phi                  );

  BRANCH(tree, genpart_Gluon_pt                  );
  BRANCH(tree, genpart_Gluon_mass                  );
  BRANCH(tree, genpart_Gluon_eta                  );
  BRANCH(tree, genpart_Gluon_phi                  );

  BRANCH(tree, genpart_MET                  );

  BRANCH(tree, genjet_index               );
  BRANCH(tree, genjet_pt                  );
  BRANCH(tree, genjet_eta                  );
  BRANCH(tree, genjet_phi                  );


}

// Insert new empty element in nested vector and returns pointer to the added element
template <typename T>
vector<T>&
make_new_element (vector< vector<T> > & vec) {
  vec.push_back( vector<T>() );
  return vec.back();
}


#endif
