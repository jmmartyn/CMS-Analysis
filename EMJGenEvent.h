#ifndef EmergingJetGenAnalysis_EMJGenAnalyzer_EMJGenEvent_h
#define EmergingJetGenAnalysis_EMJGenAnalyzer_EMJGenEvent_h

// Class describing the content of output for EmergingJetAnalyzer

#include <vector>
#include <functional>

#include "TLorentzVector.h"

#include "EmergingJetGenAnalysis/EMJGenAnalyzer/interface/OutputTree.h"
#define DEFAULTVALUE -1

using std::vector;

namespace EMJGen
{
  class Jet;
  class GenPart;
  class Event {
  public:
    Event() {}
    ~Event(){}
    void Init() {
      run                  = DEFAULTVALUE;
      lumi                 = DEFAULTVALUE;
      event                = DEFAULTVALUE;

      genpart_vector.clear();
      jet_vector.clear();
    };
    int    run                 ;
    int    lumi                ;
    int    event               ;

    vector<GenPart> genpart_vector;
    vector<Jet> jet_vector;
  }; 

  class GenPart {
  public:
    GenPart(){}
    ~GenPart(){}
    void Init(){
      index                = DEFAULTVALUE;
      pt                   = DEFAULTVALUE;
      mass                   = DEFAULTVALUE;
      eta                   = DEFAULTVALUE;
      phi                   = DEFAULTVALUE;
      pid                   = DEFAULTVALUE;
      ndau                   = DEFAULTVALUE;
      ndauch                   = DEFAULTVALUE;
      xdecay                = DEFAULTVALUE;
      ydecay                = DEFAULTVALUE;
      zdecay                = DEFAULTVALUE;

      KKPhoton_pt            = DEFAULTVALUE;
      KKPhoton_mass            = DEFAULTVALUE;
      KKPhoton_eta            = DEFAULTVALUE;
      KKPhoton_phi            = DEFAULTVALUE;

      Photon_pt            = DEFAULTVALUE;

      Radion_pt            = DEFAULTVALUE;
      Radion_mass            = DEFAULTVALUE;
      Radion_eta            = DEFAULTVALUE;
      Radion_phi            = DEFAULTVALUE;
      
    }

    void InitKKPhoton(){
      index                = DEFAULTVALUE;
      pt                   = DEFAULTVALUE;
      mass                   = DEFAULTVALUE;
      eta                   = DEFAULTVALUE;
      phi                   = DEFAULTVALUE;
      pid                   = DEFAULTVALUE;
      ndau                   = DEFAULTVALUE;
      ndauch                   = DEFAULTVALUE;
      xdecay                = DEFAULTVALUE;
      ydecay                = DEFAULTVALUE;
      zdecay                = DEFAULTVALUE;

      KKPhoton_pt            = DEFAULTVALUE;
      KKPhoton_mass            = DEFAULTVALUE;
      KKPhoton_eta            = DEFAULTVALUE;
      KKPhoton_phi            = DEFAULTVALUE;
      
    }

    void InitRadion(){
      index                = DEFAULTVALUE;
      pt                   = DEFAULTVALUE;
      mass                   = DEFAULTVALUE;
      eta                   = DEFAULTVALUE;
      phi                   = DEFAULTVALUE;
      pid                   = DEFAULTVALUE;
      ndau                   = DEFAULTVALUE;
      ndauch                   = DEFAULTVALUE;
      xdecay                = DEFAULTVALUE;
      ydecay                = DEFAULTVALUE;
      zdecay                = DEFAULTVALUE;

      Radion_pt            = DEFAULTVALUE;
      Radion_mass            = DEFAULTVALUE;
      Radion_eta            = DEFAULTVALUE;
      Radion_phi            = DEFAULTVALUE;
      
    }

    void InitPhoton(){
      index                = DEFAULTVALUE;
      pt                   = DEFAULTVALUE;
      mass                   = DEFAULTVALUE;
      eta                   = DEFAULTVALUE;
      phi                   = DEFAULTVALUE;
      pid                   = DEFAULTVALUE;
      ndau                   = DEFAULTVALUE;
      ndauch                   = DEFAULTVALUE;
      xdecay                = DEFAULTVALUE;
      ydecay                = DEFAULTVALUE;
      zdecay                = DEFAULTVALUE;

      Photon_pt            = DEFAULTVALUE;
      Photon_mass            = DEFAULTVALUE;
      Photon_eta            = DEFAULTVALUE;
      Photon_phi            = DEFAULTVALUE;
      
    }

 void InitGluon(){
      index                = DEFAULTVALUE;
      pt                   = DEFAULTVALUE;
      mass                   = DEFAULTVALUE;
      eta                   = DEFAULTVALUE;
      phi                   = DEFAULTVALUE;
      pid                   = DEFAULTVALUE;
      ndau                   = DEFAULTVALUE;
      ndauch                   = DEFAULTVALUE;
      xdecay                = DEFAULTVALUE;
      ydecay                = DEFAULTVALUE;
      zdecay                = DEFAULTVALUE;

      Gluon_pt            = DEFAULTVALUE;
      Gluon_mass            = DEFAULTVALUE;
      Gluon_eta            = DEFAULTVALUE;
      Gluon_phi            = DEFAULTVALUE;

    }

 void InitMET(){
      index                = DEFAULTVALUE;
      pt                   = DEFAULTVALUE;
      mass                   = DEFAULTVALUE;
      eta                   = DEFAULTVALUE;
      phi                   = DEFAULTVALUE;
      pid                   = DEFAULTVALUE;
      ndau                   = DEFAULTVALUE;
      ndauch                   = DEFAULTVALUE;
      xdecay                = DEFAULTVALUE;
      ydecay                = DEFAULTVALUE;
      zdecay                = DEFAULTVALUE;

      MET            = DEFAULTVALUE;
     
    }


    int    index               ;
    float  pt                  ;
    float  mass                  ;
    float  eta                  ;
    float  phi                  ;
    int  pid                  ;
    int  ndau                  ;
    int  ndauch                  ;
    float  xdecay                  ;
    float  ydecay                  ;
    float  zdecay                  ;

    float  KKPhoton_pt                  ;
    float  KKPhoton_mass                  ;
    float  KKPhoton_eta                  ;
    float  KKPhoton_phi                  ;

    float  Radion_pt                  ;
    float  Radion_mass                  ;
    float  Radion_eta                  ;
    float  Radion_phi                  ;

    float  Photon_pt                  ;
    float  Photon_mass                  ;
    float  Photon_eta                  ;
    float  Photon_phi                  ;

    float  Gluon_pt                  ;
    float  Gluon_mass                  ;
    float  Gluon_eta                  ;
    float  Gluon_phi                  ;

    float  MET                  ;

  };


  class Jet {
  public:
    Jet(){}
    ~Jet(){}
    void Init(){
      index                = DEFAULTVALUE;
      pt                   = DEFAULTVALUE;
      eta                   = DEFAULTVALUE;
      phi                   = DEFAULTVALUE;
    }
    int    index               ;
    float  pt                  ;
    float  eta                  ;
    float  phi                  ;
  };

// Turn vector of objects, into vector of member variable by calling func(object)
// Output is placed in provided vector reference
template <typename Object, typename T>
void
vectorize(const vector<Object>& input, std::function<T (const Object &)> func, vector<T>& output)
{
  output.clear();
  output.reserve(input.size()); // Doesn't reduce capacity
  for (const auto& obj : input) {
    output.push_back( func(obj) );
  }
}

// Version of vectorize() that returns a vector object
template <typename Object, typename T>
vector<T>
vectorize_new(const vector<Object>& input, std::function<T (const Object &)> func)
{
  vector<T> output;
  vectorize(input, func, output);
  return output;
}

using EMJGen::Event;
using EMJGen::GenPart;
using EMJGen::Jet;

void
WriteEventToOutput(const Event& event, EMJGen::OutputTree* otree)
{
  otree->Init(); // Reset all values and clear all vectors
  // Event-level variables, e.g. int, float, etc.
  {
    otree->run                  = event.run                 ;
    otree->lumi                 = event.lumi                ;
    otree->event                = event.event               ;
  }
  // generator-level variables
  {
    vectorize<GenPart, int   >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.index               ;}, otree->genpart_index               );
    vectorize<GenPart, float >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.pt                  ;}, otree->genpart_pt                  );
    vectorize<GenPart, float >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.mass                  ;}, otree->genpart_mass                  );
    vectorize<GenPart, float >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.eta                  ;}, otree->genpart_eta                  );
    vectorize<GenPart, float >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.phi                  ;}, otree->genpart_phi                  );
    vectorize<GenPart, int >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.pid                  ;}, otree->genpart_pid                  );
    vectorize<GenPart, int >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.ndau                  ;}, otree->genpart_ndau                  );
    vectorize<GenPart, int >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.ndauch                  ;}, otree->genpart_ndauch                  );
    vectorize<GenPart, float >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.xdecay                  ;}, otree->genpart_xdecay                  );
    vectorize<GenPart, float >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.ydecay                  ;}, otree->genpart_ydecay                  );
    vectorize<GenPart, float >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.zdecay                  ;}, otree->genpart_zdecay                  );

    vectorize<GenPart, float >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.KKPhoton_pt                  ;}, otree->genpart_KKPhoton_pt                  );    
    vectorize<GenPart, float >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.KKPhoton_mass                  ;}, otree->genpart_KKPhoton_mass                  );    
    vectorize<GenPart, float >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.KKPhoton_eta                  ;}, otree->genpart_KKPhoton_eta                  );    
    vectorize<GenPart, float >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.KKPhoton_phi                  ;}, otree->genpart_KKPhoton_phi                  ); 
   
    vectorize<GenPart, float >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.Photon_pt                  ;}, otree->genpart_Photon_pt                  );    
    vectorize<GenPart, float >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.Photon_mass                  ;}, otree->genpart_Photon_mass                  );    
    vectorize<GenPart, float >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.Photon_eta                  ;}, otree->genpart_Photon_eta                  );    
    vectorize<GenPart, float >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.Photon_phi                  ;}, otree->genpart_Photon_phi                  ); 
   
    vectorize<GenPart, float >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.Radion_pt                  ;}, otree->genpart_Radion_pt                  );    
    vectorize<GenPart, float >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.Radion_mass                  ;}, otree->genpart_Radion_mass                  );    
    vectorize<GenPart, float >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.Radion_eta                  ;}, otree->genpart_Radion_eta                  );    
    vectorize<GenPart, float >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.Radion_phi                  ;}, otree->genpart_Radion_phi                  );
     
    vectorize<GenPart, float >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.Gluon_pt                  ;}, otree->genpart_Gluon_pt                  );    
    vectorize<GenPart, float >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.Gluon_mass                  ;}, otree->genpart_Gluon_mass                  );    
    vectorize<GenPart, float >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.Gluon_eta                  ;}, otree->genpart_Gluon_eta                  );    
    vectorize<GenPart, float >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.Gluon_phi                  ;}, otree->genpart_Gluon_phi                  ); 
   
    vectorize<GenPart, float >(event.genpart_vector, [](const EMJGen::GenPart& obj ){return obj.MET                  ;}, otree->genpart_MET                  );    
}

  // jet variables
  {
    vectorize<Jet, int   >(event.jet_vector, [](const EMJGen::Jet& obj ){return obj.index               ;}, otree->genjet_index               );
    vectorize<Jet, float >(event.jet_vector, [](const EMJGen::Jet& obj ){return obj.pt                  ;}, otree->genjet_pt                  );
    vectorize<Jet, float >(event.jet_vector, [](const EMJGen::Jet& obj ){return obj.eta                  ;}, otree->genjet_eta                  );
    vectorize<Jet, float >(event.jet_vector, [](const EMJGen::Jet& obj ){return obj.phi                  ;}, otree->genjet_phi                  );
  }

}

}
#endif
