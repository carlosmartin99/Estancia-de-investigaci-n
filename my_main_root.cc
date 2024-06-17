

#include "Pythia8/Pythia.h" // Include Pythia headers.
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
using namespace Pythia8;

void FillMuonData(Event& event, int iNew, TTree* outputTree, TH2F* histogram)
{

    Particle& particle = event[iNew];

    // Check if the particle has daughters
    if (particle.daughterList().size() > 0) {

        float mu_plus_px, mu_plus_py, mu_plus_pz, mu_plus_energy; // empty variable to be filled a bit later
        float mu_minus_px, mu_minus_py, mu_minus_pz, mu_minus_energy;
        float mu_plus_pt, mu_minus_pt, mu_plus_p, mu_minus_p;

        outputTree->SetBranchAddress("mu_p_PX", &mu_plus_px); // path where to be filled
        outputTree->SetBranchAddress("mu_p_PY", &mu_plus_py);
        outputTree->SetBranchAddress("mu_p_PZ", &mu_plus_pz);
        outputTree->SetBranchAddress("mu_p_ENERGY", &mu_plus_energy);
        outputTree->SetBranchAddress("mu_m_PX", &mu_minus_px);
        outputTree->SetBranchAddress("mu_m_PY", &mu_minus_py);
        outputTree->SetBranchAddress("mu_m_PZ", &mu_minus_pz);
        outputTree->SetBranchAddress("mu_m_ENERGY", &mu_minus_energy);

        // Loop over daughters
        for (int j = 0; j < particle.daughterList().size(); ++j) {

            int daughterIndex = particle.daughterList()[j]; // get daughter particle index

            Particle& daughter = event[daughterIndex]; // get daughter particle

            // std::cout << "  Daughter " << j << ": ID = " << daughter.id() << ", Status = " << daughter.status() << std::endl;
            if (daughter.id() == 13){ // if the daughter is a muon minus
                mu_minus_px = daughter.px()*1e3; // in pythia, the units are GeV, so we convert them to MeV
                mu_minus_py = daughter.py()*1e3;
                mu_minus_pz = daughter.pz()*1e3;
                mu_minus_energy = daughter.e()*1e3;
                mu_minus_pt = sqrt(daughter.px()*daughter.px() + daughter.py()*daughter.py())*1e3;
                mu_minus_p = sqrt(daughter.px()*daughter.px() + daughter.py()*daughter.py() + daughter.pz()*daughter.pz())*1e3;
                histogram->Fill(mu_minus_p, mu_minus_pt); // fill the histogram with the data
            }
            else if (daughter.id() == -13){ // if the daughter is a muon plus do same
                mu_plus_px = daughter.px()*1e3;
                mu_plus_py = daughter.py()*1e3;
                mu_plus_pz = daughter.pz()*1e3;
                mu_plus_energy = daughter.e()*1e3;
                mu_plus_pt = sqrt(daughter.px()*daughter.px() + daughter.py()*daughter.py())*1e3;
                mu_plus_p = sqrt(daughter.px()*daughter.px() + daughter.py()*daughter.py() + daughter.pz()*daughter.pz())*1e3;
                histogram->Fill(mu_plus_p, mu_plus_pt);
            }
        }
        outputTree->Fill(); // filling output tree
    }
    // }
}

void fillParticle(auto& pythia, int id, double ee, double thetaIn, double phiIn,
  Event& event, ParticleData& pdt, Rndm& rndm,  TTree* dataTree, TTree* outputTree, TH2F* histogram, bool hasLifetime = false) {

    bool atRest = false;

    // Reset event record to allow for new event.
    event.reset(); // Reset the event record.

    // Select particle mass; where relevant according to Breit-Wigner.
    double mm = pdt.mSel(id);

    int numberOfEntries = dataTree->GetEntries(); // getting total number of events in ROOT FILE

    int entry = int(rndm.flat() * numberOfEntries); // get random event from the file
    dataTree->GetEntry(entry); // get the event

    float px, py, pz, higgs_energy; // empty variables to be filled a bit later
    dataTree->SetBranchAddress("Scalar_PX",&px); // read px value from input file
    dataTree->SetBranchAddress("Scalar_PY",&py);
    dataTree->SetBranchAddress("Scalar_PZ",&pz);

    cout << "Data from ROOT Tree for event " << entry << " : px: " << px/1e3 << " py: " << py/1e3 << " pz: " << pz/1e3 << endl;

    float px_gev = px/1e3; // convert MeV to GeV
    float py_gev = py/1e3; // convert MeV to GeV
    float pz_gev = pz/1e3;

    double calculated_energy = sqrt(px_gev*px_gev + py_gev*py_gev + pz_gev*pz_gev + mm*mm);

    int iNew = event.append( id, 1, 0, 0, px_gev, py_gev, pz_gev, calculated_energy, mm); // Add the particle to the event record

    if (hasLifetime) event[iNew].tau( event[iNew].tau0() * rndm.exp() ); // Set the lifetime of the particle

    pythia.next(); // process the event -> calculate all decays and kinematics

    FillMuonData(event, iNew, outputTree, histogram); // fill the output tree with muon data

}

int main() {
	Pythia pythia; // Declare Pythia object

    TFile *outputFile = TFile::Open("output.root", "RECREATE"); // create file, where we will write data from PYTHIA
    TTree muon_tree("MuonsTree", "A Tree with Muon Data"); // crate TTree inside the file
    float * null_float = nullptr; // create a null pointer to float to be used in the branches 

    TH2F *histogram = new TH2F("PT_vs_P_histogram", "Muon_PT_vs_P", 100, 0, 500, 100, 0, 20); // Example: 100 bins from 0 to 10 in both dimensions
    histogram->SetDirectory(outputFile); // Set the histogram to be saved in the file
    
    muon_tree.Branch("mu_p_PX", null_float); // create branches for the tree, now there are empty, but will be filled later
    muon_tree.Branch("mu_p_PY", null_float);
    muon_tree.Branch("mu_p_PZ", null_float);
    muon_tree.Branch("mu_p_ENERGY", null_float);
    muon_tree.Branch("mu_m_PX", null_float);
    muon_tree.Branch("mu_m_PY", null_float);
    muon_tree.Branch("mu_m_PZ", null_float);
    muon_tree.Branch("mu_m_ENERGY", null_float);

    muon_tree.SetDirectory(outputFile); // Set the tree to be saved in the file

    const char* filename = "/root/root_folder/pythia8311/data/ScalarBoson_2500MeV_400ps.root"; // path to the file with the Scalar data

    int idGun = 25; // Higgs boson
    double eeGun = 150; // energy of the Higgs boson
    bool hasLifetime = true; // if the particle has a lifetime

    int nEvent = 1e2; // Number of events to generate
    int nList = 1; // Number of events to list

    Event& event = pythia.event; // Declare event object
    ParticleData& pdt = pythia.particleData; // Declare particle data object

    pythia.readString("ProcessLevel:all = off"); // Switch off all process-level settings.

    pythia.readString("25:m0 = 2.5"); // Higgs mass.
    pythia.readString("25:onMode = off"); // Switch off all Higgs decays.
    pythia.readString("25:oneChannel = 1 1 21 13 -13"); // Switch on Higgs -> gg
    pythia.readString("Random:setSeed = on"); // Set seed for random number generator
    pythia.readString("Random:seed = 0"); // Set seed for random number generator

    TFile file(filename); 
    if (file.IsZombie()) { // Check if the file is open
        std::cerr << "Error: Failed to open file " << filename << std::endl;
        return 0;
    }

    // Get the tree from the file
    TTree *tree = dynamic_cast<TTree*>(file.Get("MCTuple/DecayTree")); // replace "tree_name" with the actual name of the tree
    if (!tree) {
        std::cerr << "Error: Failed to retrieve tree from file." << std::endl;
        file.Close();
        return 0;
    }


    
	pythia.init(); // Initialize.


    for (int iEvent = 0; iEvent < nEvent; ++iEvent) { // Event loop.
        fillParticle(pythia, idGun, eeGun, 5., 1., event, pdt, pythia.rndm, tree, &muon_tree, histogram, hasLifetime); // Fill the event with the particle

        if (iEvent % 50 == 0) { // Print the event number
            cout << "Event " << iEvent << " of " << nEvent << endl;
        }
      
    }
    histogram->Write(); // Write the histogram to the file
    muon_tree.Write(); // Write the tree to the file
    outputFile->Write(); // Write the file
    outputFile->Close(); // Close the output file
    file.Close(); // Close the input file

    pythia.stat();
	return 0;
}
