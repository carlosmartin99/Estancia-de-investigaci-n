#include "Pythia8/Pythia.h" // Include Pythia headers.
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include <TGraph.h>
#include <iostream>
using namespace Pythia8;

// This functions allows me to save the acceptance graph in output.file
void saveGraphToRootFile(TGraph *graph, const char *filename)
{
    // Open the existing ROOT file in "UPDATE" mode
    TFile *file = new TFile(filename, "UPDATE");

    // Write the TGraph object to the file
    graph->Write();

    // Close the file
    file->Close();

    delete file; // Delete the TFile object to free memory
}

// Function that gives true is the particle is detected and false if not
bool Detector(double vert_x, double vert_y, double vert_z, double px, double py, double pz)
{
    double posx, posy, posz, x, y;
    // Detector positions
    posz = 9415;
    posy = 2400;
    posx = 540 * 6;
    // Trajectory of my particles (z is fixed as the value where the detector is) (Define to work in any unit as px/py is dimensionless)
    x = vert_x + (px / pz) * (posz - vert_z);
    y = vert_y + (py / pz) * (posz - vert_z);

    if ((x > posx || x < -posx) && (y > posy || y < -posy))
    {
        // std::cout<< "NOT DETECTED"<<std::endl;
        return false;
    }
    else
    {
        // std::cout<< "DETECTED"<<std::endl;
        return true;
    }
}

void FillMuonData(Event &event, int iNew, TTree *outputTree, TH2F *histogram, TH1F *h1, TH1F *h2,TH1F *h3,TH1F *h4,TH1F *h5, TH1F *h6, Rndm &rndm)
{
    // This function receives an object from event class (a vector which stores the variety of events)
    Particle &particle = event[iNew];

    // Check if the particle has daughters
    if (particle.daughterList().size() > 0)
    {
        // MUONS
        float mu_plus_px, mu_plus_py, mu_plus_pz, mu_plus_energy; // empty variable to be filled a bit later
        float mu_minus_px, mu_minus_py, mu_minus_pz, mu_minus_energy;
        float mu_plus_pt, mu_minus_pt, mu_plus_p, mu_minus_p;

        // PIONS
        float pi_plus_px, pi_plus_py, pi_plus_pz, pi_plus_energy;
        float pi_minus_px, pi_minus_py, pi_minus_pz, pi_minus_energy;
        float pi_plus_pt, pi_minus_pt, pi_plus_p, pi_minus_p;

        // KAONS
        float kaon_plus_px, kaon_plus_py, kaon_plus_pz, kaon_plus_energy;
        float kaon_minus_px, kaon_minus_py, kaon_minus_pz, kaon_minus_energy;
        float kaon_plus_pt, kaon_minus_pt, kaon_plus_p, kaon_minus_p;

        // MY AUXILIARY VARIABLES
        float eventos;
        float vert_x, vert_y, vert_z;
        bool detection;
        float posz = 9415;

        outputTree->SetBranchAddress("mu_p_PX", &mu_plus_px); // path where to be filled
        outputTree->SetBranchAddress("mu_p_PY", &mu_plus_py);
        outputTree->SetBranchAddress("mu_p_PZ", &mu_plus_pz);
        outputTree->SetBranchAddress("mu_p_ENERGY", &mu_plus_energy);
        outputTree->SetBranchAddress("mu_m_PX", &mu_minus_px);
        outputTree->SetBranchAddress("mu_m_PY", &mu_minus_py);
        outputTree->SetBranchAddress("mu_m_PZ", &mu_minus_pz);
        outputTree->SetBranchAddress("mu_m_ENERGY", &mu_minus_energy);
        outputTree->SetBranchAddress("daughter_p_pid", &eventos); // This is only used to see that the branching ratio imposed is satisfied

        outputTree->SetBranchAddress("pi_p_PX", &pi_plus_px); // path where to be filled
        outputTree->SetBranchAddress("pi_p_PY", &pi_plus_py);
        outputTree->SetBranchAddress("pi_p_PZ", &pi_plus_pz);
        outputTree->SetBranchAddress("pi_p_ENERGY", &pi_plus_energy);
        outputTree->SetBranchAddress("pi_m_PX", &pi_minus_px);
        outputTree->SetBranchAddress("pi_m_PY", &pi_minus_py);
        outputTree->SetBranchAddress("pi_m_PZ", &pi_minus_pz);
        outputTree->SetBranchAddress("pi_m_ENERGY", &pi_minus_energy);

        outputTree->SetBranchAddress("kaon_p_PX", &kaon_plus_px); // path where to be filled
        outputTree->SetBranchAddress("kaon_p_PY", &kaon_plus_py);
        outputTree->SetBranchAddress("kaon_p_PZ", &kaon_plus_pz);
        outputTree->SetBranchAddress("kaon_p_ENERGY", &kaon_plus_energy);
        outputTree->SetBranchAddress("kaon_m_PX", &kaon_minus_px);
        outputTree->SetBranchAddress("kaon_m_PY", &kaon_minus_py);
        outputTree->SetBranchAddress("kaon_m_PZ", &kaon_minus_pz);
        outputTree->SetBranchAddress("kaon_m_ENERGY", &kaon_minus_energy);
        for (int j = 0; j < particle.daughterList().size(); ++j)
        {
            int daughterIndex = particle.daughterList()[j];
            Particle &daughter = event[daughterIndex];
            // std::cout<<"Madre:   "<<daughter.id()<<"  "<<daughter.zProd()<<endl;
            // std::cout << "  Daughter " << j << ": ID = " << daughter.id() << ", Status = " << daughter.status() << std::endl;
            if (daughter.id() == 35)
            {
                // Before turning Higgs into my primary I check the vertex coordinates where produced
                // as I need to calculate manually the decay vertex of my modified Higgs.
                // float tau0=300;
                // particle.tau(tau0 *rndm.exp());
                // I also define the mass arbitrarily because it does not change
                float mass = daughter.m0();
                // After this we have obtained tau modified as expected, but tau0 remains the same (should not matter i think)
                vert_x = daughter.px() * daughter.tau() / mass;
                vert_y = daughter.py() * daughter.tau() / mass;
                vert_z = daughter.pz() * daughter.tau() / mass;
                std::cout << "Tau    :" << daughter.tau() << std::endl;
                std::cout << "Mass   : " << mass << std::endl;
                std::cout << "Vertex x:  " << vert_x << std::endl;
                std::cout << "Vertex y:  " << vert_y << std::endl;
                std::cout << "Vertex z:  " << vert_z << std::endl;

                Particle &particle = daughter; // Machaco de forma que ahora la hija(Higgs) para a ser la particula referencia

                // Loop over Higgs daughters
                for (int j = 0; j < particle.daughterList().size(); ++j)
                {

                    int daughterIndex = particle.daughterList()[j]; // get daughter particle index

                    Particle &daughter = event[daughterIndex]; // get daughter particle

                    if (daughter.id() == 13)
                    {
                        mu_minus_px = daughter.px() * 1e3; // in pythia, the units are GeV, so we convert them to MeV
                        mu_minus_py = daughter.py() * 1e3;
                        mu_minus_pz = daughter.pz() * 1e3;
                        mu_minus_energy = daughter.e() * 1e3;
                        mu_minus_pt = sqrt(daughter.px() * daughter.px() + daughter.py() * daughter.py()) * 1e3;
                        mu_minus_p = sqrt(daughter.px() * daughter.px() + daughter.py() * daughter.py() + daughter.pz() * daughter.pz()) * 1e3;
                        eventos = 13;
                        histogram->Fill(mu_minus_p, mu_minus_pt); // fill the histogram with the data
                        // This condition is made so we dont count as detection a particle which decayed after the detector
                        if (posz > vert_z)
                        {
                            h1->Fill(vert_z, 1);
                        }
                        // tau=daughter.tau();
                        // m=daughter.m0();
                        detection = Detector(vert_x, vert_y, vert_z, mu_minus_px, mu_minus_py, mu_minus_pz);
                        if ((detection == true) && (posz > vert_z))
                        {
                            h2->Fill(vert_z, 1);
                        }
                    }
                    else if (daughter.id() == -13)
                    { // if the daughter is a muon plus do same
                        mu_plus_px = daughter.px() * 1e3;
                        mu_plus_py = daughter.py() * 1e3;
                        mu_plus_pz = daughter.pz() * 1e3;
                        mu_plus_energy = daughter.e() * 1e3;
                        mu_plus_pt = sqrt(daughter.px() * daughter.px() + daughter.py() * daughter.py()) * 1e3;
                        mu_plus_p = sqrt(daughter.px() * daughter.px() + daughter.py() * daughter.py() + daughter.pz() * daughter.pz()) * 1e3;
                        eventos = -13;
                        histogram->Fill(mu_plus_p, mu_plus_pt);
                        if (posz > vert_z)
                        {
                            h1->Fill(vert_z, 1);
                        }
                        detection = Detector(vert_x, vert_y, vert_z, mu_plus_px, mu_plus_py, mu_plus_pz);
                        if ((detection == true) && (posz > vert_z))
                        {
                            h2->Fill(vert_z, 1);
                        }
                    }
                    else if (daughter.id() == 211)
                    { // if the daughter is a muon plus do same
                        pi_plus_px = daughter.px() * 1e3;
                        pi_plus_py = daughter.py() * 1e3;
                        pi_plus_pz = daughter.pz() * 1e3;
                        pi_plus_energy = daughter.e() * 1e3;
                        pi_plus_pt = sqrt(daughter.px() * daughter.px() + daughter.py() * daughter.py()) * 1e3;
                        pi_plus_p = sqrt(daughter.px() * daughter.px() + daughter.py() * daughter.py() + daughter.pz() * daughter.pz()) * 1e3;
                        eventos = 211;
                        if (posz > vert_z)
                        {
                            h3->Fill(vert_z, 1);
                        }
                        detection = Detector(vert_x, vert_y, vert_z, pi_plus_px, pi_plus_py, pi_plus_pz);
                        if ((detection == true) && (posz > vert_z))
                        {
                            h4->Fill(vert_z, 1);
                        }
                    }
                    else if (daughter.id() == -211)
                    { // if the daughter is a muon plus do same
                        pi_minus_px = daughter.px() * 1e3;
                        pi_minus_py = daughter.py() * 1e3;
                        pi_minus_pz = daughter.pz() * 1e3;
                        pi_minus_energy = daughter.e() * 1e3;
                        pi_minus_pt = sqrt(daughter.px() * daughter.px() + daughter.py() * daughter.py()) * 1e3;
                        pi_minus_p = sqrt(daughter.px() * daughter.px() + daughter.py() * daughter.py() + daughter.pz() * daughter.pz()) * 1e3;
                        eventos = -211;
                         if (posz > vert_z)
                        {
                            h3->Fill(vert_z, 1);
                        }
                        detection = Detector(vert_x, vert_y, vert_z, pi_minus_px, pi_minus_py, pi_minus_pz);
                        if ((detection == true) && (posz > vert_z))
                        {
                            h4->Fill(vert_z, 1);
                        }
                    }
                    else if (daughter.id() == 321)
                    { // if the daughter is a muon plus do same
                        kaon_plus_px = daughter.px() * 1e3;
                        kaon_plus_py = daughter.py() * 1e3;
                        kaon_plus_pz = daughter.pz() * 1e3;
                        kaon_plus_energy = daughter.e() * 1e3;
                        kaon_plus_pt = sqrt(daughter.px() * daughter.px() + daughter.py() * daughter.py()) * 1e3;
                        kaon_plus_p = sqrt(daughter.px() * daughter.px() + daughter.py() * daughter.py() + daughter.pz() * daughter.pz()) * 1e3;
                        eventos = 321;
                          if (posz > vert_z)
                        {
                            h5->Fill(vert_z, 1);
                        }
                        detection = Detector(vert_x, vert_y, vert_z, pi_minus_px, pi_minus_py, pi_minus_pz);
                        if ((detection == true) && (posz > vert_z))
                        {
                            h6->Fill(vert_z, 1);
                        }
                    }
                    else if (daughter.id() == -321)
                    { // if the daughter is a muon plus do same
                        kaon_minus_px = daughter.px() * 1e3;
                        kaon_minus_py = daughter.py() * 1e3;
                        kaon_minus_pz = daughter.pz() * 1e3;
                        kaon_minus_energy = daughter.e() * 1e3;
                        kaon_minus_pt = sqrt(daughter.px() * daughter.px() + daughter.py() * daughter.py()) * 1e3;
                        kaon_minus_p = sqrt(daughter.px() * daughter.px() + daughter.py() * daughter.py() + daughter.pz() * daughter.pz()) * 1e3;
                        eventos = -321;
                          if (posz > vert_z)
                        {
                            h5->Fill(vert_z, 1);
                        }
                        detection = Detector(vert_x, vert_y, vert_z, pi_minus_px, pi_minus_py, pi_minus_pz);
                        if ((detection == true) && (posz > vert_z))
                        {
                            h6->Fill(vert_z, 1);
                        }
                    }
                }
            }
        }

        outputTree->Fill(); // filling output tree
    }
    // }
}
void fillParticle(auto &pythia, int id, double ee, double thetaIn, double phiIn,
                  Event &event, ParticleData &pdt, Rndm &rndm, TTree *dataTree, TTree *outputTree, TH2F *histogram, TH1F *h1, TH1F *h2,TH1F *h3,TH1F *h4,TH1F *h5, TH1F *h6, bool hasLifetime = false)
{

    bool atRest = false;
    // Reset event record to allow for new event.
    event.reset(); // Reset the event record.

    // Select particle mass; where relevant according to Breit-Wigner.
    double mm = pdt.mSel(id);

    int numberOfEntries = dataTree->GetEntries(); // getting total number of events in ROOT FILE (en scalar.root)

    int entry = int(rndm.flat() * numberOfEntries); // get random event from the file
    dataTree->GetEntry(entry);                      // get the event

    float px, py, pz, higgs_energy; // empty variables to be filled a bit later
    // dataTree->SetBranchAddress("Scalar_PX",&px); // read px value from input file (this was inputf file of ScalarBoson)
    // dataTree->SetBranchAddress("Scalar_PY",&py);
    // dataTree->SetBranchAddress("Scalar_PZ",&pz);
    dataTree->SetBranchAddress("B_plus_PX", &px); // read px value from input file (input where B-->ScalarBoson-->...)
    dataTree->SetBranchAddress("B_plus_PY", &py);
    dataTree->SetBranchAddress("B_plus_PZ", &pz);

    float px_gev = px / 1e3; // convert MeV to GeV
    float py_gev = py / 1e3; // convert MeV to GeV
    float pz_gev = pz / 1e3;
    float tau0;
    // La lectura del primer evento daba error (variables desmesuradas)
    // Hemos metido una cota para que estos eventos corruptos no sean procesados y asi haya decays
    // There was an error only on first event, where it diverges
    // We stablish a threshold that eliminates this kind of problem(actually only happenend for first event)
    if (std::abs(pz_gev) > 1e5)
    {
        return;
    }

    double calculated_energy = sqrt(px_gev * px_gev + py_gev * py_gev + pz_gev * pz_gev + mm * mm);

    int iNew = event.append(id, 1, 0, 0, px_gev, py_gev, pz_gev, calculated_energy, mm); // Add the particle to the event record
    tau0 = 300;                                                                          // We force our decay wanted because pythia standard declaration did not work
    if (hasLifetime)
        event[iNew].tau(tau0 * rndm.exp());                                // Set the lifetime of the particle
    pythia.next();                                                         // process the event -> calculate all decays and kinematics
    FillMuonData(event, iNew, outputTree, histogram, h1, h2, h3, h4,h5,h6, pythia.rndm); // fill the output tree with muon data
}

int main()
{
    Pythia pythia; // Declare Pythia object

    TFile *outputFile = TFile::Open("output_3000Mev_250mm.root", "RECREATE"); // create file, where we will write data from PYTHIA

    TTree muon_tree("MuonsTree", "A Tree with Muon Data"); // crate TTree inside the file
    float *null_float = nullptr;                           // create a null pointer to float to be used in the branches

    TH2F *histogram = new TH2F("PT_vs_P_histogram", "Muon_PT_vs_P", 100, 0, 500e3, 100, 0, 20e3); // Example: 100 bins from 0 to 10 in both dimensions
    histogram->SetDirectory(outputFile);

    TH1F *h1 = new TH1F("h1", "Higgs decays into muons", 50, 0, 8000); // Histogram for muons events
    h1->SetDirectory(outputFile);

    TH1F *h2 = new TH1F("h2", "Muons decayed measured by detector", 50, 0, 8000); // Histogram for detected muon events
    h2->SetDirectory(outputFile);

    TH1F *h3 = new TH1F("h3", "Higgs decays into pions", 50, 0, 8000); // Histogram for pion events
    h3->SetDirectory(outputFile);

    TH1F *h4 = new TH1F("h4", "Pions decayed measured by detector", 50, 0, 8000); // Histogram for detected pion events
    h4->SetDirectory(outputFile);

    TH1F *h5 = new TH1F("h5", "Higgs decays into kaons", 50, 0, 8000); // Histogram for kaon events
    h3->SetDirectory(outputFile);

    TH1F *h6 = new TH1F("h6", "Kaons decayed measured by detector", 50, 0, 8000); // Histogram for detected kaon events
    h4->SetDirectory(outputFile);
    // TH1F *acceptance = new TH1F("acceptance", "Acceptance",30,0,150);
    // acceptance->SetDirectory(outputFile);

    // MUON
    muon_tree.Branch("daughter_p_pid", null_float);
    muon_tree.Branch("mu_p_PX", null_float); // create branches for the tree, now there are empty, but will be filled later
    muon_tree.Branch("mu_p_PY", null_float);
    muon_tree.Branch("mu_p_PZ", null_float);
    muon_tree.Branch("mu_p_ENERGY", null_float);
    muon_tree.Branch("mu_m_PX", null_float);
    muon_tree.Branch("mu_m_PY", null_float);
    muon_tree.Branch("mu_m_PZ", null_float);
    muon_tree.Branch("mu_m_ENERGY", null_float);

    // PION
    muon_tree.Branch("pi_p_PX", null_float); // create branches for my pion tree
    muon_tree.Branch("pi_p_PY", null_float);
    muon_tree.Branch("pi_p_PZ", null_float);
    muon_tree.Branch("pi_p_ENERGY", null_float);
    muon_tree.Branch("pi_m_PX", null_float);
    muon_tree.Branch("pi_m_PY", null_float);
    muon_tree.Branch("pi_m_PZ", null_float);
    muon_tree.Branch("pi_m_ENERGY", null_float);

    // KAON
    muon_tree.Branch("kaon_p_PX", null_float); // create branches for my pion tree
    muon_tree.Branch("kaon_p_PY", null_float);
    muon_tree.Branch("kaon_p_PZ", null_float);
    muon_tree.Branch("kaon_p_ENERGY", null_float);
    muon_tree.Branch("kaon_m_PX", null_float);
    muon_tree.Branch("kaon_m_PY", null_float);
    muon_tree.Branch("kaon_m_PZ", null_float);
    muon_tree.Branch("kaon_m_ENERGY", null_float);

    muon_tree.SetDirectory(outputFile); // Set the tree to be saved in the file
    // Definition of the filename from which i read data. Its read using TFile (from ROOT)
    // const char* filename = "/root/root_folder/pythia8311/data/ScalarBoson_2500MeV_400ps.root"; // path to the file with the Scalar data (old)
    const char *filename = "/pythia8311/examples/scalar_boson/H2mumu_3000MeV_100ps.root"; //(new)
    int idGun = 521;                                                                      // B+ meson boson
    double eeGun = 150;                                                                   // energy of the B+ meson
    bool hasLifetime = true;                                                              // if the particle has a lifetime

    int nEvent = 1e6;               // Number of events to generate
    int nList = 2;                  // Number of events to list
    float higgs1, higgs2, cociente; // variables to construct acceptance

    Event &event = pythia.event;             // Declare event object
    ParticleData &pdt = pythia.particleData; // Declare particle data object

    pythia.readString("ProcessLevel:all = off");
    // IMPORTANT: Use 35, as 25 is strictly reserved in Pythia and did not allow to change particle properties
    // B+ decaying into modified Higgs (scalar boson) plus K+
    pythia.readString("521:oneChannel = 1 1 21 35 321");
    // Higgs modification: whats not zero is mass, min_mass, max_mass and tau_0 (in mm)
    pythia.readString("35:new  = Higgs0 void 0 0 0 3.0 0.0 3.0 3.0 250");
    // Random seeds (crucial)
    pythia.readString("Random:setSeed = on"); // Set seed for random number generator
    pythia.readString("Random:seed = 0");     // Set seed for random number generator
    // Eliminate all decay channels to manually add the interesting one for my study
    pythia.readString("35:onMode = off");
    // pythia.readString("35:oneChannel = 1 1 1 13 -13");

    // Channel MU MU
    pythia.readString("35:addChannel = 1 0.009 1 13 -13");
    // Channel PI PI
    pythia.readString("35:addChannel = 1 0.81 1 211 -211");
    // Channel K K
    pythia.readString("35:addChannel = 1 0.16 1 321 -321");

    // pythia.readString("321:mayDecay = on");

    // id:all = name antiName spinType chargeType colType m0 mWidth mMin mMax tau0

    TFile file(filename);
    if (file.IsZombie())
    { // Check if the file is open
        std::cerr << "Error: Failed to open file " << filename << std::endl;
        return 0;
    }

    // Get the tree from the file
    TTree *tree = dynamic_cast<TTree *>(file.Get("MCTuple/DecayTree")); // replace "tree_name" with the actual name of the tree
    if (!tree)
    {
        std::cerr << "Error: Failed to retrieve tree from file." << std::endl;
        file.Close();
        return 0;
    }

    pythia.init(); // Initialize.

    for (int iEvent = 0; iEvent < nEvent; ++iEvent)
    {                                                                                                                          // Event loop.
        fillParticle(pythia, idGun, eeGun, 5., 1., event, pdt, pythia.rndm, tree, &muon_tree, histogram, h1, h2,h3,h4,h5,h6, hasLifetime); // Fill the event with the particle

        if (iEvent % 50 == 0)
        { // Print the event number
            cout << "Event " << iEvent << " of " << nEvent << endl;
        }
    }
    TGraph *acceptance = new TGraph();
    // Acceptance bin to bin
    int nBins = h1->GetNbinsX();
    TH1F *histograma = new TH1F("histograma", "Muon acceptance factor", nBins, h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
    histograma->SetDirectory(outputFile);
    for (int i = 0; i <= nBins; ++i)
    {
        float height = h1->GetBinCenter(i);
        float higgs1 = h1->GetBinContent(i);
        float higgs2 = h2->GetBinContent(i);
        float scaleFactor = (higgs2 != 0) ? higgs2 / higgs1 : 0;
        acceptance->SetPoint(i - 1, height, scaleFactor);
        histograma->SetBinContent(i, scaleFactor);
    }
    TH1F *histograma2 = new TH1F("histograma2", "Pion acceptance factor", nBins, h3->GetXaxis()->GetXmin(), h3->GetXaxis()->GetXmax());
    histograma2->SetDirectory(outputFile);
    for (int i = 0; i <= nBins; ++i)
    {
        float height = h3->GetBinCenter(i);
        float higgs1 = h3->GetBinContent(i);
        float higgs2 = h4->GetBinContent(i);
        float scaleFactor = (higgs2 != 0) ? higgs2 / higgs1 : 0;
        histograma2->SetBinContent(i, scaleFactor);
    }
    TH1F *histograma3 = new TH1F("histograma3", "Kaon acceptance factor", nBins, h5->GetXaxis()->GetXmin(), h5->GetXaxis()->GetXmax());
    histograma3->SetDirectory(outputFile);
    for (int i = 0; i <= nBins; ++i)
    {
        float height = h5->GetBinCenter(i);
        float higgs1 = h5->GetBinContent(i);
        float higgs2 = h6->GetBinContent(i);
        float scaleFactor = (higgs2 != 0) ? higgs2 / higgs1 : 0;
        histograma3->SetBinContent(i, scaleFactor);
    }
    //ACCEPTNACE MUON
    histograma->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
    histograma->GetYaxis()->SetTitle("Scale Factor");
    histograma->Write();
    //ACCEPTANCE PION
    histograma->GetXaxis()->SetTitle(h2->GetXaxis()->GetTitle());
    histograma->GetYaxis()->SetTitle("Scale Factor");
    histograma2->Write();
     //ACCEPTANCE KAON
    histograma->GetXaxis()->SetTitle(h3->GetXaxis()->GetTitle());
    histograma->GetYaxis()->SetTitle("Scale Factor");
    histograma3->Write();



    histogram->Write(); // Write the histogram to the file
    muon_tree.Write();  // Write the tree to the file
    acceptance->Write();
    outputFile->Write(); // Write the file
    outputFile->Close(); // Close the output file
    file.Close();        // Close the input file
    saveGraphToRootFile(acceptance, "output.root");

    pythia.stat();
    return 0;
}
