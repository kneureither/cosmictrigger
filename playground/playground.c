#include <TFile.h>
#include <TROOT.h>
#include <algorithm>
#include <TTree.h>
#include <map>
#include <string>
#include <TH1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <fstream>
#include "../util/utility_functions.h"
#include "../util/custom_types.h"


// using namespace std;

using std::cout;
using std::endl;

void playground() {
    const float LEFT_BOUNDARY = 0;
    const float RIGHT_BOUNDARY = 4;
    const int BIN_COUNT = 10;

    TH1F *h2 = new TH1F("h2", "1/p reconstruction errors", BIN_COUNT, LEFT_BOUNDARY, RIGHT_BOUNDARY);
    for (int i = 0; i < 500; i++) {
        h2->Fill(1.0);
    }
    auto  *c2 = new TCanvas();
    c2->SetWindowPosition(0, 500 );

    h2->DrawClone("");
}

void playground4() {
    std::vector<int> myvec;
    std::vector<float> myfloatvec;
    srand(time(NULL));

    int sum = 0;
    float sumfloat = 0.0;

    for(int i=0; i < 20; i++) {
        int rand = std::rand();
        float float_rand = ((float) std::rand()) / (float) RAND_MAX;
        myvec.push_back(rand);
        myfloatvec.push_back(float_rand);

        sum += rand;
        sumfloat += float_rand;
    }

    cout << "FLOAT:\tVEC: " << vector_mean(myfloatvec) << " CALC: " << sumfloat / 20.0 << endl;
    cout << "INT:\tVEC: " << vector_mean(myvec) << " CALC: " << sum / 20.0 << endl;
}

void playground2() {
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(111111);
    gStyle->SetOptFit(1111);
    gStyle->SetPalette(57);
    gStyle->SetOptTitle(0);

    std::vector<float> myarray;
    for(int i=0; i < 100; i++) {
        float mynumber = ((float) rand() / RAND_MAX);
        myarray.push_back(mynumber);
        cout << mynumber << endl;
    }

    TH1F h("h", "example", 10, 0., 1.);
    for(int i=0; i < myarray.size(); i++) {
        h.Fill(myarray[i]);
    }

    auto  mycanvas = new TCanvas();
    mycanvas->SetWindowPosition(0, 400 );

    h.SetMarkerStyle(kOpenCircle);
    h.SetMarkerColor(kBlue);
    h.SetLineColor(kBlue);

//    graph.DrawClone("APE");
    h.DrawClone("");
    h.Print();
}

void playground1() {
    TFile f1("../data/mu3e_run_000006.root");
    TTree *t_mu3e;
    f1.GetObject("mu3e", t_mu3e);

    TFile f2("../data/mu3e_run_000006_trirec_cosmic2.root");
    TTree *t_segs;
    f2.GetObject("segs", t_segs);


    //data for mu3e tree in mu3eSim
    int nhits;
    int header[4];
    unsigned int entries = t_mu3e->GetEntries();

    std::vector<int> * pixelid =nullptr;
    std::vector<int> * trajids =nullptr;
    std::vector<int> * timestamp =nullptr;
    std::vector<double> * traj_px = nullptr;
    std::vector<double> * traj_py = nullptr;
    std::vector<double> * traj_pz = nullptr;

    t_mu3e->SetBranchAddress("Header", &header);
    t_mu3e->SetBranchAddress("Nhit", &nhits);

    t_mu3e->SetBranchAddress("hit_pixelid",&pixelid);
    t_mu3e->SetBranchAddress("hit_timestamp", &timestamp);
    t_mu3e->SetBranchAddress("traj_ID", &trajids);

    t_mu3e->SetBranchAddress("traj_px", &traj_px);
    t_mu3e->SetBranchAddress("traj_py", &traj_py);
    t_mu3e->SetBranchAddress("traj_pz", &traj_pz);


    //data for trirec result tree segs
    int rec_event;
    int rec_nhit;
    int rec_tid;

    float mc_p;

    float rec_p;
    float rec_r;

    t_segs->SetBranchAddress("eventId", &rec_event);
    t_segs->SetBranchAddress("nhit", &rec_nhit);
    t_segs->SetBranchAddress("mc_tid", &rec_tid);
    t_segs->SetBranchAddress("mc_p", &mc_p);
    t_segs->SetBranchAddress("p", &rec_p);
    t_segs->SetBranchAddress("r", &rec_r);
    unsigned int segs_entry = 1;
    unsigned int segs_entries = t_segs->GetEntries();

    // cout << "segs_entries = " << segs_entries << endl;

    std::vector<float> rec_errors;
    std::vector<int> rec_nhits;

    cout << std::endl;
    cout << "Entry\tevent\tNhit\tSensor\tcolumn\trow\ttrajid\t\trec id\tnhit\t";
    cout << "trajid\tp_rec\tp_mc" << endl;

    t_segs->GetEntry(segs_entry);

    for (unsigned int entry=0; entry < entries; entry++) {
        t_mu3e->GetEntry(entry);
        while(rec_event < header[0] && segs_entry < segs_entries) {
            segs_entry++;
            t_segs->GetEntry(segs_entry);
        }

        for(unsigned int hit=0; hit < (unsigned int)nhits; hit++) {
            unsigned int sensor     = (*pixelid)[hit]>>16;
            unsigned int column     = ((*pixelid)[hit]>>8)&0xFF;
            unsigned int row        = ((*pixelid)[hit])&0xFF;
            unsigned int columnaddr = column + (sensor << 8);

            unsigned int trajid = (*trajids)[hit];

            // cout << entry << "\t" << header[0] << "\t" << nhits << "\t" << sensor << "\t";
            // cout << column << "\t" << row << "\t" << (*trajids)[0];

            // if(rec_event == header[0]) {
            //   cout << entry << "\t" << header[0] << "\t" << nhits << "\t" << sensor << "\t";
            //   cout << column << "\t" << row << "\t" << (*trajids)[0] << "\t|";
            //
            //   cout << "\t" << rec_tid;
            //   cout << "\t" << rec_event;
            //   cout << "\t" << rec_nhit;
            //   // cout << "\t" << rec_p;
            //   // cout << "\t" << mc_p;
            //   cout << "\t" << abs(abs(rec_p)-abs(mc_p)) / mc_p * 100 << "%";
            //   cout << endl;
            // }

        }

        if(rec_event == header[0]) {

          float rec_error = abs(abs(rec_p)-abs(mc_p)) / mc_p;
          rec_errors.push_back(rec_error);
          rec_nhits.push_back(rec_nhit);

          cout << rec_tid;
          cout << "\t" << rec_event;
          cout << "\t" << rec_nhit;
          // cout << "\t" << rec_p;
          // cout << "\t" << mc_p;
          cout << "\t" << rec_error * 100 << "%";
          cout << endl;
        }

        if((unsigned int) nhits == 0) {
            // cout << entry << "\t" << header[0] << "\t" << nhits << "\t" << "\t";
            // cout << "\t" << "\t" << (*trajids)[0] << endl;
        }
    }

    float rec_error_mean = std::accumulate(rec_errors.begin(), rec_errors.end(),
                                            decltype(rec_errors)::value_type(0))
                                            / rec_errors.size();

    cout << endl << "the mean of errors is " << rec_error_mean * 100  <<"%" << endl;
}



// void extra() {
//   TFile in_file("mu3e_run_000006.root");
//   // TFile in_file("conductivity_experiment.root");
//   TTree* my_tuple;
//   in_file.GetObject("mu3e", my_tuple);
//
//   float weight, Nhit, hit_pixelid, traj_ID;
//   float* row_content;
//
//   cout << "Nhit\tweight\ttraj_ID\thit_pixelid" << endl;
//   for(int irow=0; irow < my_tuple->GetEntries();++irow) {
//   // for(int irow=0; irow < 2;++irow) {
//     my_tuple -> GetEntry(irow);
//     row_content = my_tuple->GetArgs();
//     Nhit = row_content[3];
//     cout << Nhit << "\t" << "noch mehr..." << endl;
//   }
// }


// void extra2() {
//   TChain in_chain("mu3e");
//   in_chain.Add("mu3e_run_000006*.root");
//
//   int nhits;
//   int events;
//   std::vector<int> * pixelid =0;
//     std::vector<int> * pixelmcid =0;
//
//
//   in_chain.SetBranchAddress("Nhit", &Nhit);
//   in_chain.SetBranchAddress("traj_ID.traj_ID", &traj_ID);
//   // in_chain.SetBranchAddress("hit_pixelid", &hit_pixelid);
//
//
//   cout << "Index\tNhit\ttraj_ID" << endl;
//
//   Int_t max_row = -1;
//
//   for (size_t irow=0; irow < in_chain.GetEntries(); ++irow) {
//     in_chain.GetEntry(irow);
//
//
//     if(Nhit >= 4) {
//       cout <<irow << "\t"<< Nhit << "\t" << traj_ID<< endl;
//     }
//
//     max_row = irow;
//   }
//
//   cout << "max row = "<< max_row << endl;
//
//
// }
