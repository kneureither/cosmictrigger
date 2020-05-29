#ifndef COSMICTRIGGER_TPLOTS_H
#define COSMICTRIGGER_TPLOTS_H

#endif //COSMICTRIGGER_TPLOTS_H

//This file contains functions for setting plot properties
//and filling plots with data.


void labelAxis(TH1 * h, const char * xtitle, const char * ytitle){
    h->GetXaxis()->SetTitle(xtitle);
    h->GetYaxis()->SetTitle(ytitle);

    h->GetXaxis()->SetTitleSize(0.04);
    h->GetYaxis()->SetTitleSize(0.04);

    h->GetXaxis()->SetTitleOffset(0.9);
    h->GetYaxis()->SetTitleOffset(1.2);

    h->GetXaxis()->SetLabelSize(0.03);
    h->GetYaxis()->SetLabelSize(0.03);
}

void labelAxis(TGraph * g, const char * xtitle, const char * ytitle) {
    g->GetXaxis()->SetTitle(xtitle);
    g->GetYaxis()->SetTitle(ytitle);

    g->GetXaxis()->SetTitleSize(0.04);
    g->GetYaxis()->SetTitleSize(0.04);

    g->GetXaxis()->SetTitleOffset(0.9);
    g->GetYaxis()->SetTitleOffset(1.8);

    g->GetXaxis()->SetLabelSize(0.03);
    g->GetYaxis()->SetLabelSize(0.03);
}

void setGraphRange(TGraph * g, float xrange_in, float xrange_out, float yrange_in, float yrange_out) {
    g->GetXaxis()->SetRangeUser(xrange_in, xrange_out);
    g->GetYaxis()->SetRangeUser(yrange_in, yrange_out);
}

template <typename T>
void fillHistWithVector(TH1F * h, std::vector<T> &data) {
    for (int i = 0; i < data.size(); i++) {
        h->Fill(data[i]);
    }
}

void SaveCanvas(TCanvas* c1, const std::string& name, const std::string& path) {
    c1->SaveAs((path + "/" + name + ".eps").c_str());
    c1->SaveAs((path + "/" + name + ".png").c_str());
    c1->SaveAs((path + "/" + name + ".pdf").c_str());
}