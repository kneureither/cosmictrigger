#ifndef COSMICTRIGGER_TPLOTS_H
#define COSMICTRIGGER_TPLOTS_H

#endif //COSMICTRIGGER_TPLOTS_H

#include <assert.h>
#include <TFile.h>
#include <TH1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <type_traits>
#include <iostream>

//This file contains functions for setting plot properties
//and filling plots with data.

using std::cout;
using std::endl;

void labelAxis(TH1 * h, const char * xtitle, const char * ytitle){
    h->GetXaxis()->SetTitle(xtitle);
    h->GetYaxis()->SetTitle(ytitle);

    h->GetXaxis()->SetTitleSize(0.04);
    h->GetYaxis()->SetTitleSize(0.04);

    h->GetXaxis()->SetTitleOffset(0.9);
    h->GetYaxis()->SetTitleOffset(1.2);

    h->GetXaxis()->SetLabelSize(0.03);
    h->GetYaxis()->SetLabelSize(0.03);

//    h->GetXaxis()->SetMaxDigits(5);
//    h->GetYaxis()->SetMaxDigits(5);
}

void labelAxis(TH1 * h, const char * xtitle, const char * ytitle, float ylabelsize){
    h->GetXaxis()->SetTitle(xtitle);
    h->GetYaxis()->SetTitle(ytitle);

    h->GetXaxis()->SetTitleSize(0.04);
    h->GetYaxis()->SetTitleSize(0.04);

    h->GetXaxis()->SetTitleOffset(0.9);
    h->GetYaxis()->SetTitleOffset(1.2);

    h->GetXaxis()->SetLabelSize(0.03);
    h->GetYaxis()->SetLabelSize(ylabelsize);

//    h->GetXaxis()->SetMaxDigits(5);
//    h->GetYaxis()->SetMaxDigits(5);
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

//    g->GetXaxis()->SetMaxDigits(5);
//    g->GetYaxis()->SetMaxDigits(5);

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

void saveCanvas(TCanvas* c1, const std::string& name, const std::string& path) {
    c1->SaveAs((path + "/" + name + ".eps").c_str());
    c1->SaveAs((path + "/" + name + ".png").c_str());
    c1->SaveAs((path + "/" + name + ".pdf").c_str());
}

template <typename RootGraph>
void makeSimpleMultiCanvas(int height, int width, int count, RootGraph ** g,
                           bool setlogy, bool setlogx, std::string plottingfile) {
    assert(height * width == count);
    TCanvas * canvas = new TCanvas("canvas", "canvas", width*300, (height > 1 ? height * 300 : 400));

    canvas->Divide(width, height);

    for (int i = 0; i < count; ++i) {
        canvas->cd(i+1);

        if(setlogx) gPad->SetLogx(1);
        if(setlogy) gPad->SetLogy(1);
        gPad->SetLeftMargin(0.15);

        if(std::is_same<RootGraph, TGraph>::value) {
            g[i]->Draw("ap");
        }
        if(std::is_same<RootGraph, TH1F>::value) {
            g[i]->Draw();
        }
    }
    canvas->Print(plottingfile.c_str(), "pdf");
}

template <typename RootGraph>
void makeSimpleMultiCanvas(int height, int width, int count, RootGraph * g, std::string plottingfile) {
    makeSimpleMultiCanvas(height, width, count, g, false, false, plottingfile);
}

template <typename RootGraph>
void makeSimpleSingleCanvas(RootGraph * g, bool setlogy, bool setlogx, std::string plottingfile) {
    TCanvas * canvas = new TCanvas("canvas", "canvas", 900, 900);

    if(setlogx) gPad->SetLogx(1);
    if(setlogy) gPad->SetLogy(1);
    gPad->SetLeftMargin(0.15);

    if(std::is_same<RootGraph, TGraph>::value) {
        g->Draw("ap");
    }
    if(std::is_same<RootGraph, TH1F>::value) {
        g->Draw();
    }
    canvas->Print(plottingfile.c_str(), "pdf");
}

template <typename RootGraph>
void makeSimpleSingleCanvas(RootGraph * g, std::string plottingfile) {
    makeSimpleSingleCanvas(g, false, false, plottingfile);
}