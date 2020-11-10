#ifndef COSMICTRIGGER_TPLOTS_H
#define COSMICTRIGGER_TPLOTS_H

#include <assert.h>
#include <TFile.h>
#include <TH1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <type_traits>
#include <iostream>
#include <TStyle.h>
#include <TMultiGraph.h>
#include <TLatex.h>

//This file contains functions for setting plot properties
//and filling plots with data.

using std::cout;
using std::endl;

static std::vector<int> setPlottingStyle() {
    gStyle->SetLabelFont(43);
    gStyle->SetLabelSize(14);
//    gStyle->SetLabelOffset(1.6);
    gStyle->SetTitleFont(53);
    gStyle->SetTitleSize(14);
    gStyle->SetLegendFont(42);
    gStyle->SetPadColor(0);
    gStyle->SetCanvasColor(0);

    gStyle->SetLegendBorderSize(0);
    gStyle->SetPalette(kThermometer);
    gStyle->SetMarkerStyle(23);

    std::vector<int> custom_color_palette = {};

    return custom_color_palette;
}

static void labelAxis(TH1 * h, const char * xtitle, const char * ytitle){
    h->GetXaxis()->SetTitle(xtitle);
    h->GetYaxis()->SetTitle(ytitle);

    h->GetXaxis()->SetTitleSize(0.04);
    h->GetYaxis()->SetTitleSize(0.04);

    h->GetXaxis()->SetTitleOffset(0.9);
    h->GetYaxis()->SetTitleOffset(0.9);

    h->GetXaxis()->SetLabelSize(0.03);
    h->GetYaxis()->SetLabelSize(0.03);

//    h->GetXaxis()->SetMaxDigits(5);
//    h->GetYaxis()->SetMaxDigits(5);
}

static void labelAxis(TH1 * h, const char * xtitle, const char * ytitle, float ylabelsize){
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

static void labelAxis(TGraph * g, const char * xtitle, const char * ytitle) {
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

static void setGraphRange(TGraph * g, float xrange_in, float xrange_out, float yrange_in, float yrange_out) {
    g->GetXaxis()->SetRangeUser(xrange_in, xrange_out);
    g->GetYaxis()->SetRangeUser(yrange_in, yrange_out);
}

template <typename T>
static void fillHistWithVector(TH1F * h, std::vector<T> &data) {
    for (int i = 0; i < data.size(); i++) {
        h->Fill(data[i]);
    }
}

static void saveCanvas(TCanvas* c1, const std::string& name, const std::string& path) {
//    c1->SaveAs((path + "/" + name + ".eps").c_str());
//    c1->SaveAs((path + "/" + name + ".png").c_str());
//    c1->SaveAs((path + "/" + name + ".pdf").c_str());
    c1->SaveAs((path + (path.substr(path.size()-1, 1) == "/" ? "" : "/") + name + ".pdf").c_str());
}

template <typename RootGraph>
static void makeSimpleMultiCanvas(int height, int width, int count, RootGraph ** g,
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
static void makeSimpleMultiCanvas(int height, int width, int count, RootGraph * g, std::string plottingfile) {
    makeSimpleMultiCanvas(height, width, count, g, false, false, plottingfile);
}

template <typename RootGraph>
static void makeSimpleSingleCanvas(RootGraph * g, bool setlogy, bool setlogx, std::string plottingfile) {
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
static void makeSimpleSingleCanvas(RootGraph * g, std::string plottingfile) {
    makeSimpleSingleCanvas(g, false, false, plottingfile);
}

static void setPlottingStyle(TH1F* hs) {
    hs->GetYaxis()->SetLabelFont(43);
    hs->GetYaxis()->SetLabelSize(16);
    hs->GetYaxis()->SetTitleFont(53);
    hs->GetYaxis()->SetTitleSize(17);
    hs->GetYaxis()->SetTitleOffset(1.4);
    hs->GetYaxis()->CenterTitle(false);

    hs->GetXaxis()->SetLabelFont(43);
    hs->GetXaxis()->SetLabelSize(16);
    hs->GetXaxis()->SetTitleFont(53);
    hs->GetXaxis()->SetTitleSize(17);
    hs->GetXaxis()->SetTitleOffset(1.5);
    hs->GetXaxis()->CenterTitle(false);
}

static void expandYaxisRange(TMultiGraph* g, float factor=0.1) {
    float min = g->GetHistogram()->GetMinimum();
    float max = g->GetHistogram()->GetMaximum();

    float range = max-min;
    g->SetMaximum(max + 0.5*factor*range);
    g->SetMinimum(min - 0.5*factor*range);
}

template <typename RootGraph>
static void expandYaxisRange(RootGraph* g, float factor=0.1) {
    float min = g->GetMinimum();
    float max = g->GetMaximum();

    float range = max-min;
    g->SetMaximum(max + 0.5*factor*range);
    g->SetMinimum(min - 0.5*factor*range);
}

static void drawAdditionalInfoBlock(TPad* pad, float x, float y, std::string l1, std::string l2, std::string l3="") {
    //not working!
    int font = 43;
    int textsize = 10;
    float spacing = 0.05;

    pad->cd();

    TLatex tline1(x,y,"hello");
    tline1.SetTextFont(font);
    tline1.SetTextSize(textsize);
    tline1.SetNDC(kTRUE);
    tline1.Draw();

//    TLatex tline2(x,y-spacing,l2.c_str());
//    tline2.SetTextFont(font);
//    tline2.SetTextSize(textsize);
//    tline2.SetNDC(kTRUE);
//    tline2.Draw();

//    TLatex tline3(x,y-2*spacing,l3.c_str());
//    tline3.SetTextFont(font);
//    tline3.SetTextSize(textsize);
//    tline3.SetNDC(kTRUE);
//    tline3.Draw();
}

#endif //COSMICTRIGGER_TPLOTS_H