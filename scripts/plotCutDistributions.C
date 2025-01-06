#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>

void plotCutDistributions() {
    // Open the root file
    gStyle->SetOptStat(0);
    TFile *file = TFile::Open("/pnfs/minerva/persistent/users/cpernas/default_analysis_loc/MC_TUPLE_Oct_23_2024_EMSHOWERSCORE.root");
    
    // Get the 8 histograms
    TH1F *h1 = (TH1F*)file->Get("EMScore_selected_signal_reco")->Clone("h1Clone");
    TH1F *h2 = (TH1F*)file->Get("EMScore_background_NuECC_QE_Like_no_proton")->Clone("h2Clone");;
    TH1F *h3 = (TH1F*)file->Get("EMScore_background_Other_NueCC")->Clone("h3Clone");
    TH1F *h4 = (TH1F*)file->Get("EMScore_background_nu___e_elastic")->Clone("h4Clone");
    TH1F *h5 = (TH1F*)file->Get("EMScore_background_NC_Coh")->Clone("h5Clone");
    TH1F *h6 = (TH1F*)file->Get("EMScore_background_Other_NC")->Clone("h6Clone");
    TH1F *h7 = (TH1F*)file->Get("EMScore_background_CC_Numu_pi0")->Clone("h7Clone");
    TH1F *h8 = (TH1F*)file->Get("EMScore_background_Other")->Clone("h8Clone");

    // Check if all histograms are retrieved
    if (!h1 || !h2 || !h3 || !h4 || !h5 || !h6 || !h7 || !h8) {
        printf("Error: Some histograms are missing.\n");
        return;
    }

    // Create a new histogram for the sum of the 7 histograms
    TH1F *hSum = (TH1F*)h2->Clone("hSum"); // Clone h2 to initialize
    hSum->Add(h3);
    hSum->Add(h4);
    hSum->Add(h5);
    hSum->Add(h6);
    hSum->Add(h7);
    hSum->Add(h8);

    // Create a canvas
    TCanvas *c1 = new TCanvas("c1", "EMShowerScore Cut", 800, 600);
    // Make room for the right axis
    gPad->SetRightMargin(0.2);

    // Set line color and style for h1 (red, not filled)
    h1->SetLineColor(kRed);
    h1->SetLineWidth(2);
    h1->SetFillStyle(0);

    // Set line color and style for hSum (blue, not filled)
    hSum->SetLineColor(kBlue);
    hSum->SetLineWidth(2);
    hSum->SetFillStyle(0);

    // Draw the first histogram (h1) with a normal y-axis on the left
    h1->Draw("HIST");

    // Get the range of the y-axis for h1 and adjust the axis range for hSum
    double max1 = h1->GetMaximum();
    double maxSum = hSum->GetMaximum();
    
    h1->GetYaxis()->SetTitle("Counts for h1");

    // Scale hSum to the same range as h1 (if necessary)
    hSum->Scale(max1 / maxSum);

    // Overlay the sum of the other histograms (hSum)
    hSum->Draw("HIST SAME");

    // Draw the second y-axis on the right side
    TGaxis *axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
                              gPad->GetUxmax(), gPad->GetUymax(),
                              0, maxSum, 510, "+L");
    axis->SetLineColor(kBlue);
    axis->SetLabelColor(kBlue);
    axis->SetTitle("Counts for hSum");
    axis->SetTitleOffset(1.5);
    axis->Draw();

    // Add a legend
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(h1, "Histogram 1 (Red)", "l");
    leg->AddEntry(hSum, "Sum of others (Blue, scaled)", "l");
    leg->Draw();

    // Save the canvas to a file (optional)
    c1->SaveAs("EMShowerScore_cut.png");

    // Clean up
    delete c1;
    delete h1;  // Free the cloned histograms
    delete h2;
    delete h3;
    delete h4;
    delete h5;
    delete h6;
    delete h7;
    delete h8;
    
    delete hSum;
    file->Close();
}
