//C++ INCLUDES
#include <iostream>
#include <map>
#include <vector>
//LOCAL INCLUDES
#include "PrintAIC.hh"
//ROOT INCLUDES
#include <TCanvas.h>
#include <TStyle.h>
#include <TPad.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TGraph.h>
//ROOFIT INCLUDES
#include <RooWorkspace.h>
#include <RooPlot.h>

void PrintAICTable(std::string MRcut,std::string RSQcut,std::map<std::string, double> delta_aic_map, std::map<std::string, double> delta_aic_map_2,std::map<std::string, double> delta_aic_map_3, std::map<std::string, double> aic_weight_map, std::map<std::string, double> aic_weight_map_2, std::map<std::string, double> aic_weight_map_3, RooWorkspace *w[]) 
{

	std::map< std::string, std::string > func_name;
	if( func_name.find("doubleExp") == func_name.end() ) func_name.insert( std::pair<std::string, std::string>("doubleExp","double exponential"));
	if( func_name.find("singleExp") == func_name.end() ) func_name.insert( std::pair<std::string, std::string>("singleExp","single exponential"));
	if( func_name.find("doublePow") == func_name.end() ) func_name.insert( std::pair<std::string, std::string>("doublePow","double power"));
	if( func_name.find("singlePow") == func_name.end() ) func_name.insert( std::pair<std::string, std::string>("singlePow","single power"));
	if( func_name.find("poly2") == func_name.end() ) func_name.insert( std::pair<std::string, std::string>("poly2","second order poly"));
	if( func_name.find("poly3") == func_name.end() ) func_name.insert( std::pair<std::string, std::string>("poly3","third order poly"));
	if( func_name.find("poly4") == func_name.end() ) func_name.insert( std::pair<std::string, std::string>("poly4","fourth order poly"));
	if( func_name.find("modExp") == func_name.end() ) func_name.insert( std::pair<std::string, std::string>("modExp","modified exponential"));

	std::map< std::string, int > num_par;
	if( num_par.find("doubleExp") == num_par.end() ) num_par.insert( std::pair<std::string, int>("doubleExp",3));
	if( num_par.find("singleExp") == num_par.end() ) num_par.insert( std::pair<std::string, int>("singleExp",1));
	if( num_par.find("doublePow") == num_par.end() ) num_par.insert( std::pair<std::string, int>("doublePow",3));
	if( num_par.find("singlePow") == num_par.end() ) num_par.insert( std::pair<std::string, int>("singlePow",1));
	if( num_par.find("poly2") == num_par.end() ) num_par.insert( std::pair<std::string, int>("poly2",3));
	if( num_par.find("poly3") == num_par.end() ) num_par.insert( std::pair<std::string, int>("poly3",4));
	if( num_par.find("poly4") == num_par.end() ) num_par.insert( std::pair<std::string, int>("poly4",5));
	if( num_par.find("modExp") == num_par.end() ) num_par.insert( std::pair<std::string, int>("modExp",2));

    //print the table to a file
	FILE* m_outfile = fopen("AIC_output/FitChoices_Table.tex", "a");
	fprintf(m_outfile,"\\begin{table}[H] \n");
	fprintf(m_outfile,"\\begin{center} \n");
	fprintf(m_outfile,"\\begin{tabular}{|c|c|cc|cc|cc|} \n");
	fprintf(m_outfile,"\\hline function & \\#P & $\\Delta A_1$ & $\\omega_1$ & $\\Delta A_2$ & $\\omega_2$ & $\\Delta A_3$ & $\\omega_3$ \\\\ \\hline \n");
	double delta_aic[8];
	double weight_aic[8];
	int idx=0;
	for ( auto tmp :func_name) 
      {
		  fprintf(m_outfile,"%s & %2d & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f \\\\ \n",func_name[tmp.first].c_str(), num_par[tmp.first], delta_aic_map[tmp.first], aic_weight_map[tmp.first], delta_aic_map_2[tmp.first], aic_weight_map_2[tmp.first], delta_aic_map_3[tmp.first], aic_weight_map_3[tmp.first]);
		  delta_aic[idx] = delta_aic_map[tmp.first];
		  weight_aic[idx] = aic_weight_map[tmp.first];
		  idx++;
	 }
	fprintf(m_outfile,"\\hline \n");
	fprintf(m_outfile,"\\end{tabular} \n");
	fprintf(m_outfile,"\\caption{$M_R$ cut = %s \\&\\& $R^2$ cut = %s.} \n", MRcut.c_str(),RSQcut.c_str());
	fprintf(m_outfile,"\\label{tab:FitChoices_%s_%s} \n", MRcut.c_str(),RSQcut.c_str());
	fprintf(m_outfile,"\\end{center} \n");
	fprintf(m_outfile,"\\end{table} \n \n \n");
	//pring the plot to a file
	FILE* m_outfile2 = fopen("AIC_output/FitChoices_Plot.tex", "a");
	fprintf(m_outfile2,"\\begin{figure}[H] \n");
	fprintf(m_outfile2,"\\begin{center} \n");
	fprintf(m_outfile2,"\\includegraphics[width=\\columnwidth]{Figure/%s_%s.pdf} \n", MRcut.c_str(),RSQcut.c_str());
	fprintf(m_outfile2,"\\caption{$M_R$ cut = %s \\&\\& $R^2$ cut = %s.} \n", MRcut.c_str(),RSQcut.c_str());
	fprintf(m_outfile2,"\\label{Fig:%s_%s} \n", MRcut.c_str(),RSQcut.c_str());
	fprintf(m_outfile2,"\\end{center} \n");
	fprintf(m_outfile2,"\\end{figure} \n \n \n");
	
	
	//plot the fit
	double delta_aic_2[8]={delta_aic[0],delta_aic[6],delta_aic[1],delta_aic[7],delta_aic[3],delta_aic[4],delta_aic[5],delta_aic[2]};
	double weight_aic_2[8]={weight_aic[0],weight_aic[6],weight_aic[1],weight_aic[7],weight_aic[3],weight_aic[4],weight_aic[5],weight_aic[2]};
	//double weight_aic_2[7]={weight_aic[0],weight_aic[5],weight_aic[1],weight_aic[6],weight_aic[3],weight_aic[4],weight_aic[2]};

	gStyle->SetTitleSize(0.15,"t");
	int n_g = 20;
        TGraph *grshade[8];
	
	TCanvas *myC= new TCanvas("c1","c1",200,10,800,600);
	myC->cd();
	gPad->SetLeftMargin(0.12);
	gPad->SetBottomMargin(0.18);

	myC->Divide(2,4,0);
	myC->cd(1);
	RooPlot * p1 = (RooPlot*)w[0]->obj("sideband_fit_doubleExp_pdf_dExp_N1N2_frame");
	p1->SetTitle("double exponential");
	p1->GetXaxis()->SetTitleSize(0.09);
	p1->GetXaxis()->SetLabelSize(0.11);
	p1->GetYaxis()->SetTitleSize(0.11);
	p1->GetYaxis()->SetLabelSize(0.11);
	p1->GetYaxis()->SetTitleOffset(0.55);
	p1->GetYaxis()->SetNdivisions(505);
	p1->Draw();
	
	double top_y = p1->GetMaximum();
	grshade[0] = new TGraph(2*(n_g+1));
	for ( int ig = 0; ig <= n_g; ig++ )
	  {
	    grshade[0]->SetPoint(ig, 120., (float)ig*(top_y/(float)n_g) );
	    grshade[0]->SetPoint( (n_g+1) + ig, 135., top_y - (float)ig*(top_y/(float)n_g) );
	  }
	//grshade->SetFillStyle(3013);
	grshade[0]->SetFillStyle(1001);  
	grshade[0]->SetFillColor(kRed-9);
	grshade[0]->Draw("f same");
	
	TLegend* leg[8];
	for(int idx=0;idx<8;idx++)
	{
	leg[idx] = new TLegend(0.75, 0.7, 0.99, 0.99);
	leg[idx]->AddEntry(p1,Form("#Delta A = %5.2f", delta_aic_2[idx]),"");
	leg[idx]->AddEntry(p1,Form("   #omega_{i} = %5.2f", weight_aic_2[idx]),"");
	leg[idx]->SetFillStyle(0);
	leg[idx]->SetBorderSize(0);
	leg[idx]->SetTextColor(kRed);
	}
	leg[0]->Draw();
		
	myC->cd(2);
	RooPlot * p2 = (RooPlot*)w[1]->obj("sideband_fit_singleExp_pdf_sExp_ext_frame");
	p2->SetTitle("single exponential");
	p2->GetXaxis()->SetTitleSize(0.09);
	p2->GetXaxis()->SetLabelSize(0.11);
	p2->GetYaxis()->SetTitleSize(0.11);
	p2->GetYaxis()->SetLabelSize(0.11);
	p2->GetYaxis()->SetTitleOffset(0.55);
	p2->GetYaxis()->SetNdivisions(505);
	p2->Draw();
	
	top_y = p1->GetMaximum();
	grshade[1] = new TGraph(2*(n_g+1));
        for ( int ig = 0; ig <= n_g; ig++ )
          {
            grshade[1]->SetPoint(ig, 120., (float)ig*(top_y/(float)n_g) );
            grshade[1]->SetPoint( (n_g+1) + ig, 135., top_y - (float)ig*(top_y/(float)n_g) );
          }
	grshade[1]->SetFillStyle(1001);
        grshade[1]->SetFillColor(kRed-9);
        grshade[1]->Draw("f same");
	leg[1]->Draw();
	myC->cd(3);
	RooPlot * p3 = (RooPlot*)w[3]->obj("sideband_fit_doublePow_dpow_ext_frame");
	p3->SetTitle("double power");
	p3->GetXaxis()->SetTitleSize(0.09);
	p3->GetXaxis()->SetLabelSize(0.11);
	p3->GetYaxis()->SetTitleSize(0.11);
	p3->GetYaxis()->SetLabelSize(0.11);
	p3->GetYaxis()->SetTitleOffset(0.55);
	p3->GetYaxis()->SetNdivisions(505);
	p3->Draw();
	//making blind band 
	top_y = p1->GetMaximum();
        grshade[2] = new TGraph(2*(n_g+1));
        for ( int ig = 0; ig <= n_g; ig++ )
          {
            grshade[2]->SetPoint(ig, 120., (float)ig*(top_y/(float)n_g) );
            grshade[2]->SetPoint( (n_g+1) + ig, 135., top_y - (float)ig*(top_y/(float)n_g) );
          }
        grshade[2]->SetFillStyle(1001);
        grshade[2]->SetFillColor(kRed-9);
        grshade[2]->Draw("f same");
	
	leg[2]->Draw();
	myC->cd(4);
	RooPlot * p4 = (RooPlot*)w[2]->obj("sideband_fit_singlePow_spow_ext_frame");
	p4->SetTitle("single power");
	p4->GetXaxis()->SetTitleSize(0.09);
	p4->GetXaxis()->SetLabelSize(0.11);
	p4->GetYaxis()->SetTitleSize(0.11);
	p4->GetYaxis()->SetLabelSize(0.11);
	p4->GetYaxis()->SetTitleOffset(0.55);
	p4->GetYaxis()->SetNdivisions(505);
	p4->Draw();
	
	top_y = p1->GetMaximum();
        grshade[3] = new TGraph(2*(n_g+1));
        for ( int ig = 0; ig <= n_g; ig++ )
          {
            grshade[3]->SetPoint(ig, 120., (float)ig*(top_y/(float)n_g) );
            grshade[3]->SetPoint( (n_g+1) + ig, 135., top_y - (float)ig*(top_y/(float)n_g) );
          }
        grshade[3]->SetFillStyle(1001);
        grshade[3]->SetFillColor(kRed-9);
        grshade[3]->Draw("f same");
	
	leg[3]->Draw();
	myC->cd(5);
	RooPlot * p5 = (RooPlot*)w[4]->obj("sideband_fit_poly2_pol2_ext_frame");
	p5->SetTitle("second order poly");
	p5->GetXaxis()->SetTitleSize(0.09);
	p5->GetXaxis()->SetLabelSize(0.11);
	p5->GetYaxis()->SetTitleSize(0.11);
	p5->GetYaxis()->SetLabelSize(0.11);
	p5->GetYaxis()->SetTitleOffset(0.55);
	p5->GetYaxis()->SetNdivisions(505);
	p5->Draw();
	
	top_y = p1->GetMaximum();
        grshade[4] = new TGraph(2*(n_g+1));
        for ( int ig = 0; ig <= n_g; ig++ )
          {
            grshade[4]->SetPoint(ig, 120., (float)ig*(top_y/(float)n_g) );
            grshade[4]->SetPoint( (n_g+1) + ig, 135., top_y - (float)ig*(top_y/(float)n_g) );
          }
        grshade[4]->SetFillStyle(1001);
        grshade[4]->SetFillColor(kRed-9);
        grshade[4]->Draw("f same");
	
	leg[4]->Draw();
	myC->cd(6);
	RooPlot * p6 = (RooPlot*)w[5]->obj("sideband_fit_poly3_pol3_ext_frame");
	p6->SetTitle("third order poly");
	p6->GetXaxis()->SetTitleSize(0.09);
	p6->GetXaxis()->SetLabelSize(0.11);
	p6->GetYaxis()->SetTitleSize(0.11);
	p6->GetYaxis()->SetLabelSize(0.11);
	p6->GetYaxis()->SetTitleOffset(0.55);
	p6->GetYaxis()->SetNdivisions(505);
	p6->Draw();

	top_y = p1->GetMaximum();
        grshade[5] = new TGraph(2*(n_g+1));
        for ( int ig = 0; ig <= n_g; ig++ )
          {
            grshade[5]->SetPoint(ig, 120., (float)ig*(top_y/(float)n_g) );
            grshade[5]->SetPoint( (n_g+1) + ig, 135., top_y - (float)ig*(top_y/(float)n_g) );
          }
        grshade[5]->SetFillStyle(1001);
        grshade[5]->SetFillColor(kRed-9);
        grshade[5]->Draw("f same");
	
	leg[5]->Draw();
	myC->cd(7);
	RooPlot * p7 = (RooPlot*)w[7]->obj("sideband_fit_poly4_pol4_ext_frame");
	p7->SetTitle("fourth order poly");
	p7->GetXaxis()->SetTitleSize(0.09);
	p7->GetXaxis()->SetLabelSize(0.11);
	p7->GetYaxis()->SetTitleSize(0.11);
	p7->GetYaxis()->SetLabelSize(0.11);
	p7->GetYaxis()->SetTitleOffset(0.55);
	p7->GetYaxis()->SetNdivisions(505);
	p7->Draw();

	top_y = p1->GetMaximum();
        grshade[6] = new TGraph(2*(n_g+1));
        for ( int ig = 0; ig <= n_g; ig++ )
          {
            grshade[6]->SetPoint(ig, 120., (float)ig*(top_y/(float)n_g) );
            grshade[6]->SetPoint( (n_g+1) + ig, 135., top_y - (float)ig*(top_y/(float)n_g) );
          }
        grshade[6]->SetFillStyle(1001);
        grshade[6]->SetFillColor(kRed-9);
        grshade[6]->Draw("f same");
	
	leg[6]->Draw();
	myC->cd(8);
	RooPlot * p8 = (RooPlot*)w[6]->obj("sideband_fit_modExp_mexp_ext_frame");
	p8->SetTitle("modified exponential");
	p8->GetXaxis()->SetTitleSize(0.09);
	p8->GetXaxis()->SetLabelSize(0.11);
	p8->GetYaxis()->SetTitleSize(0.11);
	p8->GetYaxis()->SetLabelSize(0.11);
	p8->GetYaxis()->SetTitleOffset(0.55);
	p8->GetYaxis()->SetNdivisions(505);
	p8->Draw();
	
	top_y = p1->GetMaximum();
        grshade[7] = new TGraph(2*(n_g+1));
        for ( int ig = 0; ig <= n_g; ig++ )
          {
            grshade[7]->SetPoint(ig, 120., (float)ig*(top_y/(float)n_g) );
            grshade[7]->SetPoint( (n_g+1) + ig, 135., top_y - (float)ig*(top_y/(float)n_g) );
          }
        grshade[7]->SetFillStyle(1001);
        grshade[7]->SetFillColor(kRed-9);
        grshade[7]->Draw("f same");
	
	leg[7]->Draw();
	
/*	myC->cd(8);
	TLegend* leg8 = new TLegend(0.2, 0.25, 0.8, 0.75);
	leg8->AddEntry(p1,Form("M_{R} cut: %s", MRcut.c_str()),"");
	leg8->AddEntry(p1,Form(" R^{2} cut: %s", RSQcut.c_str()),"");
	leg8->SetFillStyle(0);
	leg8->SetBorderSize(0);
	leg8->SetTextColor(kRed);
	leg8->Draw();
*/
	std::string str = "AIC_output/"+MRcut+"_"+RSQcut+".pdf";
	const char * file_Name = str.c_str();
	myC->SaveAs(file_Name);

};

void PlotSidebandFit(std::string MRcut,std::string RSQcut, RooWorkspace *w[])
{
	gStyle->SetTitleSize(0.16,"t");
	TCanvas *myC= new TCanvas("c1","c1",200,10,800,600);
	myC->cd();
    gPad->SetLeftMargin(0.12);
	gPad->SetBottomMargin(0.18);

	myC->Divide(2,4,0);
	myC->cd(1);
	RooPlot * p1 = (RooPlot*)w[0]->obj("sideband_fit_doubleExp_pdf_dExp_N1N2_frame");
	p1->SetTitle("double exponential");
	p1->GetXaxis()->SetTitleSize(0.09);
    p1->GetXaxis()->SetLabelSize(0.11);
	p1->GetYaxis()->SetTitleSize(0.11);
	p1->GetYaxis()->SetLabelSize(0.11);
	p1->GetYaxis()->SetTitleOffset(0.55);
	p1->GetYaxis()->SetNdivisions(505);
p1->Draw();
	myC->cd(2);
	RooPlot * p2 = (RooPlot*)w[1]->obj("sideband_fit_singleExp_pdf_sExp_ext_frame");
	p2->SetTitle("single exponential");
	p2->GetXaxis()->SetTitleSize(0.09);
    p2->GetXaxis()->SetLabelSize(0.11);
	p2->GetYaxis()->SetTitleSize(0.11);
	p2->GetYaxis()->SetLabelSize(0.11);
	p2->GetYaxis()->SetTitleOffset(0.55);
	p2->GetYaxis()->SetNdivisions(505);
p2->Draw();
	myC->cd(3);
	RooPlot * p3 = (RooPlot*)w[2]->obj("sideband_fit_singlePow_spow_ext_frame");
	p3->SetTitle("single power");
	p3->GetXaxis()->SetTitleSize(0.09);
    p3->GetXaxis()->SetLabelSize(0.11);
	p3->GetYaxis()->SetTitleSize(0.11);
	p3->GetYaxis()->SetLabelSize(0.11);
	p3->GetYaxis()->SetTitleOffset(0.55);
	p3->GetYaxis()->SetNdivisions(505);
p3->Draw();
	myC->cd(4);
	RooPlot * p4 = (RooPlot*)w[3]->obj("sideband_fit_doublePow_dpow_ext_frame");
	p4->SetTitle("double power");
	p4->GetXaxis()->SetTitleSize(0.09);
    p4->GetXaxis()->SetLabelSize(0.11);
	p4->GetYaxis()->SetTitleSize(0.11);
	p4->GetYaxis()->SetLabelSize(0.11);
	p4->GetYaxis()->SetTitleOffset(0.55);
	p4->GetYaxis()->SetNdivisions(505);
p4->Draw();
	myC->cd(5);
	RooPlot * p5 = (RooPlot*)w[4]->obj("sideband_fit_poly2_pol2_ext_frame");
	p5->SetTitle("second order poly");
	p5->GetXaxis()->SetTitleSize(0.09);
    p5->GetXaxis()->SetLabelSize(0.11);
	p5->GetYaxis()->SetTitleSize(0.11);
	p5->GetYaxis()->SetLabelSize(0.11);
	p5->GetYaxis()->SetTitleOffset(0.55);
	p5->GetYaxis()->SetNdivisions(505);
p5->Draw();

//	gStyle->SetTitleSize(0.14,"t");
	myC->cd(6);
	RooPlot * p6 = (RooPlot*)w[5]->obj("sideband_fit_poly3_pol3_ext_frame");
	p6->SetTitle("third order poly");
	p6->SetTitleSize(0.1);
	p6->GetXaxis()->SetTitleSize(0.09);
    p6->GetXaxis()->SetLabelSize(0.11);
	p6->GetYaxis()->SetTitleSize(0.11);
	p6->GetYaxis()->SetLabelSize(0.11);
	p6->GetYaxis()->SetTitleOffset(0.55);
	p6->GetYaxis()->SetNdivisions(505);
p6->Draw();
	myC->cd(7);
	RooPlot * p7 = (RooPlot*)w[6]->obj("sideband_fit_modExp_mexp_ext_frame");
	p7->SetTitle("modified exponential");
	p7->SetTitleSize(0.1);
	p7->GetXaxis()->SetTitleSize(0.09);
    p7->GetXaxis()->SetLabelSize(0.11);
	p7->GetYaxis()->SetTitleSize(0.11);
	p7->GetYaxis()->SetLabelSize(0.11);
	p7->GetYaxis()->SetTitleOffset(0.55);
	p7->GetYaxis()->SetNdivisions(505);
p7->Draw();
	std::string str = "AIC_output/"+MRcut+"_"+RSQcut+".pdf";
	const char * file_Name = str.c_str();
	myC->SaveAs(file_Name);
};
