void plot(){
gStyle->SetTitleYOffset(0.45);
gStyle->SetTitleXOffset(0.45);
gStyle->SetTitleYSize(0.1);
gStyle->SetTitleXSize(0.1);
gStyle->SetLabelSize(0.045,"Y");
gStyle->SetLabelSize(0.045,"X");
gStyle->SetOptStat(0);
TCanvas *c1=new TCanvas("Raw scaler data","Raw scaler data",1000,1000);
TPad *pad11=new TPad("pad11","pad11",0.05,0.92,0.97,0.97);
pad11->Draw();
pad11->SetFillColor(kGray);
pad11->cd();
TText *ttext= new TText(0.5,0.5,"Anguler Distribution  ");
ttext->SetTextAlign(22);
ttext->SetTextSize(0.5);
ttext->Draw();
c1->cd();
TPad *pad22=new TPad("pad22","pad22",0.05,0.05,0.97,0.90);
pad22->Draw();
pad22->cd();
pad22->Divide(5,5);
Double_t chi2,s2,a,v,w,rv,av,k;
ifstream input_file;
input_file.open("elastic_data.txt");
Int_t np =10;
const int NbinsX= 50;
const int NbinsY= 50;
Int_t Xbins2,Ybins2,Xbina,Ybina,Xbinv,Ybinv,Xbinw,Ybinw,Xbinrv,Ybinrv,Xbinav,Ybinav;
Double_t s2MIN = 1;
Double_t s2MAX = 1.9;
Double_t aMIN = 11.3;
Double_t aMAX = 12.4;
Double_t vMIN = 140;
Double_t vMAX = 220;
Double_t wMIN = 5;
Double_t wMAX = 54;
Double_t rvMIN = 1.2;
Double_t rvMAX = 1.6;
Double_t avMIN = 0.4;
Double_t avMAX = 0.7;

TH2F * hMinChi2rvav=new TH2F("hMinChi2rvav","",NbinsX,avMIN,avMAX,NbinsY,rvMIN,rvMAX);
TH2F * hMinChi2s2a=new TH2F("hMinChi2s2a","",NbinsX,aMIN,aMAX,NbinsY,s2MIN,s2MAX);
TH2F * hMinChi2s2v=new TH2F("hMinChi2s2v","",NbinsX,vMIN,vMAX,NbinsY,s2MIN,s2MAX);
TH2F * hMinChi2s2w=new TH2F("hMinChi2s2w","",NbinsX,wMIN,wMAX,NbinsY,s2MIN,s2MAX);
TH2F * hMinChi2s2rv=new TH2F("hMinChi2s2rv","",NbinsX,rvMIN,rvMAX,NbinsY,s2MIN,s2MAX);
TH2F * hMinChi2s2av=new TH2F("hMinChi2s2av","",NbinsX,avMIN,avMAX,NbinsY,s2MIN,s2MAX);
TH2F * hMinChi2av=new TH2F("hMinChi2av","",NbinsX,vMIN,vMAX,NbinsY,aMIN,aMAX);
TH2F * hMinChi2aw=new TH2F("hMinChi2aw","",NbinsX,wMIN,wMAX,NbinsY,aMIN,aMAX);
TH2F * hMinChi2arv=new TH2F("hMinChi2arv","",NbinsX,rvMIN,rvMAX,NbinsY,aMIN,aMAX);
TH2F * hMinChi2aav=new TH2F("hMinChi2aav","",NbinsX,avMIN,avMAX,NbinsY,aMIN,aMAX);
TH2F * hMinChi2vw=new TH2F("hMinChi2vw","",NbinsX,wMIN,wMAX,NbinsY,vMIN,vMAX);
TH2F * hMinChi2vrv=new TH2F("hMinChi2vrv","",NbinsX,rvMIN,rvMAX,NbinsY,vMIN,vMAX);
TH2F * hMinChi2vav=new TH2F("hMinChi2vav","",NbinsX,avMIN,avMAX,NbinsY,vMIN,vMAX);
TH2F * hMinChi2wrv=new TH2F("hMinChi2wrv","",NbinsX,rvMIN,rvMAX,NbinsY,wMIN,wMAX);
TH2F * hMinChi2wav=new TH2F("hMinChi2wav","",NbinsX,avMIN,avMAX,NbinsY,wMIN,wMAX);
Double_t MinChi2rvav[NbinsX][NbinsY];
Double_t MinChi2s2a[NbinsX][NbinsY];
Double_t MinChi2s2v[NbinsX][NbinsY];
Double_t MinChi2s2w[NbinsX][NbinsY];
Double_t MinChi2s2rv[NbinsX][NbinsY];
Double_t MinChi2s2av[NbinsX][NbinsY];
Double_t MinChi2av[NbinsX][NbinsY];
Double_t MinChi2aw[NbinsX][NbinsY];
Double_t MinChi2arv[NbinsX][NbinsY];
Double_t MinChi2aav[NbinsX][NbinsY];
Double_t MinChi2vw[NbinsX][NbinsY];
Double_t MinChi2vrv[NbinsX][NbinsY];
Double_t MinChi2vav[NbinsX][NbinsY];
Double_t MinChi2wrv[NbinsX][NbinsY];
Double_t MinChi2wav[NbinsX][NbinsY];
Double_t s2ForBinX[NbinsX];
Double_t s2ForBinY[NbinsY];
Double_t aForBinX[NbinsX];
Double_t aForBinY[NbinsY];
Double_t vForBinX[NbinsX];
Double_t vForBinY[NbinsY];
Double_t wForBinX[NbinsX];
Double_t wForBinY[NbinsY];
Double_t rvForBinX[NbinsX];
Double_t rvForBinY[NbinsY];
Double_t avForBinX[NbinsX];
Double_t avForBinY[NbinsY];

    
for(int m=0; m < NbinsX ; m++){
for (int n=0; n < NbinsY ; n++){
MinChi2rvav[m][n] = 999999999.;
MinChi2s2a[m][n] = 999999999.;
MinChi2s2v[m][n] = 999999999.;
MinChi2s2w[m][n] = 999999999.;
MinChi2s2rv[m][n] = 999999999.;
MinChi2s2av[m][n] = 999999999.;
MinChi2av[m][n] = 999999999.;
MinChi2aw[m][n] = 999999999.;
MinChi2arv[m][n] = 999999999.;
MinChi2aav[m][n] = 999999999.;
MinChi2vw[m][n] = 999999999.;
MinChi2vrv[m][n] = 999999999.;
MinChi2vav[m][n] = 999999999.;
MinChi2wrv[m][n] = 999999999.;
MinChi2wav[m][n] = 999999999.;
}
}
    
    
while(input_file.good()){
input_file >>k>>chi2>>s2>>a>>v>>w>>rv>>av;
Xbins2 = (int) ceil( NbinsX * (s2 - s2MIN) / (s2MAX - s2MIN) ) -1;
Ybins2 = (int) ceil( NbinsY * (s2 - s2MIN) / (s2MAX - s2MIN) ) -1;
Xbina = (int) ceil( NbinsX * (a - aMIN) / (aMAX - aMIN) ) -1;
Ybina = (int) ceil( NbinsY * (a - aMIN) / (aMAX - aMIN) ) -1;
Xbinv = (int) ceil( NbinsX * (v - vMIN) / (vMAX - vMIN) ) -1;
Ybinv = (int) ceil( NbinsY * (v - vMIN) / (vMAX - vMIN) ) -1;
Xbinw = (int) ceil( NbinsX * (w - wMIN) / (wMAX - wMIN) ) -1;
Ybinw = (int) ceil( NbinsY * (w - wMIN) / (wMAX - wMIN) ) -1;
Xbinrv = (int) ceil( NbinsX * (rv - rvMIN) / (rvMAX - rvMIN) ) -1;
Ybinrv = (int) ceil( NbinsY * (rv - rvMIN) / (rvMAX - rvMIN) ) -1;
Xbinav = (int) ceil( NbinsX * (av - avMIN) / (avMAX - avMIN) ) -1;
Ybinav = (int) ceil( NbinsY * (av - avMIN) / (avMAX - avMIN) ) -1;
s2ForBinX[Xbins2] = s2;
s2ForBinY[Ybins2] = s2;
aForBinX[Xbina] = a;
aForBinY[Ybina] = a;
vForBinX[Xbinv] = v;
vForBinY[Ybinv] = v;
wForBinX[Xbinw] = w;
wForBinY[Ybinw] = w;
rvForBinX[Xbinrv] = rv;
rvForBinY[Ybinrv] = rv;
avForBinX[Xbinav] = av;
avForBinY[Ybinav] = av;
///check if the current chi-square is less than any others that fell in this bin
if(chi2 < MinChi2rvav[Xbinav][Ybinrv]){
MinChi2rvav[Xbinav][Ybinrv] = chi2;
}
///check if the current chi-square is less than any others that fell in this bin
if(chi2 < MinChi2s2a[Xbina][Ybins2]){
MinChi2s2a[Xbina][Ybins2] = chi2;
}
///check if the current chi-square is less than any others that fell in this bin
if(chi2 < MinChi2s2v[Xbinv][Ybins2]){
MinChi2s2v[Xbinv][Ybins2] = chi2;
}
///check if the current chi-square is less than any others that fell in this bin
if(chi2 < MinChi2s2w[Xbinw][Ybins2]){
MinChi2s2w[Xbinw][Ybins2] = chi2;
}
///check if the current chi-square is less than any others that fell in this bin
if(chi2 < MinChi2s2rv[Xbinrv][Ybins2]){
MinChi2s2rv[Xbinrv][Ybins2] = chi2;
}
///check if the current chi-square is less than any others that fell in this bin
if(chi2 < MinChi2s2av[Xbinav][Ybins2]){
MinChi2s2av[Xbinav][Ybins2] = chi2;
}
///check if the current chi-square is less than any others that fell in this bin
if(chi2 < MinChi2av[Xbinv][Ybina]){
MinChi2av[Xbinv][Ybina] = chi2;
}
///check if the current chi-square is less than any others that fell in this bin
if(chi2 < MinChi2aw[Xbinw][Ybina]){
MinChi2aw[Xbinw][Ybina] = chi2;
}
///check if the current chi-square is less than any others that fell in this bin
if(chi2 < MinChi2arv[Xbinrv][Ybina]){
MinChi2arv[Xbinrv][Ybina] = chi2;
}
///check if the current chi-square is less than any others that fell in this bin
if(chi2 < MinChi2aav[Xbinav][Ybina]){
MinChi2aav[Xbinav][Ybina] = chi2;
}
///check if the current chi-square is less than any others that fell in this bin
if(chi2 < MinChi2vw[Xbinw][Ybinv]){
MinChi2vw[Xbinw][Ybinv] = chi2;
}
///check if the current chi-square is less than any others that fell in this bin
if(chi2 < MinChi2vrv[Xbinrv][Ybinv]){
MinChi2vrv[Xbinrv][Ybinv] = chi2;
}
///check if the current chi-square is less than any others that fell in this bin
if(chi2 < MinChi2vav[Xbinav][Ybinv]){
MinChi2vav[Xbinav][Ybinv] = chi2;
}
///check if the current chi-square is less than any others that fell in this bin
if(chi2 < MinChi2wrv[Xbinrv][Ybinw]){
MinChi2wrv[Xbinrv][Ybinw] = chi2;
}
///check if the current chi-square is less than any others that fell in this bin
if(chi2 < MinChi2wav[Xbinav][Ybinw]){
MinChi2wav[Xbinav][Ybinw] = chi2;
}
    
}//end looping through the data file
    
for(int m=0 ; m < NbinsX ; m++){
for(int n=0 ; n < NbinsY ; n++){
hMinChi2rvav->Fill(avForBinX[m],rvForBinY[n],MinChi2rvav[m][n]);
hMinChi2s2a->Fill(aForBinX[m],s2ForBinY[n],MinChi2s2a[m][n]);
hMinChi2s2v->Fill(vForBinX[m],s2ForBinY[n],MinChi2s2v[m][n]);
hMinChi2s2w->Fill(wForBinX[m],s2ForBinY[n],MinChi2s2w[m][n]);
hMinChi2s2rv->Fill(rvForBinX[m],s2ForBinY[n],MinChi2s2rv[m][n]);
hMinChi2s2av->Fill(avForBinX[m],s2ForBinY[n],MinChi2s2av[m][n]);
hMinChi2av->Fill(vForBinX[m],aForBinY[n],MinChi2av[m][n]);
hMinChi2aw->Fill(wForBinX[m],aForBinY[n],MinChi2aw[m][n]);
hMinChi2arv->Fill(rvForBinX[m],aForBinY[n],MinChi2arv[m][n]);
hMinChi2aav->Fill(avForBinX[m],aForBinY[n],MinChi2aav[m][n]);
hMinChi2vw->Fill(wForBinX[m],vForBinY[n],MinChi2vw[m][n]);
hMinChi2vrv->Fill(rvForBinX[m],vForBinY[n],MinChi2vrv[m][n]);
hMinChi2vav->Fill(avForBinX[m],vForBinY[n],MinChi2vav[m][n]);
hMinChi2wrv->Fill(rvForBinX[m],wForBinY[n],MinChi2wrv[m][n]);
hMinChi2wav->Fill(avForBinX[m],wForBinY[n],MinChi2wav[m][n]);
cout <<n<< "  " << m<<"   " << avForBinX[m]<<"   "  << s2ForBinY[n]<<  "    "<<m <<  "     "<<MinChi2s2av[m][n]<<endl;
}
}
    
pad22->cd(1);
c1->cd(1);
hMinChi2rvav->Smooth();
hMinChi2rvav->SetMaximum(500);
hMinChi2rvav->GetXaxis()->SetTitle("av");
hMinChi2rvav->GetYaxis()->SetTitle("rv");
hMinChi2rvav->DrawCopy("colz");
Double_t ang_contrvav = hMinChi2rvav->GetMinimum();
cout<< ang_contrvav << endl;
double contours1[2];
contours1[0] = ang_contrvav + 22.44;
contours1[1] = ang_contrvav + 28.41;
hMinChi2rvav->SetContour(2,contours1);
hMinChi2rvav->Draw("cont3 same");
hMinChi2rvav->SetLineColor(kRed);
pad22->cd(21);
c1->cd(21);
hMinChi2s2a->Smooth();
//hMinChi2s2a->SetMaximum(1000);
hMinChi2s2a->GetXaxis()->SetTitle("a");
hMinChi2s2a->GetYaxis()->SetTitle("s2");
hMinChi2s2a->DrawCopy("colz");
Double_t ang_conts2a = hMinChi2s2a->GetMinimum();
cout<< ang_conts2a << endl;
double contours2[2];
contours2[0] = ang_conts2a + 22.44;
contours2[1] = ang_conts2a + 28.41;
hMinChi2s2a->SetContour(2,contours2);
hMinChi2s2a->Draw("cont3 same");
hMinChi2s2a->SetLineColor(kRed);
pad22->cd(22);
c1->cd(22);
hMinChi2s2v->Smooth();
//hMinChi2s2a->SetMaximum(1000);
hMinChi2s2v->GetXaxis()->SetTitle("v");
hMinChi2s2v->GetYaxis()->SetTitle("s2");
hMinChi2s2v->DrawCopy("colz");
Double_t ang_conts2v = hMinChi2s2v->GetMinimum();
cout<< ang_conts2v << endl;
double contours3[2];
contours3[0] = ang_conts2v + 22.44;
contours3[1] = ang_conts2v + 28.41;
hMinChi2s2v->SetContour(2,contours3);
hMinChi2s2v->Draw("cont3 same");
hMinChi2s2v->SetLineColor(kRed);
pad22->cd(23);
c1->cd(23);
hMinChi2s2w->Smooth();
//hMinChi2s2a->SetMaximum(1000);
hMinChi2s2w->GetXaxis()->SetTitle("w");
hMinChi2s2w->GetYaxis()->SetTitle("s2");
hMinChi2s2w->DrawCopy("colz");
Double_t ang_conts2w = hMinChi2s2w->GetMinimum();
cout<< ang_conts2w << endl;
double contours4[2];
contours4[0] = ang_conts2w + 22.44;
contours4[1] = ang_conts2w + 28.41;
hMinChi2s2w->SetContour(2,contours4);
hMinChi2s2w->Draw("cont3 same");
hMinChi2s2w->SetLineColor(kRed);
pad22->cd(24);
c1->cd(24);
hMinChi2s2rv->Smooth();
hMinChi2s2rv->SetMaximum(500);
hMinChi2s2rv->GetXaxis()->SetTitle("rv");
hMinChi2s2rv->GetYaxis()->SetTitle("s2");
hMinChi2s2rv->DrawCopy("colz");
Double_t ang_conts2rv = hMinChi2s2rv->GetMinimum();
cout<< ang_conts2rv << endl;
double contours5[2];
contours5[0] = ang_conts2rv + 22.44;
contours5[1] = ang_conts2rv + 28.41;
hMinChi2s2rv->SetContour(2,contours5);
hMinChi2s2rv->Draw("cont3 same");
hMinChi2s2rv->SetLineColor(kRed);
pad22->cd(25);
c1->cd(25);
hMinChi2s2av->Smooth();
hMinChi2s2av->SetMaximum(500);
hMinChi2s2av->GetXaxis()->SetTitle("av");
hMinChi2s2av->GetYaxis()->SetTitle("s2");
hMinChi2s2av->DrawCopy("colz");
Double_t ang_conts2av = hMinChi2s2av->GetMinimum();
cout<< ang_conts2av << endl;
double contours6[2];
contours6[0] = ang_conts2av + 22.44;
contours6[1] = ang_conts2av + 28.41;
hMinChi2s2av->SetContour(2,contours6);
hMinChi2s2av->Draw("cont3 same");
hMinChi2s2av->SetLineColor(kRed);
pad22->cd(16);
c1->cd(16);
hMinChi2av->Smooth();
// hMinChi2av->SetMaximum(500);
hMinChi2av->GetXaxis()->SetTitle("v");
hMinChi2av->GetYaxis()->SetTitle("a");
hMinChi2av->DrawCopy("colz");
Double_t ang_contav = hMinChi2av->GetMinimum();
cout<< ang_contav << endl;
double contours7[2];
contours7[0] = ang_contav + 22.44;
contours7[1] = ang_contav + 28.41;
hMinChi2av->SetContour(2,contours7);
hMinChi2av->Draw("cont3 same");
hMinChi2av->SetLineColor(kRed);
pad22->cd(17);
c1->cd(17);
hMinChi2aw->Smooth();
//hMinChi2aw->SetMaximum(500);
hMinChi2aw->GetXaxis()->SetTitle("w");
hMinChi2aw->GetYaxis()->SetTitle("a");
hMinChi2aw->DrawCopy("colz");
Double_t ang_contaw = hMinChi2aw->GetMinimum();
cout<< ang_contaw << endl;
double contours8[2];
contours8[0] = ang_contaw + 22.44;
contours8[1] = ang_contaw + 28.41;
hMinChi2aw->SetContour(2,contours8);
hMinChi2aw->Draw("cont3 same");
hMinChi2aw->SetLineColor(kRed);
pad22->cd(18);
c1->cd(18);
hMinChi2arv->Smooth();
hMinChi2arv->SetMaximum(500);
hMinChi2arv->GetXaxis()->SetTitle("rv");
hMinChi2arv->GetYaxis()->SetTitle("a");
hMinChi2arv->DrawCopy("colz");
Double_t ang_contarv = hMinChi2arv->GetMinimum();
cout<< ang_contarv << endl;
double contours9[2];
contours9[0] = ang_contarv + 22.44; //22.44;
contours9[1] = ang_contarv + 28.41; //22.44;
hMinChi2arv->SetContour(2,contours9);
hMinChi2arv->Draw("cont3 same");
hMinChi2arv->SetLineColor(kRed);
pad22->cd(19);
c1->cd(19);
hMinChi2aav->Smooth();
hMinChi2aav->SetMaximum(500);
hMinChi2aav->GetXaxis()->SetTitle("av");
hMinChi2aav->GetYaxis()->SetTitle("a");
hMinChi2aav->DrawCopy("colz");
Double_t ang_contaav = hMinChi2aav->GetMinimum();
cout<< ang_contaav << endl;
double contours10[2];
contours10[0] = ang_contaav + 22.44; //22.44;
contours10[1] = ang_contaav + 28.41; //22.44;
hMinChi2aav->SetContour(2,contours10);
hMinChi2aav->Draw("cont3 same");
hMinChi2aav->SetLineColor(kRed);
pad22->cd(11);
c1->cd(11);
hMinChi2vw->Smooth();
//hMinChi2vw->SetMaximum(400);
hMinChi2vw->GetXaxis()->SetTitle("w");
hMinChi2vw->GetYaxis()->SetTitle("v");
hMinChi2vw->DrawCopy("colz");
Double_t ang_contvw = hMinChi2vw->GetMinimum();
cout<< ang_contvw << endl;
double contours11[2];
contours11[0] = ang_contvw + 22.44;
contours11[1] = ang_contvw + 28.41;
hMinChi2vw->SetContour(2,contours11);
hMinChi2vw->Draw("cont3 same");
hMinChi2vw->SetLineColor(kRed);
pad22->cd(12);
c1->cd(12);
hMinChi2vrv->Smooth();
hMinChi2vrv->SetMaximum(500);
hMinChi2vrv->GetXaxis()->SetTitle("rv");
hMinChi2vrv->GetYaxis()->SetTitle("v");
hMinChi2vrv->DrawCopy("colz");
Double_t ang_contvrv = hMinChi2vrv->GetMinimum();
cout<< ang_contvrv << endl;
double contours12[2];
contours12[0] = ang_contvrv + 22.44;
contours12[1] = ang_contvrv + 28.41;
hMinChi2vrv->SetContour(2,contours12);
hMinChi2vrv->Draw("cont3 same");
hMinChi2vrv->SetLineColor(kRed);
pad22->cd(13);
c1->cd(13);
hMinChi2vav->Smooth();
hMinChi2vav->SetMaximum(500);
hMinChi2vav->GetXaxis()->SetTitle("av");
hMinChi2vav->GetYaxis()->SetTitle("v");
hMinChi2vav->DrawCopy("colz");
Double_t ang_contvav = hMinChi2vav->GetMinimum();
cout<< ang_contvav << endl;
double contours13[2];
contours13[0] = ang_contvav + 22.44;
contours13[1] = ang_contvav + 28.41;
hMinChi2vav->SetContour(2,contours13);
hMinChi2vav->Draw("cont3 same");
hMinChi2vav->SetLineColor(kRed);
pad22->cd(6);
c1->cd(6);
hMinChi2wrv->Smooth();
hMinChi2wrv->SetMaximum(500);
hMinChi2wrv->GetXaxis()->SetTitle("rv");
hMinChi2wrv->GetYaxis()->SetTitle("w");
hMinChi2wrv->DrawCopy("colz");
Double_t ang_contwrv = hMinChi2wrv->GetMinimum();
cout<< ang_contwrv << endl;
double contours14[2];
contours14[0] = ang_contwrv + 22.44;
contours14[1] = ang_contwrv + 28.41;
hMinChi2wrv->SetContour(2,contours14);
hMinChi2wrv->Draw("cont3 same");
hMinChi2wrv->SetLineColor(kRed);
pad22->cd(7);
c1->cd(7);
hMinChi2wav->Smooth();
hMinChi2wav->SetMaximum(500);
hMinChi2wav->GetXaxis()->SetTitle("av");
hMinChi2wav->GetYaxis()->SetTitle("w");
hMinChi2wav->DrawCopy("colz");
Double_t ang_contwav = hMinChi2wav->GetMinimum();
cout<< ang_contwav << endl;
double contours15[2];
contours15[0] = ang_contwav + 22.44;
contours15[1] = ang_contwav + 28.41;
hMinChi2wav->SetContour(2,contours15);
hMinChi2wav->Draw("cont3 same");
hMinChi2wav->SetLineColor(kRed);
    
}







