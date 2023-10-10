void test_OTPCMap()
{

   gSystem->Load("libAtTpcMap.so");
   AtMap *otpc = new AtOTPCMap();
   otpc->GeneratePadPlane();

   TH2Poly *padplane = otpc->GetPadPlane();

   for (auto i = 0; i < 20; ++i) {

      Int_t bin = padplane->Fill(i + 0.5, -i - 0.5, i * 100);
      std::cout << " Bin : " << bin << "\n";
   }

   TArc *arc = new TArc(0.0, 0.0, 25.0);
   arc->SetLineColor(kRed);
   arc->SetLineWidth(4);
   arc->SetFillStyle(0);

   padplane->Draw("ZCOL L");
   arc->Draw();
}
