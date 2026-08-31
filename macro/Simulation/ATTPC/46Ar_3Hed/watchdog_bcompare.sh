#!/usr/bin/env bash
# Guard the overnight run: check the FIRST fitted sample and kill the campaign if it produced no
# fitted tracks. Twice tonight a par that parsed cleanly yielded zero fits (SamplingRate 5 MHz ->
# TBTime -1 ns), each time after ~85 min. This turns that into a 1-sample loss instead of a 9-hour
# one. A separate process on purpose: never edit a script that is already running.
set -eo pipefail
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
ARMS="/mnt/f/ar46_3hed_mb_B285 /mnt/f/ar46_3hed_mb_B38"
STATUS=/mnt/f/AR46_OVERNIGHT_STATUS.txt
say() { echo "[$(date +%F' '%H:%M:%S)] WATCHDOG: $*" | tee -a "$STATUS"; }

say "armed on the 2.85 and 3.8 T arms (2.0 T already validated by hand: 11807 fitted tracks)"
for OUT in $ARMS; do
  F=$OUT/gs_s3001_genfitter_d.root
  say "waiting for $(basename "$OUT") first sample"
  while ! { [ -s "$F" ] && grep -q "Done" "$OUT/gs_s3001_fit.log" 2>/dev/null; }; do
    if ! pgrep -f "[r]un_3Hed_bcompare" >/dev/null; then
      say "campaign ended before $(basename "$OUT") produced a sample -- stopping watch"; exit 0
    fi
    sleep 120
  done
  n=$(cd "$REPO" && export VMCWORKDIR=$PWD && set +u && source build/config.sh >/dev/null 2>&1 && set -u
      root -b -q -l -e "
        gSystem->Load(\"libAtReconstruction.so\");
        TFile*f=TFile::Open(\"$F\"); TTree*t=(TTree*)f->Get(\"cbmsim\");
        TClonesArray*te=nullptr; t->SetBranchAddress(\"AtTrackingEvent\",&te);
        long n=0; for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i);
          if(!te||!te->GetEntriesFast())continue; auto*ev=(AtTrackingEvent*)te->At(0); if(!ev)continue;
          for(auto&ft:ev->GetFittedTracks()) if(ft && ft->GetKinematics().kineticEnergy>0) ++n;}
        printf(\"NFIT %ld\n\",n);" 2>/dev/null | awk '/^NFIT /{print $2}')
  n=${n:-0}
  if [ "$n" -lt 500 ]; then
    say "*** $(basename "$OUT") FIRST SAMPLE HAS ONLY $n FITTED TRACKS -- KILLING ***"
    say "    check TBTime in its reco log; -1 ns means the SamplingRate is not a GET rate"
    pkill -f "[r]un_3Hed_bcompare" 2>/dev/null || true
    pkill -f "[a]ccumulate_3Hed" 2>/dev/null || true
    sleep 2
    pkill -f "[r]un_reco_Ar46_TC" 2>/dev/null || true
    exit 1
  fi
  say "$(basename "$OUT") OK: $n fitted tracks"
done
say "all arms validated; campaign left to finish"
