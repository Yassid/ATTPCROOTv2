#!/usr/bin/env python
"""Self-contained HTML viewer to navigate events and compare ATTPCROOT vs Spyral point clouds.
ATTPCROOT side = AtPSAMultiFit aligned to Spyral's peak selection (threshold 40, prominence 20, sep 50).
ATTPCROOT z is flipped to Spyral's convention (z -> 1150 - z) so the same event aligns spatially.
In: mfs_cloud.csv (ATTPCROOT), spyral_cloud.csv (Spyral).  Out: compare_psa_viewer.html
Run: ~/gnn_env/bin/python build_psa_viewer.py
"""
import pandas as pd, json
ZFLIP = 1150.0
a = pd.read_csv("mfs_cloud.csv")        # ATTPCROOT MultiFit, Spyral-aligned
s = pd.read_csv("spyral_cloud.csv")     # Spyral
a["z"] = ZFLIP - a["z"]                  # to Spyral z convention
events = []
for eid in sorted(set(a.eid) & set(s.eid)):
    ga, gs = a[a.eid == eid], s[s.eid == eid]
    events.append({
        "id": int(eid),
        "ax": [round(float(v), 1) for v in ga.x], "ay": [round(float(v), 1) for v in ga.y],
        "az": [round(float(v), 1) for v in ga.z],
        "sx": [round(float(v), 1) for v in gs.x], "sy": [round(float(v), 1) for v in gs.y],
        "sz": [round(float(v), 1) for v in gs.z],
    })
DATA = json.dumps(events)

HTML = """<!DOCTYPE html><html><head><meta charset="utf-8">
<title>ATTPCROOT vs Spyral PSA — run_0300</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
 body{font-family:sans-serif;margin:0;background:#111;color:#eee}
 #bar{padding:8px 12px;background:#1c1c1c;position:sticky;top:0;display:flex;gap:10px;align-items:center;flex-wrap:wrap}
 button{background:#333;color:#eee;border:1px solid #555;padding:6px 12px;border-radius:5px;cursor:pointer;font-size:14px}
 button:hover{background:#444}
 #plots{display:flex;gap:6px}#at,#sp{width:49vw;height:82vh}
 input{width:70px;background:#222;color:#eee;border:1px solid #555;padding:5px;border-radius:4px}
 #cnt{color:#9cf} .lab{font-weight:bold}
</style></head><body>
<div id="bar">
 <button onclick="step(-1)">&#9664; Prev</button>
 <button onclick="step(1)">Next &#9654;</button>
 event <input id="jump" type="number" onchange="goto(+this.value)">
 <span id="cnt"></span>
 <button onclick="setProj('zy')">z-y</button>
 <button onclick="setProj('xy')">x-y</button>
 <span class="lab" style="margin-left:auto;color:#8cf">left: ATTPCROOT (Spyral-aligned: thr40, prom20)</span>
 <span class="lab" style="color:#fc8">right: Spyral</span>
</div>
<div id="plots"><div id="at"></div><div id="sp"></div></div>
<script>
const DATA = __DATA__;
let idx=0, proj='zy';
function trace(ev,pre,color){
 const hx = proj=='zy'? ev[pre+'z'] : ev[pre+'x'];
 return [{x:hx,y:ev[pre+'y'],mode:'markers',type:'scattergl',marker:{size:4,color:color}}];
}
function layout(title){return{title:{text:title,font:{color:'#eee'}},paper_bgcolor:'#111',plot_bgcolor:'#181818',
 font:{color:'#bbb'},xaxis:{title:proj=='zy'?'z (mm)':'x (mm)',gridcolor:'#333',range:proj=='zy'?[50,1100]:[-280,280]},
 yaxis:{title:'y (mm)',gridcolor:'#333',range:[-280,280]},margin:{l:50,r:10,t:40,b:45},showlegend:false};}
function render(){
 const ev=DATA[idx];
 document.getElementById('cnt').textContent=`#${idx+1}/${DATA.length}  (event ${ev.id})  ATTPCROOT:${ev.ay.length}  Spyral:${ev.sy.length} hits`;
 document.getElementById('jump').value=ev.id;
 Plotly.react('at',trace(ev,'a','#5ab4ff'),layout('ATTPCROOT — event '+ev.id),{responsive:true});
 Plotly.react('sp',trace(ev,'s','#ffae5a'),layout('SPYRAL — event '+ev.id),{responsive:true});
}
function step(d){idx=(idx+d+DATA.length)%DATA.length;render();}
function goto(id){const i=DATA.findIndex(e=>e.id==id);if(i>=0){idx=i;render();}}
function setProj(p){proj=p;render();}
document.addEventListener('keydown',e=>{if(e.key=='ArrowLeft')step(-1);if(e.key=='ArrowRight')step(1);});
render();
</script></body></html>"""
open("compare_psa_viewer.html", "w").write(HTML.replace("__DATA__", DATA))
print(f"wrote compare_psa_viewer.html ({len(events)} events)")
