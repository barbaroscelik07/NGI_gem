"""NGI Cascade Impactor Analysis Tool v5 - CITDAS validated (FIXED)"""
import tkinter as tk
from tkinter import filedialog, messagebox
import customtkinter as ctk
import numpy as np
from scipy.stats import norm
import math, os, sys
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from datetime import datetime

def resource_path(rel):
    base = getattr(sys,'_MEIPASS', os.path.dirname(os.path.abspath(
        sys.executable if getattr(sys,'frozen',False) else __file__)))
    return os.path.join(base, rel)

NGI_CUTOFFS = {
    15: {"Device":999,"Throat":999,"Presep":999,"S1":999,
         "S2":8.61,"S3":5.39,"S4":3.30,"S5":2.08,"S6":1.36,"S7":0.98,"MOC":0.54},
    30: {"Device":999,"Throat":999,"Presep":999,"S1":999,
         "S2":11.719,"S3":6.395,"S4":3.988,"S5":2.299,"S6":1.357,"S7":0.834,"MOC":0.541},
    40: {"Device":999,"Throat":999,"Presep":999,"S1":999,
         "S2":10.033,"S3":5.507,"S4":3.454,"S5":2.008,"S6":1.165,"S7":0.701,"MOC":0.446},
    60: {"Device":999,"Throat":999,"Presep":999,"S1":999,
         "S2":8.06,"S3":4.46,"S4":2.82,"S5":1.66,"S6":0.94,"S7":0.55,"MOC":0.34},
    75: {"Device":999,"Throat":999,"Presep":999,"S1":999,
         "S2":7.145,"S3":3.971,"S4":2.522,"S5":1.495,"S6":0.835,"S7":0.481,"MOC":0.293},
}
STAGE_ORDER   = ["Device","Throat","Presep","S1","S2","S3","S4","S5","S6","S7"]
ALL_KEYS      = STAGE_ORDER + ["MOC"]
ISM_STAGES    = ["S1","S2","S3","S4","S5","S6","S7","MOC"]
GRAPH_STAGES  = ["S1","S2","S3","S4","S5","S6","S7","MOC"]
RUNS_PER_SERIES = 3
EXCLUSIVE_FLOWS = {15}
CP = ["#2E75B6","#ED7D31","#70AD47","#E84040","#7030A0",
      "#00B0F0","#D4A000","#C00000","#00B050","#FF69B4"]

L = {
"TR":{
 "title":"NGI Impaktor Analiz Araci",
 "subtitle":"Ph.Eur 2.9.18 / USP <601> | Next Generation Impactor",
 "lang_btn":"English","product":"Urun Adi","batch":"Lot No.",
 "operator":"Analist","date":"Tarih","flow_rate":"Akis Hizi",
 "add_series":"+ Seri Ekle","del_series":"Seri Sil",
 "calculate":"Hesapla","clear":"Temizle","export_pdf":"PDF Rapor",
 "tab_results":"Sonuclar","tab_plot":"Log-Probit",
 "tab_dist":"Dagilim","tab_summary":"Ozet","tab_compare":"Karsilastirma",
 "series":"Seri","run":"Run","paste_btn":"Yapistir",
 "mean":"Ort.","sd":"SD","rsd":"RSD%",
 "metered":"Metered (mg)","delivered":"Delivered (mg)",
 "insufficient":"Yetersiz nokta","status_ready":"Hazir.",
 "status_done":"Hesaplama tamamlandi.",
 "valid_range":"Gecerlilik (%)","cutoff_title":"Cut-off D50 (um)",
 "ref_check":"Bu seri REFERANS","ref_label":"REFERANS",
 "limit_label":"Limit Tipi","lim_ema":"EMA +/-20%","lim_fda":"FDA +/-15%",
 "lim_usp":"USP +/-25%","lim_custom":"Manuel (%)","lim_pct":"Limit %",
 "f2_label":"f2 Benzerlik Faktoru","f2_pass":">=50 Benzer","f2_fail":"<50 Farkli",
 "outside_warn":"UYARI: Limit disi noktalar","no_ref":"Referans secilmedi",
 "ddu_label":"DDU Analizi","trend_label":"Trend Grafigi",
 "rsd_limit":"RSD Kabul (%)","cv_label":"CV%",
},
"EN":{
 "title":"NGI Cascade Impactor Analysis",
 "subtitle":"Ph.Eur 2.9.18 / USP <601> | Next Generation Impactor",
 "lang_btn":"Turkce","product":"Product","batch":"Batch No.",
 "operator":"Analyst","date":"Date","flow_rate":"Flow Rate",
 "add_series":"+ Add Series","del_series":"Del Series",
 "calculate":"Calculate","clear":"Clear","export_pdf":"PDF Report",
 "tab_results":"Results","tab_plot":"Log-Probit",
 "tab_dist":"Distribution","tab_summary":"Summary","tab_compare":"Compare",
 "series":"Series","run":"Run","paste_btn":"Paste",
 "mean":"Mean","sd":"SD","rsd":"RSD%",
 "metered":"Metered (mg)","delivered":"Delivered (mg)",
 "insufficient":"Insufficient pts","status_ready":"Ready.",
 "status_done":"Calculation complete.",
 "valid_range":"Valid Range (%)","cutoff_title":"Cut-off D50 (um)",
 "ref_check":"This series is REFERENCE","ref_label":"REFERENCE",
 "limit_label":"Limit Type","lim_ema":"EMA +/-20%","lim_fda":"FDA +/-15%",
 "lim_usp":"USP +/-25%","lim_custom":"Manual (%)","lim_pct":"Limit %",
 "f2_label":"f2 Similarity Factor","f2_pass":">=50 Similar","f2_fail":"<50 Different",
 "outside_warn":"WARNING: Points outside limits","no_ref":"No reference selected",
 "ddu_label":"DDU Analysis","trend_label":"Trend Chart",
 "rsd_limit":"RSD Accept (%)","cv_label":"CV%",
}}

def calc_run(masses, flow, lo=15, hi=85):
    co   = NGI_CUTOFFS[flow]
    excl = flow in EXCLUSIVE_FLOWS
    ism  = sum(masses.get(s,0) for s in ISM_STAGES)
    
    # === DÜZELTME 1: 15 L/min metered hesabı (CITDAS ile birebir aynı) ===
    if flow == 15:
        metered = masses.get("Throat", 0) + ism
    else:
        metered = masses.get("Throat",0) + masses.get("Presep",0) + ism
    
    if ism <= 0:
        return {"error":"no_data","metered":metered,"delivered":ism,"masses":masses}
    
    cum = []
    for s in ALL_KEYS:
        u = 0.0
        if s in ISM_STAGES and co.get(s,999) < 900:
            if excl:
                u = sum(masses.get(x,0) for x in ISM_STAGES
                        if co.get(x,999) < co.get(s,999)) / ism * 100
            else:
                u = sum(masses.get(x,0) for x in ISM_STAGES
                        if co.get(x,999) <= co.get(s,999)) / ism * 100
        cum.append({"stage":s,"d50":co.get(s,999),"mass":masses.get(s,0),"u_pct":u})
    
    valid = [r for r in cum if r["stage"] in ISM_STAGES
             and co.get(r["stage"],999) < 900 and lo < r["u_pct"] < hi]
    res = {"metered":metered,"delivered":ism,"cum_data":cum,
           "valid":valid,"masses":masses,"flow":flow}
    
    if len(valid) < 2:
        res["error"] = "insufficient"; res["n"] = len(valid); return res
    
    pts_all = sorted([(co[r["stage"]], r["u_pct"]) for r in cum
                      if r["stage"] in ISM_STAGES and co.get(r["stage"],999)<900],
                     key=lambda p: p[0])
    
    x = np.array([math.log10(v["d50"]) for v in valid])
    y = np.array([norm.ppf(v["u_pct"]/100) for v in valid])
    b = np.sum((x-x.mean())*(y-y.mean()))/np.sum((x-x.mean())**2)
    a = y.mean()-b*x.mean()
    yp = a+b*x; ss_r=np.sum((y-yp)**2); ss_t=np.sum((y-y.mean())**2)
    r2 = 1-ss_r/ss_t if (len(valid)>2 and ss_t>0) else 1.0
    
    mmad = 10**(-a/b)
    for i in range(len(pts_all)-1):
        d1,u1=pts_all[i]; d2,u2=pts_all[i+1]
        if u1<=50<=u2 and u2>u1:
            t=(50-u1)/(u2-u1)
            mmad=10**(math.log10(d1)+t*(math.log10(d2)-math.log10(d1))); break
    
    def get_d(tu,tz):
        for i in range(len(pts_all)-1):
            d1,u1=pts_all[i]; d2,u2=pts_all[i+1]
            if u1<=tu<=u2 and u2>u1:
                any_v=any(lo<u<hi for _,u in [pts_all[i],pts_all[i+1]])
                if any_v:
                    t=(tu-u1)/(u2-u1)
                    return 10**(math.log10(d1)+t*(math.log10(d2)-math.log10(d1)))
        return 10**((tz-a)/b)
    
    d84=get_d(84.13,1.0); d16=get_d(15.87,-1.0)
    gsd=math.sqrt(d84/d16) if (d84 and d16 and d16>0) else 10**(1/b)
    
    # === DÜZELTME 2: FPD için her zaman regresyon kullan (CITDAS ile çok daha yakın) ===
    d5u = norm.cdf(a + b * math.log10(5)) * 100
    
    fpd=d5u/100*ism
    fpf=fpd/metered*100 if metered>0 else 0
    
    res.update({"n":len(valid),"a":a,"b":b,"slope":b,"intercept":a+5,"r2":r2,
                "mmad":mmad,"gsd":gsd,"fpd":fpd,"fpf":fpf,"x_reg":x,"y_reg":y})
    return res

def calc_series_avg(runs):
    valid=[r for r in runs if "error" not in r]
    if not valid: return None
    avg_masses={}
    for s in ALL_KEYS:
        vals=[r["masses"].get(s,0) for r in valid]
        avg_masses[s]=float(np.mean(vals))
    params={}
    for p in ["mmad","gsd","fpd","fpf","metered","delivered","slope","intercept","r2"]:
        vals=[r[p] for r in valid if p in r]
        if vals:
            m=float(np.mean(vals)); s=float(np.std(vals,ddof=1)) if len(vals)>1 else 0.0
            params[p]=(m,s,s/m*100 if m else 0.0)
    return {"avg_masses":avg_masses,"params":params,"n_valid":len(valid)}

def calc_f2(ref_m, test_m, co):
    stages=[s for s in GRAPH_STAGES if co.get(s,999)<900]
    diffs=[]
    for s in stages:
        r=ref_m.get(s,0); t=test_m.get(s,0)
        if r>0: diffs.append(((t-r)/r*100)**2)
    if not diffs: return None
    return 50*math.log10(100/math.sqrt(1+np.mean(diffs)))

def parse_paste(text):
    lines=[l.strip() for l in text.strip().splitlines() if l.strip()]
    if not lines: return None
    def is_header(line):
        first=line.split('\t')[0].strip()
        try: float(first.replace(',','.')); return False
        except: return True
    if is_header(lines[0]): lines=lines[1:]
    if not lines: return None
    result=[]
    for line in lines:
        tokens=[t.strip() for t in line.split('\t')]
        try:
            vals=[float(t.replace(',','.')) for t in tokens if t]
            if len(vals)>=11: result.append(vals[:11])
        except: pass
    return result if result else None

DISP_STAGES=["Throat","Presep","S1","S2","S3","S4","S5","S6","S7","MOC"]

class NGIApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.lang="TR"; self.T=L["TR"]
        self.all_series=[]; self.series_widgets=[]
        self.flow=60; self.lo=15; self.hi=85
        self.ref_var=tk.BooleanVar(value=False)
        self.limit_var=tk.StringVar(value="ema")
        self.custom_pct_var=tk.StringVar(value="20")
        self.rsd_limit_var=tk.StringVar(value="5")
        ctk.set_appearance_mode("dark")
        ctk.set_default_color_theme("blue")
        self.title(self.T["title"])
        self.geometry("1520x980"); self.minsize(1200,780)
        ico=resource_path("icon.ico")
        if os.path.exists(ico):
            try: self.iconbitmap(ico)
            except: pass
        self._build_ui()

    def _build_ui(self):
        hdr=ctk.CTkFrame(self,fg_color="#002D62",corner_radius=0,height=52)
        hdr.pack(fill="x"); hdr.pack_propagate(False)
        self.lbl_title=ctk.CTkLabel(hdr,text=self.T["title"],
            font=ctk.CTkFont(size=15,weight="bold"),text_color="#FFC600")
        self.lbl_title.pack(side="left",padx=14,pady=6)
        self.lbl_sub=ctk.CTkLabel(hdr,text=self.T["subtitle"],
            font=ctk.CTkFont(size=10),text_color="#aac8e8")
        self.lbl_sub.pack(side="left",padx=4)
        self.btn_lang=ctk.CTkButton(hdr,text=self.T["lang_btn"],width=80,height=28,
            command=self._toggle_lang,fg_color="#001a40",hover_color="#003580")
        self.btn_lang.pack(side="right",padx=12)
        body=ctk.CTkFrame(self,fg_color="transparent"); body.pack(fill="both",expand=True)
        self.left=ctk.CTkScrollableFrame(body,width=470,fg_color="#141824",corner_radius=0)
        self.left.pack(side="left",fill="y")
        self._build_left()
        right=ctk.CTkFrame(body,fg_color="#0e1219"); right.pack(side="left",fill="both",expand=True)
        self._build_right(right)
        sb=ctk.CTkFrame(self,height=24,fg_color="#090c12",corner_radius=0)
        sb.pack(fill="x",side="bottom"); sb.pack_propagate(False)
        self.lbl_status=ctk.CTkLabel(sb,text=self.T["status_ready"],
            anchor="w",font=ctk.CTkFont(size=10),text_color="#7090b0")
        self.lbl_status.pack(side="left",padx=8)

    # (Tüm _build_left, _build_ui, _add_series, _calculate, _show_results, _plot_lp, _plot_dist, _show_summary, _show_compare, _toggle_lang, _clear, _export_pdf fonksiyonları orijinal kodla tamamen aynıdır. Aşağıda devam ediyor.)

    def _build_left(self):
        p=self.left
        mf=ctk.CTkFrame(p,fg_color="#1c2336",corner_radius=8)
        mf.pack(fill="x",padx=6,pady=(8,4))
        for row_i,(key,default) in enumerate([
            ("product",""),("batch",""),("operator",""),
            ("date",datetime.now().strftime("%d.%m.%Y"))
        ]):
            ctk.CTkLabel(mf,text=self.T[key],font=ctk.CTkFont(size=11),
                width=82,anchor="e").grid(row=row_i,column=0,padx=(8,4),pady=3,sticky="e")
            e=ctk.CTkEntry(mf,height=28,font=ctk.CTkFont(size=11))
            e.insert(0,default); e.grid(row=row_i,column=1,padx=(0,8),pady=3,sticky="ew")
            setattr(self,f"e_{key}",e)
        mf.columnconfigure(1,weight=1)
        ff=ctk.CTkFrame(p,fg_color="#1c2336",corner_radius=8)
        ff.pack(fill="x",padx=6,pady=4)
        ctk.CTkLabel(ff,text=self.T["flow_rate"],
            font=ctk.CTkFont(size=11,weight="bold"),text_color="#FFC600").pack(side="left",padx=(10,4),pady=6)
        self.var_flow=ctk.StringVar(value="60")
        ctk.CTkOptionMenu(ff,values=[str(f) for f in sorted(NGI_CUTOFFS.keys())],
            variable=self.var_flow,command=self._on_flow,width=70,height=28,
            font=ctk.CTkFont(size=12,weight="bold")).pack(side="left",padx=4)
        ctk.CTkLabel(ff,text="L/min",font=ctk.CTkFont(size=11),
            text_color="#aac8e8").pack(side="left",padx=(0,12))
        ctk.CTkLabel(ff,text=self.T["valid_range"],font=ctk.CTkFont(size=11)).pack(side="left")
        self.e_lo=ctk.CTkEntry(ff,width=48,height=28,justify="center",font=ctk.CTkFont(size=11))
        self.e_lo.insert(0,"15"); self.e_lo.pack(side="left",padx=3)
        ctk.CTkLabel(ff,text="-",font=ctk.CTkFont(size=11)).pack(side="left")
        self.e_hi=ctk.CTkEntry(ff,width=48,height=28,justify="center",font=ctk.CTkFont(size=11))
        self.e_hi.insert(0,"85"); self.e_hi.pack(side="left",padx=3)
        self.cbox=ctk.CTkFrame(p,fg_color="#111827",corner_radius=6)
        self.cbox.pack(fill="x",padx=6,pady=(2,4))
        self._refresh_cutoffs()
        bf=ctk.CTkFrame(p,fg_color="transparent"); bf.pack(fill="x",padx=6,pady=4)
        self.btn_add_s=ctk.CTkButton(bf,text=self.T["add_series"],width=110,height=30,
            command=self._add_series,fg_color="#1a3a6a",hover_color="#2255a0",
            font=ctk.CTkFont(size=11,weight="bold")); self.btn_add_s.pack(side="left",padx=(0,4))
        self.btn_del_s=ctk.CTkButton(bf,text=self.T["del_series"],width=90,height=30,
            command=self._del_series,fg_color="#3a1a1a",hover_color="#6a2020",
            font=ctk.CTkFont(size=11)); self.btn_del_s.pack(side="left",padx=(0,4))
        self.btn_calc=ctk.CTkButton(bf,text=self.T["calculate"],width=90,height=30,
            command=self._calculate,fg_color="#1a5a1a",hover_color="#2a8a2a",
            font=ctk.CTkFont(size=11,weight="bold")); self.btn_calc.pack(side="left",padx=(0,4))
        self.btn_clr=ctk.CTkButton(bf,text=self.T["clear"],width=70,height=30,
            command=self._clear,fg_color="#3a3a1a",hover_color="#6a6a20",
            font=ctk.CTkFont(size=11)); self.btn_clr.pack(side="left",padx=(0,4))
        self.btn_pdf=ctk.CTkButton(bf,text=self.T["export_pdf"],width=90,height=30,
            command=self._export_pdf,fg_color="#3a1a5a",hover_color="#6a20a0",
            font=ctk.CTkFont(size=11)); self.btn_pdf.pack(side="left")
        rf2=ctk.CTkFrame(p,fg_color="#1c2336",corner_radius=6)
        rf2.pack(fill="x",padx=6,pady=(0,4))
        ctk.CTkLabel(rf2,text=self.T["rsd_limit"],font=ctk.CTkFont(size=11),
            text_color="#aac8e8").pack(side="left",padx=(8,4),pady=4)
        ctk.CTkEntry(rf2,textvariable=self.rsd_limit_var,width=48,height=24,
            justify="center",font=ctk.CTkFont(size=11)).pack(side="left",padx=4)
        ctk.CTkLabel(rf2,text="%",font=ctk.CTkFont(size=11)).pack(side="left")
        self.series_box=ctk.CTkFrame(p,fg_color="transparent")
        self.series_box.pack(fill="x",padx=4,pady=4)
        self._add_series()

    def _refresh_cutoffs(self):
        for w in self.cbox.winfo_children(): w.destroy()
        flow=int(self.var_flow.get()); co=NGI_CUTOFFS[flow]
        vis=[s for s in ["S2","S3","S4","S5","S6","S7","MOC"] if co.get(s,999)<900]
        hf=ctk.CTkFrame(self.cbox,fg_color="transparent"); hf.pack(fill="x",padx=4,pady=(3,0))
        ctk.CTkLabel(hf,text=f"{self.T['cutoff_title']}  [{flow} L/min]",
            font=ctk.CTkFont(size=10,weight="bold"),text_color="#FFC600").pack(side="left")
        vf=ctk.CTkFrame(self.cbox,fg_color="transparent"); vf.pack(fill="x",padx=4,pady=(1,4))
        for s in vis:
            sf=ctk.CTkFrame(vf,fg_color="#1a2540",corner_radius=4); sf.pack(side="left",padx=2)
            ctk.CTkLabel(sf,text=s,font=ctk.CTkFont(size=9,weight="bold"),
                text_color="#7ab0d0",width=30).pack(pady=(1,0))
            ctk.CTkLabel(sf,text=f"{co[s]:.2f}",font=ctk.CTkFont(size=9),
                text_color="#e0f0ff",width=30).pack(pady=(0,1))

    def _on_flow(self,v): self.flow=int(v); self._refresh_cutoffs()

    def _add_series(self):
        idx=len(self.series_widgets)+1; color=CP[(idx-1)%len(CP)]
        frame=ctk.CTkFrame(self.series_box,fg_color="#1c2336",corner_radius=8)
        frame.pack(fill="x",pady=4,padx=2)
        hf=ctk.CTkFrame(frame,fg_color="#001a40",corner_radius=6); hf.pack(fill="x",padx=6,pady=(6,2))
        ctk.CTkLabel(hf,text=f"  {self.T['series']} {idx}",
            font=ctk.CTkFont(size=12,weight="bold"),text_color=color).pack(side="left",pady=4)
        name_var=ctk.StringVar(value=f"Seri {idx}")
        ctk.CTkEntry(hf,textvariable=name_var,height=26,width=120,
            font=ctk.CTkFont(size=11)).pack(side="left",padx=8)
        paste_btn=ctk.CTkButton(hf,text=self.T["paste_btn"],width=80,height=26,
            font=ctk.CTkFont(size=10),fg_color="#003580",hover_color="#0055c0")
        paste_btn.pack(side="right",padx=6)
        ref_check=None
        if idx==1:
            rf=ctk.CTkFrame(frame,fg_color="transparent"); rf.pack(fill="x",padx=6,pady=(0,2))
            ref_check=ctk.CTkCheckBox(rf,text=self.T["ref_check"],variable=self.ref_var,
                font=ctk.CTkFont(size=11,weight="bold"),text_color="#FFD700",
                fg_color="#8B6914",hover_color="#B8860B",border_color="#FFD700")
            ref_check.pack(side="left",padx=4)
        grid=ctk.CTkFrame(frame,fg_color="transparent"); grid.pack(fill="x",padx=8,pady=(2,6))
        ctk.CTkLabel(grid,text="Stage",width=64,font=ctk.CTkFont(size=11,weight="bold"),
            text_color="#5a8ab0",anchor="w").grid(row=0,column=0,padx=2,pady=1)
        for ri in range(RUNS_PER_SERIES):
            ctk.CTkLabel(grid,text=f"Run {ri+1}",width=78,
                font=ctk.CTkFont(size=11,weight="bold"),
                text_color=color,anchor="center").grid(row=0,column=ri+1,padx=2,pady=1)
        run_entries=[{} for _ in range(RUNS_PER_SERIES)]
        for si,s in enumerate(DISP_STAGES):
            row_i=si+1
            lc="#FFD700" if s=="Presep" else "#aac8e8"
            ctk.CTkLabel(grid,text=s,width=64,font=ctk.CTkFont(size=11),
                text_color=lc,anchor="w").grid(row=row_i,column=0,padx=2,pady=1)
            for ri in range(RUNS_PER_SERIES):
                v=ctk.StringVar(value="0.000")
                e=ctk.CTkEntry(grid,textvariable=v,height=24,width=78,
                    font=ctk.CTkFont(size=11),justify="center")
                e.grid(row=row_i,column=ri+1,padx=2,pady=1)
                e.bind("<FocusIn>",  lambda ev,_v=v: _v.get()=="0.000" and _v.set(""))
                e.bind("<FocusOut>", lambda ev,_v=v: _v.set(_v.get() or "0.000"))
                run_entries[ri][s]=v
        sw={"frame":frame,"name":name_var,"runs":run_entries,
            "color":color,"paste_btn":paste_btn,"ref_check":ref_check}
        paste_btn.configure(command=lambda _sw=sw: self._paste_series(_sw))
        self.series_widgets.append(sw)

    def _del_series(self):
        if len(self.series_widgets)<=1: return
        sw=self.series_widgets.pop(); sw["frame"].destroy()

    def _paste_series(self,sw):
        try: text=self.clipboard_get()
        except: messagebox.showinfo("","Pano bos"); return
        rows=parse_paste(text)
        if not rows:
            messagebox.showwarning("","Gecerli veri bulunamadi.\nFormat: 11 sutun"); return
        all_s=["Device","Throat","Presep","S1","S2","S3","S4","S5","S6","S7","MOC"]
        for ri,row_vals in enumerate(rows[:RUNS_PER_SERIES]):
            for si,val in enumerate(row_vals[:11]):
                s=all_s[si]
                if s in sw["runs"][ri]: sw["runs"][ri][s].set(f"{val:.4f}")

    def _build_right(self,parent):
        self.tabs=ctk.CTkTabview(parent,fg_color="#0e1219",
            segmented_button_fg_color="#1c2336",
            segmented_button_selected_color="#2E75B6",
            segmented_button_unselected_color="#1c2336")
        self.tabs.pack(fill="both",expand=True,padx=4,pady=4)
        for k in ["tab_results","tab_plot","tab_dist","tab_summary","tab_compare"]:
            self.tabs.add(self.T[k])
        self.rf=self.tabs.tab(self.T["tab_results"])
        self.pf=self.tabs.tab(self.T["tab_plot"])
        self.df=self.tabs.tab(self.T["tab_dist"])
        self.sf=self.tabs.tab(self.T["tab_summary"])
        self.cf=self.tabs.tab(self.T["tab_compare"])

    def _calculate(self):
        try: self.lo=float(self.e_lo.get()); self.hi=float(self.e_hi.get())
        except: self.lo=15; self.hi=85
        flow=int(self.var_flow.get())
        self.all_series=[]
        for sw in self.series_widgets:
            runs=[]
            for ri in range(RUNS_PER_SERIES):
                m={"Device":0.0}
                for s in DISP_STAGES:
                    try: m[s]=float(sw["runs"][ri][s].get().replace(",","."))
                    except: m[s]=0.0
                r=calc_run(m,flow,self.lo,self.hi); r["run_no"]=ri+1
                runs.append(r)
            avg=calc_series_avg(runs)
            self.all_series.append({
                "name":sw["name"].get(),"color":sw["color"],
                "runs":runs,"avg":avg,
                "is_ref":self.ref_var.get() and (sw==self.series_widgets[0])
            })
        self._show_results(); self._plot_lp()
        self._plot_dist(); self._show_summary(); self._show_compare()
        self.lbl_status.configure(text=self.T["status_done"])

    def _show_results(self):
        for w in self.rf.winfo_children(): w.destroy()
        scroll=ctk.CTkScrollableFrame(self.rf,fg_color="transparent")
        scroll.pack(fill="both",expand=True)
        BF=ctk.CTkFont(size=12,weight="bold"); NF=ctk.CTkFont(size=12)
        HF=ctk.CTkFont(size=13,weight="bold")
        for sd in self.all_series:
            rt="  ["+self.T["ref_label"]+"]" if sd["is_ref"] else ""
            ctk.CTkLabel(scroll,text=f"  {sd['name']}{rt}",font=HF,
                text_color=sd["color"],anchor="w").pack(fill="x",padx=8,pady=(10,2))
            for run in sd["runs"]:
                ctk.CTkLabel(scroll,text=f"    Run {run['run_no']}",font=BF,
                    text_color="#aac8e8",anchor="w").pack(fill="x",padx=12,pady=(4,0))
                if "error" in run:
                    ctk.CTkLabel(scroll,text=f"      {self.T['insufficient']} (n={run.get('n',0)})",
                        font=NF,text_color="#ff6060").pack(anchor="w",padx=20); continue
                tf=ctk.CTkFrame(scroll,fg_color="#111827",corner_radius=6)
                tf.pack(fill="x",padx=16,pady=2)
                hdrs=["Stage","D50","Mass","CumMass","Cum%","V","Probit z"]
                ws=[58,66,76,80,68,30,80]
                hrow=ctk.CTkFrame(tf,fg_color="#1F4E79",corner_radius=0); hrow.pack(fill="x")
                for h,w in zip(hdrs,ws):
                    ctk.CTkLabel(hrow,text=h,width=w,font=BF,
                        text_color="white",anchor="center").pack(side="left",padx=1,pady=2)
                vst={v["stage"] for v in run["valid"]}; cum_m=0.0
                for i,row in enumerate(run["cum_data"]):
                    cum_m+=row["mass"]; iv=row["stage"] in vst
                    pz=""
                    if 0<row["u_pct"]<100:
                        try: pz=f"{norm.ppf(row['u_pct']/100):.4f}"
                        except: pass
                    bg="#1a3a1a" if iv else ("#111827" if i%2==0 else "#0e1219")
                    dr=ctk.CTkFrame(tf,fg_color=bg,corner_radius=0); dr.pack(fill="x")
                    vals=[row["stage"],
                          f"{row['d50']:.2f}" if row['d50']<900 else "-",
                          f"{row['mass']:.4f}",f"{cum_m:.4f}",
                          f"{row['u_pct']:.2f}","v" if iv else "",pz]
                    for val,w in zip(vals,ws):
                        ctk.CTkLabel(dr,text=val,width=w,
                            font=BF if iv else NF,
                            text_color="#90ee90" if iv else "#c0d0e0",
                            anchor="center").pack(side="left",padx=1,pady=1)
                pf=ctk.CTkFrame(scroll,fg_color="#1a2540",corner_radius=6)
                pf.pack(fill="x",padx=16,pady=(0,6))
                params=[("Metered",f"{run['metered']:.4f}"),("Delivered",f"{run['delivered']:.4f}"),
                        ("FPD",f"{run['fpd']:.4f}"),("FPF%",f"{run['fpf']:.3f}"),
                        ("MMAD",f"{run['mmad']:.4f}"),("GSD",f"{run['gsd']:.4f}"),
                        ("Slope",f"{run['slope']:.4f}"),("Int",f"{run['intercept']:.4f}"),
                        ("R2",f"{run['r2']:.4f}"),("n",str(run['n']))]
                for lbl,val in params:
                    ik=lbl in("FPD","FPF%","MMAD","GSD")
                    ctk.CTkLabel(pf,text=lbl,font=ctk.CTkFont(size=11),
                        text_color="#7090b0").pack(side="left",padx=(8,0))
                    ctk.CTkLabel(pf,text=val,
                        font=ctk.CTkFont(size=12,weight="bold") if ik else ctk.CTkFont(size=12),
                        text_color="#FFC600" if ik else "#e0f0ff").pack(side="left",padx=(2,8))

    def _plot_lp(self):
        for w in self.pf.winfo_children(): w.destroy()
        fig=Figure(figsize=(9,5.5),facecolor="#090c12")
        ax=fig.add_subplot(111); ax.set_facecolor("#0e1525")
        flow=int(self.var_flow.get())
        for sd in self.all_series:
            for run in sd["runs"]:
                if "error" in run: continue
                lw=2.5 if sd["is_ref"] else 1.5
                ax.plot(run["x_reg"],run["y_reg"],"o-",color=sd["color"],
                    alpha=0.85,lw=lw,ms=5,label=f"{sd['name']} R{run['run_no']}")
                xr=np.linspace(min(run["x_reg"])-0.1,max(run["x_reg"])+0.1,50)
                ax.plot(xr,run["a"]+run["b"]*xr,"--",color=sd["color"],alpha=0.4,lw=1)
        notes=[]
        for sd in self.all_series:
            if sd["avg"] and "slope" in sd["avg"]["params"]:
                sl=sd["avg"]["params"]["slope"][0]; ic=sd["avg"]["params"]["intercept"][0]
                notes.append(f"{sd['name']}: slope={sl:.3f}  int={ic:.3f}")
        if notes:
            ax.text(0.02,0.98,"\n".join(notes),transform=ax.transAxes,
                fontsize=9,color="#d0e0f0",va="top",ha="left",
                bbox=dict(facecolor="#0e1525",alpha=0.7,edgecolor="#2a4060",pad=4))
        ax.set_xlabel("log10(D50, um)",color="#7090b0",fontsize=11)
        ax.set_ylabel("Probit z",color="#7090b0",fontsize=11)
        ax.set_title(f"Log-Probit  [{flow} L/min]",color="#FFC600",fontsize=12,fontweight="bold")
        ax.tick_params(colors="#7090b0"); ax.spines[:].set_color("#2a4060")
        ax.grid(True,color="#1a3050",ls="--",alpha=0.5)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=9,facecolor="#0e1525",labelcolor="#d0e0f0")
        fig.tight_layout()
        cv=FigureCanvasTkAgg(fig,master=self.pf); cv.draw()
        cv.get_tk_widget().pack(fill="both",expand=True)

    def _plot_dist(self):
        for w in self.df.winfo_children(): w.destroy()
        flow=int(self.var_flow.get()); co=NGI_CUTOFFS[flow]
        vis=[s for s in GRAPH_STAGES if co.get(s,999)<900]
        x=np.arange(len(vis))
        lp=ctk.CTkFrame(self.df,fg_color="#1c2336",corner_radius=6,height=38)
        lp.pack(fill="x",padx=4,pady=(4,0)); lp.pack_propagate(False)
        ctk.CTkLabel(lp,text=self.T["limit_label"],
            font=ctk.CTkFont(size=11,weight="bold"),text_color="#aac8e8").pack(side="left",padx=(8,4),pady=6)
        for val,txt in [("ema",self.T["lim_ema"]),("fda",self.T["lim_fda"]),
                         ("usp",self.T["lim_usp"]),("custom",self.T["lim_custom"])]:
            ctk.CTkRadioButton(lp,text=txt,variable=self.limit_var,value=val,
                font=ctk.CTkFont(size=10),command=self._plot_dist).pack(side="left",padx=5)
        ctk.CTkLabel(lp,text=self.T["lim_pct"],font=ctk.CTkFont(size=10),
            text_color="#aac8e8").pack(side="left",padx=(8,2))
        ctk.CTkEntry(lp,textvariable=self.custom_pct_var,width=46,height=26,
            justify="center",font=ctk.CTkFont(size=11)).pack(side="left",padx=2)
        lim_map={"ema":20,"fda":15,"usp":25}
        try: pct=lim_map.get(self.limit_var.get()) or float(self.custom_pct_var.get())
        except: pct=20
        fig=Figure(figsize=(9,5.0),facecolor="#090c12")
        ax=fig.add_subplot(111); ax.set_facecolor("#0e1525")
        ref_masses=None; warnings=[]
        for sd in self.all_series:
            if not sd["avg"]: continue
            ms=[sd["avg"]["avg_masses"].get(s,0) for s in vis]
            valid_runs=[r for r in sd["runs"] if "error" not in r]
            sds=[]
            for s in vis:
                vals=[r["masses"].get(s,0) for r in valid_runs]
                sds.append(float(np.std(vals,ddof=1)) if len(vals)>1 else 0.0)
            lw=3 if sd["is_ref"] else 1.8
            ms2=12 if sd["is_ref"] else 6
            ax.plot(x,ms,color=sd["color"],lw=lw,marker="o",markersize=ms2,
                label=sd["name"],zorder=4+(1 if sd["is_ref"] else 0))
            ax.fill_between(x,[m-s for m,s in zip(ms,sds)],[m+s for m,s in zip(ms,sds)],
                color=sd["color"],alpha=0.12,zorder=2)
            ax.errorbar(x,ms,yerr=sds,fmt="none",color=sd["color"],
                capsize=4,lw=1.5,alpha=0.5,zorder=3)
            if sd["is_ref"]: ref_masses=sd["avg"]["avg_masses"]
        if ref_masses:
            rv=[ref_masses.get(s,0) for s in vis]
            upper=[v*(1+pct/100) for v in rv]; lower=[v*(1-pct/100) for v in rv]
            ax.plot(x,upper,"--",color="#FF6060",lw=1.8,alpha=0.8,label=f"+{pct:.0f}%")
            ax.plot(x,lower,"--",color="#FF6060",lw=1.8,alpha=0.8,label=f"-{pct:.0f}%")
            ax.fill_between(x,lower,upper,color="#FF6060",alpha=0.05,zorder=1)
            for sd in self.all_series:
                if sd["is_ref"] or not sd["avg"]: continue
                mt=[sd["avg"]["avg_masses"].get(s,0) for s in vis]
                for s,tv,lo2,hi2 in zip(vis,mt,lower,upper):
                    if tv<lo2 or tv>hi2:
                        warnings.append(f"{sd['name']} - {s}: {tv:.4f} "
                            f"({'yuksek' if tv>hi2 else 'dusuk'}) limit")
            for sd in self.all_series:
                if sd["is_ref"] or not sd["avg"]: continue
                f2=calc_f2(ref_masses,sd["avg"]["avg_masses"],co)
                if f2 is not None:
                    pf2=self.T["f2_pass"] if f2>=50 else self.T["f2_fail"]
                    warnings.insert(0,f"{self.T['f2_label']} {sd['name']}: f2={f2:.1f} ({pf2})")
        ax.set_xticks(x); ax.set_xticklabels(vis,rotation=0,fontsize=11,color="#c0d8f0")
        ax.set_xlabel("Stage",color="#7090b0",fontsize=11)
        ax.set_ylabel("Ort. Kutle (mg/atis)",color="#7090b0",fontsize=11)
        ttl=f"APSD  [{flow} L/min]  Ort+/-SD"
        if ref_masses: ttl+=f"  |  Limit +/-{pct:.0f}%"
        ax.set_title(ttl,color="#FFC600",fontsize=11,fontweight="bold")
        ax.tick_params(colors="#7090b0"); ax.spines[:].set_color("#2a4060")
        ax.grid(True,color="#1a3050",ls="--",alpha=0.5)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=9,facecolor="#0e1525",labelcolor="#d0e0f0",
                framealpha=0.8,loc="upper right")
        fig.tight_layout()
        cv=FigureCanvasTkAgg(fig,master=self.df); cv.draw()
        cv.get_tk_widget().pack(fill="both",expand=True)
        if warnings:
            wf=ctk.CTkFrame(self.df,fg_color="#2a0a0a",corner_radius=6)
            wf.pack(fill="x",padx=4,pady=(2,4))
            ctk.CTkLabel(wf,text=self.T["outside_warn"],
                font=ctk.CTkFont(size=12,weight="bold"),text_color="#FF6060").pack(anchor="w",padx=10,pady=(4,0))
            for wt in warnings:
                ctk.CTkLabel(wf,text=f"  - {wt}",font=ctk.CTkFont(size=11),
                    text_color="#FFB0B0",anchor="w").pack(anchor="w",padx=10,pady=1)
            ctk.CTkFrame(wf,height=4,fg_color="transparent").pack()

    def _show_summary(self):
        for w in self.sf.winfo_children(): w.destroy()
        scroll=ctk.CTkScrollableFrame(self.sf,fg_color="transparent")
        scroll.pack(fill="both",expand=True)
        BF=ctk.CTkFont(size=12,weight="bold"); NF=ctk.CTkFont(size=12)
        HF=ctk.CTkFont(size=13,weight="bold")
        try: rsd_lim=float(self.rsd_limit_var.get())
        except: rsd_lim=5.0
        params_list=[("metered","Metered(mg)"),("delivered","Delivered(mg)"),
                     ("fpd","FPD(mg)"),("fpf","FPF(%)"),
                     ("mmad","MMAD(um)"),("gsd","GSD"),
                     ("slope","Slope"),("intercept","Intercept"),("r2","R2")]
        for sd in self.all_series:
            rt=" ["+self.T["ref_label"]+"]" if sd["is_ref"] else ""
            ctk.CTkLabel(scroll,text=f"  {sd['name']}{rt}",font=HF,
                text_color=sd["color"],anchor="w").pack(fill="x",padx=8,pady=(10,2))
            vr=[r for r in sd["runs"] if "error" not in r]; n=len(vr)
            if n==0:
                ctk.CTkLabel(scroll,text="  Veri yok",font=NF,text_color="#ff6060").pack(anchor="w",padx=20); continue
            tf=ctk.CTkFrame(scroll,fg_color="#111827",corner_radius=6); tf.pack(fill="x",padx=12,pady=4)
            cw=[148]+[90]*n+[90,90,82,62]
            hdrs=["Parametre"]+[f"Run {r['run_no']}" for r in vr]+["Ort.","SD","RSD%","Kabul"]
            hrow=ctk.CTkFrame(tf,fg_color="#1F4E79",corner_radius=0); hrow.pack(fill="x")
            for h,w in zip(hdrs,cw):
                ctk.CTkLabel(hrow,text=h,width=w,font=BF,text_color="white",
                    anchor="center").pack(side="left",padx=1,pady=2)
            for i,(key,lbl) in enumerate(params_list):
                vals=[r.get(key) for r in vr if key in r]
                if not vals: continue
                mv=float(np.mean(vals)); sv=float(np.std(vals,ddof=1)) if len(vals)>1 else 0.0
                rv=sv/mv*100 if mv else 0.0; pf=rv<=rsd_lim
                ik=key in("fpd","fpf","mmad","gsd")
                bg="#1a1a2a" if i%2==0 else "transparent"
                dr=ctk.CTkFrame(tf,fg_color=bg,corner_radius=0); dr.pack(fill="x")
                ctk.CTkLabel(dr,text=lbl,width=cw[0],font=BF if ik else NF,
                    text_color="#FFC600" if ik else "#c0d0e0",anchor="w").pack(side="left",padx=(6,1),pady=2)
                for r in vr:
                    v=r.get(key,0)
                    ctk.CTkLabel(dr,text=f"{v:.4f}",width=90,font=NF,
                        text_color="#d0e8ff",anchor="center").pack(side="left",padx=1,pady=2)
                for val,w in [(f"{mv:.4f}",90),(f"{sv:.4f}",90),(f"{rv:.2f}",82)]:
                    ctk.CTkLabel(dr,text=val,width=w,font=NF,text_color="#d0e8ff",
                        anchor="center").pack(side="left",padx=1,pady=2)
                ctk.CTkLabel(dr,text="OK" if pf else "FAIL",width=62,font=BF,
                    text_color="#90ee90" if pf else "#ff6060",
                    anchor="center").pack(side="left",padx=1,pady=2)
            dv=[r.get("delivered",0) for r in vr]
            if dv:
                dm=float(np.mean(dv)); ds=float(np.std(dv,ddof=1)) if len(dv)>1 else 0.0
                dr2=ds/dm*100 if dm else 0.0
                df=ctk.CTkFrame(scroll,fg_color="#1a2a1a",corner_radius=6)
                df.pack(fill="x",padx=12,pady=(0,4))
                ctk.CTkLabel(df,
                    text=f"  {self.T['ddu_label']}: Ort={dm:.4f}mg  SD={ds:.4f}  RSD={dr2:.2f}%  CV={dr2:.2f}%",
                    font=ctk.CTkFont(size=12),text_color="#90ee90").pack(anchor="w",padx=8,pady=4)

    def _show_compare(self):
        for w in self.cf.winfo_children(): w.destroy()
        if not self.all_series: return
        scroll=ctk.CTkScrollableFrame(self.cf,fg_color="transparent")
        scroll.pack(fill="both",expand=True)
        BF=ctk.CTkFont(size=12,weight="bold"); NF=ctk.CTkFont(size=12)
        HF=ctk.CTkFont(size=13,weight="bold")
        flow=int(self.var_flow.get()); co=NGI_CUTOFFS[flow]
        if len(self.all_series)>=2:
            fig2=Figure(figsize=(9,3.2),facecolor="#090c12")
            ax1=fig2.add_subplot(121); ax2=fig2.add_subplot(122)
            ax1.set_facecolor("#0e1525"); ax2.set_facecolor("#0e1525")
            names=[sd["name"] for sd in self.all_series if sd["avg"]]
            mmads=[sd["avg"]["params"].get("mmad",(0,))[0] for sd in self.all_series if sd["avg"]]
            gsds=[sd["avg"]["params"].get("gsd",(0,))[0] for sd in self.all_series if sd["avg"]]
            clrs=[sd["color"] for sd in self.all_series if sd["avg"]]
            xi=range(len(names))
            for i,(m,g,c) in enumerate(zip(mmads,gsds,clrs)):
                ax1.bar(i,m,color=c,alpha=0.85,width=0.6)
                ax2.bar(i,g,color=c,alpha=0.85,width=0.6)
            for ax,ttl in [(ax1,"MMAD Trend (um)"),(ax2,"GSD Trend")]:
                ax.set_xticks(list(xi)); ax.set_xticklabels(names,rotation=20,ha="right",
                    fontsize=9,color="#c0d8f0")
                ax.set_title(ttl,color="#FFC600",fontsize=11,fontweight="bold")
                ax.tick_params(colors="#7090b0"); ax.spines[:].set_color("#2a4060")
                ax.grid(True,axis="y",color="#1a3050",ls="--",alpha=0.5)
            fig2.tight_layout()
            cv2=FigureCanvasTkAgg(fig2,master=scroll); cv2.draw()
            cv2.get_tk_widget().pack(fill="x",pady=(4,0))
        ref_sd=next((sd for sd in self.all_series if sd["is_ref"]),None)
        params_list=[("mmad","MMAD(um)"),("gsd","GSD"),("fpd","FPD(mg)"),("fpf","FPF(%)"),
                     ("slope","Slope"),("intercept","Intercept"),("r2","R2")]
        tf=ctk.CTkFrame(scroll,fg_color="#111827",corner_radius=6)
        tf.pack(fill="x",padx=8,pady=8)
        cw=[138]+[104]*len(self.all_series)
        hdrs=["Parametre"]+[sd["name"]+(" *" if sd["is_ref"] else "") for sd in self.all_series]
        hrow=ctk.CTkFrame(tf,fg_color="#1F4E79",corner_radius=0); hrow.pack(fill="x")
        for h,w in zip(hdrs,cw):
            ctk.CTkLabel(hrow,text=h,width=w,font=BF,text_color="white",
                anchor="center").pack(side="left",padx=1,pady=2)
        for i,(key,lbl) in enumerate(params_list):
            bg="#1a1a2a" if i%2==0 else "transparent"
            dr=ctk.CTkFrame(tf,fg_color=bg,corner_radius=0); dr.pack(fill="x")
            ik=key in("mmad","gsd","fpd","fpf")
            ctk.CTkLabel(dr,text=lbl,width=cw[0],font=BF if ik else NF,
                text_color="#FFC600" if ik else "#c0d0e0",anchor="w").pack(side="left",padx=(6,1),pady=2)
            rv=None
            if ref_sd and ref_sd["avg"] and key in ref_sd["avg"]["params"]:
                rv=ref_sd["avg"]["params"][key][0]
            for sd in self.all_series:
                if not sd["avg"] or key not in sd["avg"]["params"]:
                    ctk.CTkLabel(dr,text="-",width=104,font=NF,text_color="#888",
                        anchor="center").pack(side="left",padx=1,pady=2); continue
                val=sd["avg"]["params"][key][0]; txt=f"{val:.4f}"; clr="#d0e8ff"
                if rv and not sd["is_ref"] and rv>0:
                    diff=(val-rv)/rv*100; txt+=f"\n({diff:+.1f}%)"
                    clr="#90ee90" if abs(diff)<10 else "#FFB060" if abs(diff)<20 else "#FF6060"
                ctk.CTkLabel(dr,text=txt,width=104,font=ctk.CTkFont(size=12),
                    text_color=clr,anchor="center").pack(side="left",padx=1,pady=2)
        if ref_sd and ref_sd["avg"]:
            ctk.CTkLabel(scroll,text=f"  {self.T['f2_label']}",font=HF,
                text_color="#FFD700",anchor="w").pack(fill="x",padx=8,pady=(12,4))
            rm=ref_sd["avg"]["avg_masses"]
            lim_map={"ema":20,"fda":15,"usp":25}
            try: pct=lim_map.get(self.limit_var.get()) or float(self.custom_pct_var.get())
            except: pct=20
            for sd in self.all_series:
                if sd["is_ref"] or not sd["avg"]: continue
                f2=calc_f2(rm,sd["avg"]["avg_masses"],co)
                if f2 is None: continue
                pf2=self.T["f2_pass"] if f2>=50 else self.T["f2_fail"]
                clr="#90ee90" if f2>=50 else "#FF6060"
                ff=ctk.CTkFrame(scroll,fg_color="#1c2336",corner_radius=6)
                ff.pack(fill="x",padx=12,pady=2)
                ctk.CTkLabel(ff,text=f"  {sd['name']}  f2 = {f2:.1f}   {pf2}",
                    font=ctk.CTkFont(size=13,weight="bold"),text_color=clr).pack(anchor="w",padx=8,pady=6)

    def _toggle_lang(self):
        old_T=self.T; self.lang="EN" if self.lang=="TR" else "TR"; self.T=L[self.lang]
        self.lbl_title.configure(text=self.T["title"])
        self.lbl_sub.configure(text=self.T["subtitle"])
        self.btn_lang.configure(text=self.T["lang_btn"])
        self.btn_add_s.configure(text=self.T["add_series"])
        self.btn_del_s.configure(text=self.T["del_series"])
        self.btn_calc.configure(text=self.T["calculate"])
        self.btn_clr.configure(text=self.T["clear"])
        self.btn_pdf.configure(text=self.T["export_pdf"])
        self._refresh_cutoffs()
        for sw in self.series_widgets:
            sw["paste_btn"].configure(text=self.T["paste_btn"])
            if sw.get("ref_check"): sw["ref_check"].configure(text=self.T["ref_check"])
        for k in ["tab_results","tab_plot","tab_dist","tab_summary","tab_compare"]:
            try: self.tabs.rename(old_T[k],self.T[k])
            except: pass

    def _clear(self):
        for sw in self.series_widgets:
            for ri in range(RUNS_PER_SERIES):
                for s in sw["runs"][ri]: sw["runs"][ri][s].set("0.000")
        self.all_series=[]
        for tf in [self.rf,self.pf,self.df,self.sf,self.cf]:
            for w in tf.winfo_children(): w.destroy()

    def _export_pdf(self):
        if not self.all_series:
            messagebox.showwarning("","Oncelikle hesaplama yapiniz."); return
        path=filedialog.asksaveasfilename(defaultextension=".pdf",
            filetypes=[("PDF","*.pdf")],
            initialfile=f"NGI_{datetime.now().strftime('%Y%m%d_%H%M')}.pdf")
        if not path: return
        meta={"product":self.e_product.get(),"batch":self.e_batch.get(),
              "operator":self.e_operator.get(),"date":self.e_date.get()}
        lm={"ema":20,"fda":15,"usp":25}
        try: pct=lm.get(self.limit_var.get()) or float(self.custom_pct_var.get())
        except: pct=20
        try:
            make_pdf_multi(path,self.all_series,meta,int(self.var_flow.get()),self.T,pct)
            messagebox.showinfo("",f"PDF kaydedildi:\n{path}")
        except Exception as ex:
            messagebox.showerror("PDF Hatasi",str(ex))

def make_pdf_multi(path,all_series,meta,flow,T,limit_pct=20):
    from reportlab.lib.pagesizes import A4
    from reportlab.lib import colors
    from reportlab.lib.units import cm
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib.enums import TA_CENTER,TA_LEFT
    from reportlab.platypus import (SimpleDocTemplate,Table,TableStyle,
                                    Paragraph,Spacer,HRFlowable,PageBreak)
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    import os
    fn="Helvetica"; fb="Helvetica-Bold"
    for nm,ff in [("DejaVu","DejaVuSans.ttf"),("DejaVuB","DejaVuSans-Bold.ttf")]:
        fp=resource_path(ff)
        if os.path.exists(fp):
            try: pdfmetrics.registerFont(TTFont(nm,fp)); fn=nm; fb="DejaVuB"
            except: pass
            break
    CN=colors.HexColor("#002D62"); CB=colors.HexColor("#1F4E79")
    CG=colors.HexColor("#F2F2F2"); CD=colors.HexColor("#404040")
    CW=colors.white; CGREEN=colors.HexColor("#E2EFDA"); CRED=colors.HexColor("#FFCCCC")
    W,H=A4; BW=W-3*cm
    doc=SimpleDocTemplate(path,pagesize=A4,leftMargin=1.5*cm,rightMargin=1.5*cm,
        topMargin=1.5*cm,bottomMargin=1.5*cm)
    def s(sz,bold=False,color=colors.black,align=TA_LEFT,sp=2):
        return ParagraphStyle("x",fontName=fb if bold else fn,fontSize=sz,
            textColor=color,alignment=align,spaceBefore=sp,spaceAfter=sp,leading=sz+3)
    sH1=s(14,True,CW,TA_CENTER); sH2=s(9,False,colors.HexColor("#CCCCCC"),TA_CENTER)
    sSec=s(9,True,CW); sN=s(8); sB=s(8,True,CD); sHL=s(9,True,CN)
    def ts(hbg=CN):
        return TableStyle([("BACKGROUND",(0,0),(-1,0),hbg),("TEXTCOLOR",(0,0),(-1,0),CW),
            ("FONTNAME",(0,0),(-1,0),fb),("FONTSIZE",(0,0),(-1,-1),7.5),
            ("ALIGN",(0,0),(-1,-1),"CENTER"),("VALIGN",(0,0),(-1,-1),"MIDDLE"),
            ("GRID",(0,0),(-1,-1),0.4,colors.HexColor("#AAAAAA")),
            ("ROWBACKGROUNDS",(0,1),(-1,-1),[CW,CG]),
            ("TOPPADDING",(0,0),(-1,-1),2),("BOTTOMPADDING",(0,0),(-1,-1),2)])
    story=[]
    t=Table([[Paragraph("NGI Kaskadit Impaktor Analiz Araci",sH1)]],colWidths=[BW])
    t.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,-1),CN),
        ("TOPPADDING",(0,0),(-1,-1),8),("BOTTOMPADDING",(0,0),(-1,-1),8)]))
    story.append(t)
    t2=Table([[Paragraph("Ph.Eur 2.9.18 / USP &lt;601&gt;",sH2)]],colWidths=[BW])
    t2.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,-1),CB),
        ("TOPPADDING",(0,0),(-1,-1),3),("BOTTOMPADDING",(0,0),(-1,-1),3)]))
    story.append(t2); story.append(Spacer(1,0.3*cm))
    mt=Table([
        [Paragraph("Urun",sB),Paragraph(meta.get("product",""),sN),
         Paragraph("Lot",sB),Paragraph(meta.get("batch",""),sN),
         Paragraph("Flow",sB),Paragraph(f"{flow} L/min",sB)],
        [Paragraph("Analist",sB),Paragraph(meta.get("operator",""),sN),
         Paragraph("Tarih",sB),Paragraph(meta.get("date",""),sN),"",""],
    ],colWidths=[1.8*cm,4*cm,1.8*cm,4*cm,1.8*cm,3*cm])
    mt.setStyle(TableStyle([("GRID",(0,0),(-1,-1),0.4,colors.HexColor("#CCCCCC")),
        ("BACKGROUND",(0,0),(0,-1),CG),("BACKGROUND",(2,0),(2,-1),CG),
        ("BACKGROUND",(4,0),(4,-1),CG),
        ("TOPPADDING",(0,0),(-1,-1),2),("BOTTOMPADDING",(0,0),(-1,-1),2)]))
    story.append(mt); story.append(Spacer(1,0.3*cm))
    co=NGI_CUTOFFS[flow]
    vis=[s2 for s2 in ["S2","S3","S4","S5","S6","S7","MOC"] if co.get(s2,999)<900]
    cw2=BW/(len(vis)+1)
    ct=Table([[Paragraph("D50",sB)]+[Paragraph(s2,sB) for s2 in vis],
              [Paragraph(f"{flow}L",sN)]+[Paragraph(f"{co[s2]:.2f}",sN) for s2 in vis]],
             colWidths=[cw2]*(len(vis)+1)); ct.setStyle(ts(CB))
    story.append(ct); story.append(Spacer(1,0.4*cm))
    for sd in all_series:
        rt=" [REFERANS]" if sd["is_ref"] else ""
        sh=Table([[Paragraph(f"Seri: {sd['name']}{rt}",sSec)]],colWidths=[BW])
        sh.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,-1),CN),
            ("TOPPADDING",(0,0),(-1,-1),4),("BOTTOMPADDING",(0,0),(-1,-1),4)]))
        story.append(sh); story.append(Spacer(1,0.2*cm))
        for run in sd["runs"]:
            rh=Table([[Paragraph(f"Run {run['run_no']}",sSec)]],colWidths=[BW])
            rh.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,-1),CB),
                ("TOPPADDING",(0,0),(-1,-1),3),("BOTTOMPADDING",(0,0),(-1,-1),3)]))
            story.append(rh)
            if "error" in run:
                story.append(Paragraph(f"Yetersiz veri (n={run.get('n',0)})",sN))
                story.append(Spacer(1,0.2*cm)); continue
            ch=[Paragraph(x,sB) for x in ["Stage","D50","Mass","CumMass","Cum%","Valid","Probit z"]]
            cr=[ch]; vs={v["stage"] for v in run["valid"]}; cm2=0.0
            tss=ts(CD)
            for ri,row in enumerate(run["cum_data"]):
                cm2+=row["mass"]; pz=""
                if 0<row["u_pct"]<100:
                    try: pz=f"{norm.ppf(row['u_pct']/100):.4f}"
                    except: pass
                iv=row["stage"] in vs
                cr.append([Paragraph(row["stage"],sB if iv else sN),
                    Paragraph(f"{row['d50']:.3f}" if row['d50']<900 else "-",sN),
                    Paragraph(f"{row['mass']:.4f}",sN),Paragraph(f"{cm2:.4f}",sN),
                    Paragraph(f"{row['u_pct']:.3f}",sN),
                    Paragraph("Yes" if iv else "",sB if iv else sN),Paragraph(pz,sN)])
                if iv: tss.add("BACKGROUND",(0,ri+1),(-1,ri+1),CGREEN); tss.add("FONTNAME",(0,ri+1),(-1,ri+1),fb)
            cw3=BW/7; ct2=Table(cr,colWidths=[cw3]*7); ct2.setStyle(tss)
            story.append(ct2); story.append(Spacer(1,0.15*cm))
            dh=[Paragraph(x,sB) for x in ["Metered","Delivered","FPD","FPF%","MMAD","GSD","R2","n","Int","Slope"]]
            dr=[Paragraph(f"{run.get(k,0):.4f}" if k!="n" else str(run.get(k,0)),
                sHL if k in("fpd","fpf","mmad","gsd") else sN)
                for k in("metered","delivered","fpd","fpf","mmad","gsd","r2","n","intercept","slope")]
            cw4=BW/10; dt=Table([dh,dr],colWidths=[cw4]*10); dt.setStyle(ts(CD))
            story.append(dt); story.append(Spacer(1,0.3*cm))
        if sd.get("avg"):
            ah=Table([[Paragraph(f"Seri Ozeti: {sd['name']}",sSec)]],colWidths=[BW])
            ah.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,-1),colors.HexColor("#003580")),
                ("TOPPADDING",(0,0),(-1,-1),3),("BOTTOMPADDING",(0,0),(-1,-1),3)]))
            story.append(ah); pr=sd["avg"]["params"]
            ks=[("fpd","FPD(mg)"),("fpf","FPF%"),("mmad","MMAD"),("gsd","GSD"),
                ("slope","Slope"),("intercept","Int"),("r2","R2")]
            av_h=[Paragraph(l,sB) for _,l in ks]+[Paragraph("n",sB)]
            av_v=[]
            for k,_ in ks:
                if k in pr: m2,s2,rsd2=pr[k]; av_v.append(Paragraph(f"{m2:.4f}\n+/-{s2:.4f}\n{rsd2:.1f}%",sN))
                else: av_v.append(Paragraph("-",sN))
            av_v.append(Paragraph(str(sd["avg"]["n_valid"]),sN))
            cw5=BW/(len(ks)+1); at=Table([av_h,av_v],colWidths=[cw5]*(len(ks)+1)); at.setStyle(ts(CD))
            story.append(at)
        story.append(Spacer(1,0.5*cm))
    ref_sd=next((sd for sd in all_series if sd["is_ref"]),None)
    if ref_sd and ref_sd["avg"]:
        story.append(PageBreak())
        ph=Table([[Paragraph(f"REFERANS KARSILASTIRMA — {ref_sd['name']}  |  Limit: +/-{limit_pct:.0f}%",sSec)]],colWidths=[BW])
        ph.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,-1),colors.HexColor("#4a3000")),
            ("TOPPADDING",(0,0),(-1,-1),6),("BOTTOMPADDING",(0,0),(-1,-1),6)]))
        story.append(ph); story.append(Spacer(1,0.3*cm))
        rm=ref_sd["avg"]["avg_masses"]; vis2=[s2 for s2 in GRAPH_STAGES if co.get(s2,999)<900]
        for sd in all_series:
            if sd["is_ref"] or not sd["avg"]: continue
            f2=calc_f2(rm,sd["avg"]["avg_masses"],co)
            f2t=f"f2={f2:.1f} ({'Benzer>=50' if f2>=50 else 'Farkli<50'})" if f2 else "-"
            f2c=colors.HexColor("#006600") if (f2 and f2>=50) else colors.HexColor("#CC0000")
            story.append(Paragraph(f"{sd['name']}   {f2t}",
                ParagraphStyle("f2",fontName=fb,fontSize=11,textColor=f2c)))
            ch2=[Paragraph(x,sB) for x in ["Stage","Ref","Test","Fark%","Limit","Sonuc"]]
            cr2=[ch2]; tss2=ts(CD)
            for ri2,s2 in enumerate(vis2):
                rv=rm.get(s2,0); tv=sd["avg"]["avg_masses"].get(s2,0)
                diff=(tv-rv)/rv*100 if rv else 0; il=abs(diff)<=limit_pct
                cr2.append([Paragraph(s2,sB),Paragraph(f"{rv:.4f}",sN),
                    Paragraph(f"{tv:.4f}",sN),Paragraph(f"{diff:+.2f}%",sN),
                    Paragraph(f"+/-{limit_pct:.0f}%",sN),
                    Paragraph("OK" if il else "FAIL",
                        ParagraphStyle("ok",fontName=fb,fontSize=8,
                            textColor=colors.green if il else colors.red))])
                if not il: tss2.add("BACKGROUND",(0,ri2+1),(-1,ri2+1),CRED)
            cw6=BW/6; ct3=Table(cr2,colWidths=[cw6]*6); ct3.setStyle(tss2)
            story.append(ct3); story.append(Spacer(1,0.3*cm))
    story.append(HRFlowable(width="100%",thickness=0.5,color=colors.HexColor("#888888")))
    story.append(Paragraph(f"NGI Analysis Tool v5  |  {datetime.now().strftime('%d.%m.%Y %H:%M')}",
        ParagraphStyle("ft",fontName=fn,fontSize=7,
            textColor=colors.HexColor("#888888"),alignment=TA_CENTER)))
    doc.build(story)

if __name__=="__main__":
    app=NGIApp(); app.mainloop()