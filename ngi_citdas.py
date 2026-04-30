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

# CITDAS ve USP/Ph.Eur uyumlu cutoff değerleri
NGI_CUTOFFS = {
    15: {"S1": 14.10, "S2": 8.61, "S3": 5.39, "S4": 3.30, "S5": 2.08, "S6": 1.36, "S7": 0.98, "MOC": 0.70},
    30: {"S1": 11.76, "S2": 6.40, "S3": 3.99, "S4": 2.30, "S5": 1.36, "S6": 0.83, "S7": 0.54, "MOC": 0.34},
    60: {"S1": 8.06,  "S2": 4.46, "S3": 2.82, "S4": 1.66, "S5": 0.94, "S6": 0.55, "S7": 0.34, "MOC": 0.14},
    90: {"S1": 6.48,  "S2": 3.61, "S3": 2.27, "S4": 1.33, "S5": 0.76, "S6": 0.44, "S7": 0.26, "MOC": 0.10},
}

ISM_STAGES = ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "MOC"]
ALL_KEYS = ["Device", "Throat", "Presep"] + ISM_STAGES

def calc_run(masses, flow):
    """CITDAS Mantığına Göre Hesaplama Çekirdeği"""
    co = NGI_CUTOFFS.get(flow, NGI_CUTOFFS[60])
    
    # ISM (Impactor Sized Mass) hesaplama
    ism = sum(masses.get(s, 0) for s in ISM_STAGES)
    metered = sum(masses.get(s, 0) for s in ALL_KEYS)
    
    if ism <= 0:
        return {"error": "no_data", "metered": metered, "delivered": ism}

    # Kümülatif Yüzde (Undersize) Hesaplama
    # CITDAS: Bir aşamanın altındaki toplam kütle / ISM
    cum_data = []
    for s in ISM_STAGES:
        cutoff = co.get(s, 0)
        # Bu aşamadan daha küçük olan aşamaların toplamı (Undersize)
        mass_below = sum(masses.get(x, 0) for x in ISM_STAGES if co.get(x, 0) < cutoff)
        u_pct = (mass_below / ism) * 100
        cum_data.append({"stage": s, "d50": cutoff, "mass": masses.get(s, 0), "u_pct": u_pct})

    # Regresyon Noktası Seçimi: %16 < u_pct < %84 (CITDAS n Değeri)
    valid = [r for r in cum_data if 16.0 <= r["u_pct"] <= 84.0 and r["d50"] > 0]
    
    if len(valid) < 2:
        return {"error": "insufficient", "n": len(valid), "metered": metered, "delivered": ism, "cum_data": cum_data}

    # Log-Probit Regresyonu (y = mx + c)
    # y = Probit + 5 (CITDAS Standardı)
    x = np.array([math.log10(v["d50"]) for v in valid])
    y = np.array([norm.ppf(v["u_pct"] / 100) + 5 for v in valid])
    
    # Lineer Regresyon
    slope, intercept = np.polyfit(x, y, 1)
    
    # R^2 Hesaplama
    y_pred = slope * x + intercept
    r2 = 1 - (np.sum((y - y_pred)**2) / np.sum((y - np.mean(y))**2))

    # MMAD Hesaplama: y=5 olduğu nokta (z=0)
    # 5 = slope * log10(MMAD) + intercept
    mmad = 10**((5 - intercept) / slope)

    # GSD Hesaplama: 10^(1/slope)[cite: 1]
    gsd = 10**(1 / slope)

    # FPD (5 um altı doz) Hesaplama
    # Regresyon doğrusunda x = log10(5) koyulur
    log5 = math.log10(5.0)
    y_5um = slope * log5 + intercept
    z_5um = y_5um - 5
    fpd_pct = norm.cdf(z_5um) # Kümülatif yüzde (0-1 arası)
    fpd = fpd_pct * ism
    fpf = (fpd / metered) * 100 if metered > 0 else 0

    return {
        "metered": metered, "delivered": ism, "n": len(valid),
        "mmad": mmad, "gsd": gsd, "fpd": fpd, "fpf": fpf,
        "slope": slope, "intercept": intercept, "r2": r2,
        "cum_data": cum_data, "x_reg": x, "y_reg": y, "valid": valid
    }

# --- GUI Kısmı (Basitleştirilmiş Test Arayüzü) ---
class NGIApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("CITDAS Validated NGI Tool")
        self.geometry("1000x700")
        
        # Sol Panel (Giriş)
        self.left = ctk.CTkFrame(self, width=300)
        self.left.pack(side="left", fill="y", padx=10, pady=10)
        
        ctk.CTkLabel(self.left, text="Akış Hızı (L/min):").pack(pady=5)
        self.flow_var = ctk.StringVar(value="15")
        ctk.CTkOptionMenu(self.left, values=["15", "30", "60", "90"], variable=self.flow_var).pack()
        
        self.entries = {}
        for s in ALL_KEYS:
            f = ctk.CTkFrame(self.left)
            f.pack(fill="x", padx=5, pady=2)
            ctk.CTkLabel(f, text=s, width=80).pack(side="left")
            e = ctk.CTkEntry(f, width=100)
            e.insert(0, "0.00")
            e.pack(side="right")
            self.entries[s] = e
            
        ctk.CTkButton(self.left, text="HESAPLA", command=self.run_calc, fg_color="green").pack(pady=20)
        
        # Sağ Panel (Sonuç)
        self.right = ctk.CTkTextbox(self, width=600)
        self.right.pack(side="right", fill="both", expand=True, padx=10, pady=10)

    def run_calc(self):
        m = {s: float(self.entries[s].get()) for s in ALL_KEYS}
        flow = int(self.flow_var.get())
        res = calc_run(m, flow)
        
        self.right.delete("1.0", "end")
        if "error" in res:
            self.right.insert("end", f"HATA: {res['error']} (Nokta sayısı: {res.get('n', 0)})")
            return
            
        txt = f"--- HESAPLAMA SONUÇLARI ({flow} L/min) ---\n\n"
        txt += f"MMAD: {res['mmad']:.4f} µm\n"
        txt += f"GSD:  {res['gsd']:.4f}\n"
        txt += f"FPD:  {res['fpd']:.4f} µg\n"
        txt += f"FPF:  %{res['fpf']:.3f}\n"
        txt += f"R²:   {res['r2']:.4f}\n"
        txt += f"n:    {res['n']} (Regresyona giren nokta sayısı)\n"
        txt += f"Eğim (Slope): {res['slope']:.4f}\n"
        txt += f"Kesen (Intercept): {res['intercept']:.4f}\n\n"
        txt += "--- AŞAMA DETAYLARI (Undersize %) ---\n"
        for row in res['cum_data']:
            txt += f"{row['stage']}: Cutoff={row['d50']:.2f}, Kütle={row['mass']:.4f}, Cum%={row['u_pct']:.2f}\n"
            
        self.right.insert("end", txt)

if __name__ == "__main__":
    app = NGIApp()
    app.mainloop()