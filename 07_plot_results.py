"""
07_plot_results.py
Mereproduksi preprocessing pipeline & menghasilkan grafik statistik lengkap:
 1. Distribusi pIC50 (Histogram)
 2. Eksperimental vs Prediksi (R2 Plot)
 3. Residual Plot
 4. Williams Plot (Applicability Domain)
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn.model_selection import train_test_split, cross_val_predict, KFold
from sklearn.feature_selection import VarianceThreshold, RFE
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

# ─── LOAD & PREPROCESSING (Replikasi pipeline 04) ────────────────────────────
print("=== Mereplikasi preprocessing pipeline ===")
df = pd.read_csv('h1_antagonists_descriptors.csv').dropna()
print(f"  Data dimuat: {len(df)} baris")

X = df.drop(columns=['canonical_smiles', 'pIC50'], errors='ignore')
y = df['pIC50']

# Step 1: Variance Threshold
vt = VarianceThreshold(threshold=0.01)
X_vt = vt.fit_transform(X)
X = pd.DataFrame(X_vt, columns=X.columns[vt.get_support()])

# Step 2: Hapus korelasi tinggi (>0.85)
corr_matrix = X.corr().abs()
upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
to_drop = [col for col in upper.columns if any(upper[col] > 0.85)]
X = X.drop(columns=to_drop)
print(f"  Fitur tersisa setelah seleksi: {X.shape[1]}")

# Step 3: Train-Test Split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Step 4: Standarisasi
scaler = StandardScaler()
X_train_scaled = pd.DataFrame(scaler.fit_transform(X_train), columns=X_train.columns)
X_test_scaled  = pd.DataFrame(scaler.transform(X_test),  columns=X_test.columns)

# Step 5: RFE - pilih 5 fitur terbaik
selector = RFE(LinearRegression(), n_features_to_select=5, step=1)
selector.fit(X_train_scaled, y_train)
selected_features = X_train.columns[selector.support_].tolist()
print(f"  Fitur terpilih: {selected_features}")

X_train_f = X_train_scaled[selected_features]
X_test_f  = X_test_scaled[selected_features]

# Step 6: Model MLR
model = LinearRegression()
model.fit(X_train_f, y_train)

# ─── HITUNG METRIK ───────────────────────────────────────────────────────────
y_train_pred = model.predict(X_train_f)
y_test_pred  = model.predict(X_test_f)

r2_train  = r2_score(y_train, y_train_pred)
rmse_train = np.sqrt(mean_squared_error(y_train, y_train_pred))

cv = KFold(n_splits=5, shuffle=True, random_state=42)
y_cv_pred = cross_val_predict(model, X_train_f, y_train, cv=cv)
q2_cv = r2_score(y_train, y_cv_pred)

r2_test   = r2_score(y_test, y_test_pred)
rmse_test = np.sqrt(mean_squared_error(y_test, y_test_pred))

print("\n=== HASIL STATISTIK MLR ===")
print(f"  R2 Training     : {r2_train:.3f}")
print(f"  RMSE Training   : {rmse_train:.3f}")
print(f"  Q2 (5-Fold CV)  : {q2_cv:.3f}")
print(f"  R2 Test (Ext.)  : {r2_test:.3f}")
print(f"  RMSE Test       : {rmse_test:.3f}")

# Persamaan MLR
coef = model.coef_
intercept = model.intercept_
eq_parts = [f"pIC50 = {intercept:.3f}"]
for feat, c in zip(selected_features, coef):
    sign = "+" if c >= 0 else "-"
    eq_parts.append(f"{sign} {abs(c):.3f}×[{feat}]")
equation = " ".join(eq_parts)
print(f"\n  Persamaan MLR:\n  {equation}")

# ─── GAMBAR 4 PANEL ──────────────────────────────────────────────────────────
BLUE   = '#2563EB'
GREEN  = '#16A34A'
RED    = '#DC2626'
ORANGE = '#EA580C'
GRAY   = '#6B7280'

fig = plt.figure(figsize=(16, 14))
fig.suptitle(
    "Hasil Statistik Model QSAR – Multiple Linear Regression (MLR)\n"
    "Antagonis Reseptor Histamin H1 (ChEMBL231)",
    fontsize=14, fontweight='bold', y=0.98
)
gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.40, wspace=0.35)

# ── Panel 1: Distribusi pIC50 ─────────────────────────────────────────────
ax1 = fig.add_subplot(gs[0, 0])
ax1.hist(y, bins=30, color=BLUE, edgecolor='white', alpha=0.85)
ax1.axvline(y.mean(), color=RED, linestyle='--', linewidth=1.8, label=f'Mean = {y.mean():.2f}')
ax1.set_xlabel('pIC50', fontsize=11)
ax1.set_ylabel('Jumlah Senyawa', fontsize=11)
ax1.set_title('Distribusi Nilai Aktivitas (pIC50)', fontsize=12, fontweight='bold')
ax1.legend(fontsize=10)
ax1.grid(True, linestyle=':', alpha=0.6)
# Statistik deskriptif di pojok
stats_text = (f"n = {len(y)}\n"
              f"Min = {y.min():.2f}\n"
              f"Max = {y.max():.2f}\n"
              f"SD  = {y.std():.2f}")
ax1.text(0.97, 0.97, stats_text, transform=ax1.transAxes,
         fontsize=8.5, va='top', ha='right',
         bbox=dict(boxstyle='round,pad=0.4', facecolor='lightyellow', alpha=0.8))

# ── Panel 2: Eksperimental vs Prediksi ───────────────────────────────────
ax2 = fig.add_subplot(gs[0, 1])
ax2.scatter(y_train, y_train_pred, alpha=0.6, color=BLUE,
            edgecolor='white', s=40, label=f'Training (n={len(y_train)})')
ax2.scatter(y_test, y_test_pred, alpha=0.85, color=ORANGE,
            edgecolor='white', s=55, marker='D', label=f'Test (n={len(y_test)})')
min_val = min(y.min(), y_train_pred.min(), y_test_pred.min()) - 0.3
max_val = max(y.max(), y_train_pred.max(), y_test_pred.max()) + 0.3
ax2.plot([min_val, max_val], [min_val, max_val], 'k--', lw=1.5, label='Ideal (y = x)')
ax2.set_xlabel('pIC50 Eksperimental', fontsize=11)
ax2.set_ylabel('pIC50 Prediksi', fontsize=11)
ax2.set_title('Eksperimental vs Prediksi (pIC50)', fontsize=12, fontweight='bold')
ax2.legend(fontsize=9)
ax2.grid(True, linestyle=':', alpha=0.6)
stats_text2 = (f"R² Train  = {r2_train:.3f}\n"
               f"Q² (CV)   = {q2_cv:.3f}\n"
               f"R² Test   = {r2_test:.3f}\n"
               f"RMSE Test = {rmse_test:.3f}")
ax2.text(0.03, 0.97, stats_text2, transform=ax2.transAxes,
         fontsize=8.5, va='top', ha='left',
         bbox=dict(boxstyle='round,pad=0.4', facecolor='lightyellow', alpha=0.8))

# ── Panel 3: Residual Plot ────────────────────────────────────────────────
ax3 = fig.add_subplot(gs[1, 0])
res_train = y_train.values - y_train_pred
res_test  = y_test.values  - y_test_pred
ax3.scatter(y_train_pred, res_train, alpha=0.55, color=BLUE,
            edgecolor='white', s=40, label='Training')
ax3.scatter(y_test_pred,  res_test,  alpha=0.85, color=ORANGE,
            edgecolor='white', s=55, marker='D', label='Test')
ax3.axhline(0,   color='black', lw=1.2)
ax3.axhline( 2,  color=RED, lw=1.2, linestyle='--', alpha=0.7, label='±2 Batas Residual')
ax3.axhline(-2,  color=RED, lw=1.2, linestyle='--', alpha=0.7)
ax3.set_xlabel('pIC50 Prediksi', fontsize=11)
ax3.set_ylabel('Residual (Eksperimental − Prediksi)', fontsize=11)
ax3.set_title('Residual Plot', fontsize=12, fontweight='bold')
ax3.legend(fontsize=9)
ax3.grid(True, linestyle=':', alpha=0.6)

# ── Panel 4: Williams Plot (Applicability Domain) ────────────────────────
ax4 = fig.add_subplot(gs[1, 1])
X_all_scaled = pd.DataFrame(scaler.transform(X), columns=X.columns)
X_all_f = X_all_scaled[selected_features].values
H = X_all_f @ np.linalg.pinv(X_all_f.T @ X_all_f) @ X_all_f.T
hat = np.diag(H)

y_all_pred = model.predict(X_all_f)
residuals_all = y.values - y_all_pred
s = np.std(residuals_all, ddof=X_all_f.shape[1])
std_res = residuals_all / s if s > 0 else residuals_all

h_star = 3 * (X_all_f.shape[1] + 1) / X_all_f.shape[0]

# Warna berdasarkan zona
inside  = (hat <= h_star) & (np.abs(std_res) <= 3)
outside = ~inside

ax4.scatter(hat[inside],  std_res[inside],  alpha=0.55, color=BLUE,
            edgecolor='white', s=30, label='Dalam AD')
ax4.scatter(hat[outside], std_res[outside], alpha=0.85, color=RED,
            edgecolor='white', s=50, marker='^', label='Outlier / Luar AD')
ax4.axvline(h_star, color=GREEN, lw=1.8, linestyle='--',
            label=f'h* = {h_star:.3f}')
ax4.axhline( 3, color='gray', lw=1.2, linestyle='--')
ax4.axhline(-3, color='gray', lw=1.2, linestyle='--')
ax4.set_xlabel('Leverage (h)', fontsize=11)
ax4.set_ylabel('Residual Terstandarisasi (σ)', fontsize=11)
ax4.set_title('Williams Plot – Applicability Domain', fontsize=12, fontweight='bold')
ax4.legend(fontsize=9)
ax4.grid(True, linestyle=':', alpha=0.6)
pct_inside = inside.sum() / len(inside) * 100
ax4.text(0.97, 0.97, f"{pct_inside:.1f}% senyawa\ndalam AD",
         transform=ax4.transAxes, fontsize=8.5, va='top', ha='right',
         bbox=dict(boxstyle='round,pad=0.4', facecolor='lightgreen', alpha=0.7))

plt.savefig('hasil_statistik_qsar.png', dpi=300, bbox_inches='tight')
print("\nGrafik berhasil disimpan sebagai 'hasil_statistik_qsar.png'")
plt.close()
