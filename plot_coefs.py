import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# df = pd.read_csv("adjoint_test/alm_p3_b2_eps1_interp/coefs.csv")
df = pd.read_csv("results/neko/alm/coefs.csv")
df.columns = df.columns.str.strip()

n = 20
print(f"Drag: {df['Drag'].iloc[-1]}")
print(f"Convergence delta (n={n}): {abs(df['Drag'].iloc[-1] - df['Drag'].iloc[-1 - n])}")

plt.figure()
plt.xlabel("Time step [-]")
plt.ylabel("Drag [-]")
plt.plot(df["t"], df["Drag"], label="Neko ALM (eps=b/8)")
# plt.plot([0, len(df["Drag"])], [0.00548 / 5.0, 0.00548 / 5.0], "r--", label=r"Matlab ALM (eps=b/8)")
plt.legend()
plt.title("Drag convergence")
plt.show()
