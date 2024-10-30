import numpy as np
import subprocess 

etas = np.arange(0.05,1.5,0.1)

for eta in etas:
    subprocess.Popen(["./antimips",f"params_eta_{eta:.2f}.dat"])
    print(f'eta {eta:.2f} running')
