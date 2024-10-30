import numpy as np


etas = np.arange(0.05,1.5,0.1)
# etas = [1]

labels = ['box', 'cutoff', 'Nparticles', 'Vo', 'eta', 'sigma', 'R',
       'epsilon', 'aligstr', 'T', 'Dt']

def generate_script(eta):
    script = f'''box 80
cutoff 0.5
Nparticles 6400
Vo 1
eta {eta:.2f}
sigma 0.5
R 1
epsilon 0.1
aligstr 1
T 10000
Dt 0.1
'''

    return script


for eta in etas:
    file = f'params_eta_{eta:.2f}.dat'
    with open(file,'w') as i:
        script = generate_script(eta)
        i.write(script)