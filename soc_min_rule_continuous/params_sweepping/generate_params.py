import numpy as np
import subprocess


aligns = np.arange(3,4.5,0.5)



labels = ['box', 'cutoff', 'Nparticles', 'Vo', 'eta', 'sigma', 'R',
       'epsilon', 'aligstr', 'T', 'Dt']

def generate_script(alignstr):
    if alignstr.is_integer():
        alignstr = int(alignstr)
    script = f'''box 32
cutoff 0.5
Nparticles 1000
Vo 0.5
eta 0.9
sigma 0.5
R 1
epsilon 0.1
aligstr {alignstr}
T 10000
Dt 0.1
'''

    return script


for alignstr in aligns:
    file = f"params_align_{alignstr}.dat"
    with open(file,'w') as i:
        script = generate_script(alignstr)
        i.write(script)
    subprocess.Popen(["./soc_min_rule",file])
    print(f'align {alignstr:.1f} running')
