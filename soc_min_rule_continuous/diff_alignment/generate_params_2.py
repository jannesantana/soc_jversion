import numpy as np
import subprocess


sizes = [1002,2000,5000,8000]
rho = 1



labels = ['box', 'cutoff', 'Nparticles', 'Vo', 'eta', 'sigma', 'R',
       'epsilon', 'aligstr', 'aligstrmin','T', 'Dt']

def generate_script(size):
    # if alignstr.is_integer():
    #     alignstr = int(alignstr)
    script = f'''box {int(np.sqrt(size/rho))}
cutoff 0.5
Nparticles {size}
Vo 0.5
eta 0.6
sigma 0.5
R 1
epsilon 0.1
aligstr 1
aligstrmin 3
T 20000
Dt 0.1
'''

    return script

def generate_bash_script(size):
    script = f'''\
#$ -N soc_min_rule_diff_align_size_{size}
#$ -q seq_medium
## Quantité de mémoire réservée
#$ -l m_mem_free=80G
#$ -cwd 
#$ -j y

./soc_diff_align params_size_{size}.dat
'''
    return script


for size in sizes:
    file = f"params_size_{size}.dat"
    with open(file,'w') as i:
        script = generate_script(size)
        i.write(script)
    bashname=f"job_{size}.sh"
    with open(bashname,'w') as j:
        sub_script = generate_bash_script(size)
        j.write(sub_script)
    subprocess.run(["qsub",bashname])
    print(f'size {size} running')
