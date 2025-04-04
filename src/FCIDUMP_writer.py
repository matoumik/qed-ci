import numpy as np
from math import sqrt

def make_FCIDUMP(PFgen):
    print("Generating FCIDUMP FOR DMRG")
    
    # determine orbital range for the active space
    if(PFgen.n_act_orb == 0):
        act_start = 0;
        num_act = PFgen.nmo
        print("Full space dump")
    else:
        act_start = PFgen.ndocc-PFgen.n_act_el//2;
        act_end = act_start+PFgen.n_act_orb;
        num_act = PFgen.n_act_orb
    
    #dip_two_el = np.tensordot(PFgen.d_cmo,PFgen.d_cmo,axes=0)

    #dip_dot = PFgen.nuclear_dipole_moment[0]*PFgen.lambda_vector[0] + PFgen.nuclear_dipole_moment[1]*PFgen.lambda_vector[1] + PFgen.nuclear_dipole_moment[2]*PFgen.lambda_vector[2] 
    dip_dot = -PFgen.d_exp
    dip_one_el = PFgen.q_mo + PFgen.d_cmo*dip_dot

    #ints_two_el  = PFgen.eri_dump + dip_two_el

    ints_one_el  = PFgen.T_p_V_mo + dip_one_el
    omega = PFgen.omega


    dump = open("DMRGDUMP", "w")
    eldump = open("FCIDUMP_el", "w")
    phdump = open("FCIDUMP_ph", "w")
    intdump = open("FCIDUMP_int", "w")
    
    header = "&FCI NORB=" +str(num_act)+""",NELEC=0,MS2=0
ORBSYM=1,1,1,1,
ISYM=1,
/
"""

    eldump.write(header)
    phdump.write(header)
    intdump.write(header)
    
    dump.write(f' # Printing electron repulsion integrals in spatial MO basis\n')
    for i in range(num_act):
        for j in range(0,i+1):
            for k in range(num_act):
                for l in range(0,k+1):
                    teint = PFgen.eri_dump[act_start+i,act_start+j,act_start+k,act_start+l]
                    if np.abs(teint)>0.0:
                        dump.write(f'  {teint:30.20e} \t {i+1} {j+1} {k+1} {l+1}\n')
		    
                    #teint = ints_two_el[act_start+i,act_start+j,act_start+k,act_start+l]
                    teint += PFgen.d_cmo[act_start+i,act_start+j]*PFgen.d_cmo[act_start+k,act_start+l]
                    if np.abs(teint)>0.0:    
                     	eldump.write(f'  {teint:30.20e} \t {i+1} {j+1} {k+1} {l+1}\n')
                        

    dump.write(f' # Printing T + V integrals in spatial MO basis\n')
    
    for i in range(num_act):
        for j in range(num_act):
            #oeint = PFgen.T_p_V_mo[act_start+i,act_start+j]
            oeint = ints_one_el[act_start+i,act_start+j]
            for o in range(act_start):
                oeint+=2*PFgen.eri_dump[o,o,act_start+i,act_start+j]
                oeint-=PFgen.eri_dump[o,act_start+i,o,act_start+j]
                oeint+=2*PFgen.d_cmo[o,o]*PFgen.d_cmo[act_start+i,act_start+j]
                oeint-=PFgen.d_cmo[o,act_start+i]*PFgen.d_cmo[o,act_start+j]
            
            eldump.write(f'  {oeint:30.20e} \t {i+1} {j+1} {0} {0}\n')
            dump.write(f'  {oeint:30.20e} \t {i+1} {j+1} {0} {0}\n')

    dump.write(f' # Printing -1/2 \lambda \lambda q integrals in spatial MO basis\n')
    for i in range(num_act):
        for j in range(num_act):
            oeint = PFgen.q_mo[act_start+i,act_start+j]
            dump.write(f'  {oeint:30.20e} \t {i+1} {j+1} {0} {0}\n')

    dump.write(f' # Printing \lambda mu integrals in spatial MO basis\n')
    for i in range(num_act):
        for j in range(num_act):
            oeint = PFgen.d_cmo[act_start+i,act_start+j]
            dump.write(f'  {oeint:30.20e} \t {i+1} {j+1} {0} {0}\n')
            oeint *= -sqrt(omega/2)
            intdump.write(f'  {oeint:30.20e} \t {i+1} {j+1} {1} {1}\n')
            
    intdump.write(f'  {0.0:30.20e} \t {0} {0} {0} {0}\n')           


    dump.write(f' # Printing nuclear repulsion energy\n')
    Enuc = 0
    Enuc += PFgen.Enuc
    for o in range(act_start):
        Enuc += 2*PFgen.T_p_V_mo[o,o]
        Enuc += 2*dip_one_el[o,o]
        Enuc += PFgen.eri_dump[o,o,o,o]
        Enuc += PFgen.d_cmo[o,o]*PFgen.d_cmo[o,o]
        for p in range(o):
            Enuc+=4*PFgen.eri_dump[o,o,p,p]
            Enuc-=2*PFgen.eri_dump[o,p,o,p]
            Enuc+=4*PFgen.d_cmo[o,o]*PFgen.d_cmo[p,p]
            Enuc-=2*PFgen.d_cmo[o,p]*PFgen.d_cmo[o,p]
    dump.write(f'  {Enuc:30.20e} \t {0} {0} {0} {0}\n')
    #now add the dipole contributions to the energy
    Enuc += 0.5*dip_dot*dip_dot 
    eldump.write(f'  {Enuc:30.20e} \t {0} {0} {0} {0}\n')

    dump.write(f' # Printing photon energy\n')
    dump.write(f'  {PFgen.omega:30.20e} \t {0} {0} {0} {0}\n')
    phdump.write(f'  {PFgen.omega:30.20e} \t {1} {1} {1} {1}\n')

    
    dump.write(f' # Printing lambda vector\n')
    dump.write(f'  {PFgen.lambda_vector[0]:16.10e}, {PFgen.lambda_vector[1]:16.10e}, {PFgen.lambda_vector[2]:16.10e} \t {0} {0} {0} {0}\n')
    
    dump.write(f' # Printing RHF Electronic Dipole Moment\n')
    dump.write(f'  {PFgen.electronic_dipole_moment[0]:16.10e}, {PFgen.electronic_dipole_moment[1]:16.10e}, {PFgen.electronic_dipole_moment[2]:16.10e} \t {0} {0} {0} {0}\n')
    
    dump.write(f' # Printing Nuclear Dipole Moment\n')
    
    Dnuc = 0
    for o in range(act_start):
        #Dnuc += 2*PFgen.q_mo[o,o]
        Dnuc += 2*PFgen.d_cmo[o,o]
        
    Dnuc+=dip_dot
    
    Dnuc*=-sqrt(omega/2)
    
    phdump.write(f'  {Dnuc:16.10e} \t {1} {1} {0} {0}\n')
    phdump.write(f'  {0.0:16.10e} \t {0} {0} {0} {0}\n')
    
    lambda_norm = PFgen.lambda_vector[0]*PFgen.lambda_vector[0] + PFgen.lambda_vector[1]*PFgen.lambda_vector[1] + PFgen.lambda_vector[2]*PFgen.lambda_vector[2]
    #lambda_norm = math.sqrt(lambda_norm)
    #lambda_norm = 0;
    if(lambda_norm>1e-12):
        Dnuc_vec = [PFgen.nuclear_dipole_moment[0] + Dnuc * PFgen.lambda_vector[0]/lambda_norm, PFgen.nuclear_dipole_moment[1] + Dnuc * PFgen.lambda_vector[1]/lambda_norm, PFgen.nuclear_dipole_moment[2] + Dnuc * PFgen.lambda_vector[2]/lambda_norm ]
    else:
        Dnuc_vec = [PFgen.nuclear_dipole_moment[0], PFgen.nuclear_dipole_moment[1], PFgen.nuclear_dipole_moment[2]] 

    dump.write(f'  {Dnuc_vec[0]:16.10e}, {Dnuc_vec[1]:16.10e}, {Dnuc_vec[2]:16.10e} \t {0} {0} {0} {0}\n')
   
