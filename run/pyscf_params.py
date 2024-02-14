Li_631g_star_star_trunc={
     # The units the geometry of the molecule is set up in
    'units':'Bohr',
     # The geometry of the molecule being investigated
    'atoms': 'Li 0 0 0; Li 0 0 6',
    # The type of basis used to generate the 1 and 2 electron integrals
    'bs' : '6-31g**',
    # How verbose do you want the PyScf output to be in your terminal?
    'verbosity' : 2,
    'symmetry' :True,
    'spin':0,
    'charge':0,
    'symmetry_subgroup' : 0, #0 is code for A1 point group
    'norb': 5,
    'nel':6
}

Li2_sto3g={
     # The units the geometry of the molecule is set up in
    'units':'Angstrom',
     # The geometry of the molecule being investigated
    'atoms': 'Li 0 0 0; Li 0 0 2.673',
    # The type of basis used to generate the 1 and 2 electron integrals
    'bs' : 'sto-3g',
    # How verbose do you want the PyScf output to be in your terminal?
    'verbosity' : 2,
    'symmetry' :True,
    'spin':0,
    'charge':0,
    'symmetry_subgroup' : 0, #0 is code for A1 point group
    'norb': 10,
    'nel':6
}

Li2_ccpvdz={
     # The units the geometry of the molecule is set up in
    'units':'Angstrom',
     # The geometry of the molecule being investigated
    'atoms': 'Li 0 0 0; Li 0 0 2.673',
    # The type of basis used to generate the 1 and 2 electron integrals
    'bs' : 'cc-pVDZ',
    # How verbose do you want the PyScf output to be in your terminal?
    'verbosity' : 2,
    'symmetry' :True,
    'spin':0,
    'charge':0,
    'symmetry_subgroup' : 0, #0 is code for A1 point group
    'norb': 28,
    'nel':6
}

Li_atom_ccpvdz={
     # The units the geometry of the molecule is set up in
    'units':'atom',
     # The geometry of the molecule being investigated
    'atoms': 'Li 0 0 0',
    # The type of basis used to generate the 1 and 2 electron integrals
    'bs' : 'cc-pVDZ',
    # How verbose do you want the PyScf output to be in your terminal?
    'verbosity' : 2,
    'symmetry' :True,
    'spin':1,
    'charge':0,
    'norb': 14,
    'nel':3
}

Be_atom_ccpvdz={
     # The units the geometry of the molecule is set up in
    'units':'atom',
     # The geometry of the molecule being investigated
    'atoms': 'Be 0 0 0',
    # The type of basis used to generate the 1 and 2 electron integrals
    'bs' : 'cc-pVDZ',
    # How verbose do you want the PyScf output to be in your terminal?
    'verbosity' : 2,
    'symmetry' :True,
    'spin':0,
    'charge':0,
    'norb': 14,
    'nel':4
}

N_atom_ccpvdz={
     # The units the geometry of the molecule is set up in
    'units':'atom',
     # The geometry of the molecule being investigated
    'atoms': 'N 0 0 0',
    # The type of basis used to generate the 1 and 2 electron integrals
    'bs' : 'cc-pVDZ',
    # How verbose do you want the PyScf output to be in your terminal?
    'verbosity' : 2,
    'symmetry' :True,
    'spin':1,
    'charge':0,
    'norb': 14,
    'nel':7
}

F_atom_ccpvdz={
     # The units the geometry of the molecule is set up in
    'units':'atom',
     # The geometry of the molecule being investigated
    'atoms': 'F 0 0 0',
    # The type of basis used to generate the 1 and 2 electron integrals
    'bs' : 'cc-pVDZ',
    # How verbose do you want the PyScf output to be in your terminal?
    'verbosity' : 2,
    'symmetry' :True,
    'spin':1,
    'charge':0,
    'norb': 14,
    'nel':9
}

Be2_ccpvdz={
     # The units the geometry of the molecule is set up in
    'units':'Angstrom',
     # The geometry of the molecule being investigated
    'atoms': 'Be 0 0 0; Be 0 0 2.45',
    # The type of basis used to generate the 1 and 2 electron integrals
    'bs' : 'cc-pVDZ',
    # How verbose do you want the PyScf output to be in your terminal?
    'verbosity' : 2,
    'symmetry' :True,
    'spin':0,
    'charge':0,
    'symmetry_subgroup' : 0, #0 is code for A1 point group
    'norb': 28,
    'nel':8
}

BH_631g_star_star_stretch={
     # The units the geometry of the molecule is set up in
    'units':'Angstrom',
     # The geometry of the molecule being investigated
    'atoms': 'B 0 0 0; H 0 0 4.0',
    # The type of basis used to generate the 1 and 2 electron integrals
    'bs' : '6-31g**',
    # How verbose do you want the PyScf output to be in your terminal?
    'verbosity' : 2,
    'symmetry' :True,
    'spin':0,
    'charge':0,
    'symmetry_subgroup' : 0, #0 is code for A1 point group
    'norb': 19,
    'nel':6
}

BH_631g_star_star_stretch={
     # The units the geometry of the molecule is set up in
    'units':'Angstrom',
     # The geometry of the molecule being investigated
    'atoms': 'B 0 0 0; H 0 0 1.2324',
    # The type of basis used to generate the 1 and 2 electron integrals
    'bs' : '6-31g**',
    # How verbose do you want the PyScf output to be in your terminal?
    'verbosity' : 2,
    'symmetry' :True,
    'spin':0,
    'charge':0,
    'symmetry_subgroup' : 0, #0 is code for A1 point group
    'norb':19,
    'nel':6
}

BH_321g_stretch={
     # The units the geometry of the molecule is set up in
    'units':'Angstrom',
     # The geometry of the molecule being investigated
    'atoms': 'B 0 0 0; H 0 0 1.2324',
    # The type of basis used to generate the 1 and 2 electron integrals
    'bs' : '3-21G', #'6-31g**',
    # How verbose do you want the PyScf output to be in your terminal?
    'verbosity' : 2,
    'symmetry' :True,
    'spin':0,
    'charge':0,
    'symmetry_subgroup' : 0, #0 is code for A1 point group
    'norb': 11,
    'nel':6
}

Li2_321g={
     # The units the geometry of the molecule is set up in
    'units':'Angstrom',
     # The geometry of the molecule being investigated
    'atoms': 'Li 0 0 0; Li 0 0 2.673',
    # The type of basis used to generate the 1 and 2 electron integrals
    'bs' : '3-21G', #'6-31g**',
    # How verbose do you want the PyScf output to be in your terminal?
    'verbosity' : 2,
    'symmetry' :True,
    'spin':0,
    'charge':0,
    'symmetry_subgroup' : 0, #0 is code for A1 point group
    'norb': 18,
    'nel':6
}


N2_ccpvdz={
     # The units the geometry of the molecule is set up in
    'units':'Angstrom',
     # The geometry of the molecule being investigated
    'atoms': 'N 0 0 0; N 0 0 1.094',
    # The type of basis used to generate the 1 and 2 electron integrals
    'bs' : 'cc-pVDZ',
    # How verbose do you want the PyScf output to be in your terminal?
    'verbosity' : 2,
    'symmetry' :True,
    'spin':0,
    'charge':0,
    'symmetry_subgroup' : 0, #0 is code for A1 point group
    'norb': 28,
    'nel':14
}

H2_ccpvdz={
     # The units the geometry of the molecule is set up in
    'units':'Angstrom',
     # The geometry of the molecule being investigated
    'atoms': 'H 0 0 0; H 0 0 0.741',
    # The type of basis used to generate the 1 and 2 electron integrals
    'bs' : 'cc-pVDZ',
    # How verbose do you want the PyScf output to be in your terminal?
    'verbosity' : 2,
    'symmetry' :True,
    'spin':0,
    'charge':0,
    # 'symmetry_subgroup' : 0, #0 is code for A1 point group
    'norb': 10,
    'nel':2
}