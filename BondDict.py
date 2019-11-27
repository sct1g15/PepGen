aa321 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
         'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

aa123 = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
         'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN',
         'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP',
         'A': 'ALA', 'V': 'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}

Lib_class = {"NonPGIV": ("CYS", "ASP", "SER", "GLN", "LYS",
                         "THR", "PHE", "ASN", "HIS", "LEU",
                         "ARG", "TRP", "ALA", "GLU", "TYR", "MET"),
             "Ile_Val": ("ILE", "VAL"),
             "Pro": "PRO",
             "Gly": "GLY"}
x_class = {"_nonxpro": ("CYS", "ASP", "SER", "GLN", "LYS",
                        "THR", "PHE", "ASN", "HIS", "LEU",
                        "ARG", "TRP", "ALA", "GLU", "TYR",
                        "MET", "ILE", "Val", "GLY"),
           "_xpro": "Pro"}

bond_length = ()




# def shorten(x):
#     if len(x) % 3 != 0:
#         raise ValueError('Input length should be a multiple of three')
#
#     y = ''
#     for i in range(len(x)/3):
#             y += d[x[3*i:3*i+3]]
#     return y