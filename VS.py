from time import process_time
from threading import Thread
from numpy import array, zeros, argsort
from scipy.spatial.distance import cdist
from rdkit import rdBase, DataStructs
from rdkit.Chem import SDMolSupplier, MolToSmiles, rdmolops, rdMolDescriptors, MACCSkeys
from rdkit.ML.Descriptors import MoleculeDescriptors
rdBase.DisableLog('rdApp.*') #rdApp.error means suppress the standard error output from rdkit 

from tkinter import Tk, Text, Radiobutton, Checkbutton, messagebox, StringVar, IntVar, Menu
from tkinter.messagebox import showinfo, showwarning
from tkinter.ttk import Label, Entry, Button, LabelFrame, Combobox
from tkinter.scrolledtext import ScrolledText
from tkinter.filedialog import askopenfilename, asksaveasfilename

def cal_Descriptor(mols, desc_List):
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(desc_List)
    desc_Set = list(map(calc.CalcDescriptors, mols))
    desc_Set = array(desc_Set)
    return desc_Set

def cal_FP(mols, fptype):
    def toArray(fp):
        fp_arr = zeros((1,))
        DataStructs.ConvertToNumpyArray(fp,fp_arr)
        return array(fp_arr)
    if fptype == 'RDKitFingerprint':
        fingerprint = map(lambda mol: rdmolops.RDKFingerprint(mol), mols)
    elif fptype == 'Morgan2':
        fingerprint = map(lambda mol: rdMolDescriptors.GetMorganFingerprintAsBitVect(mol,2), mols)
    elif fptype == 'MACCSkeys':
        fingerprint = map(lambda mol: MACCSkeys.GenMACCSKeys(mol), mols)
    elif fptype == 'AtomPairs':
        fingerprint = map(lambda mol: rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol), mols)
    elif fptype == 'TopoTorsion':
        fingerprint = map(lambda mol: rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(mol), mols)
    fingerprint1 = list(map(toArray, fingerprint))
    return array(fingerprint1)

def reset():
    workspace_path1.set('')
    workspace_path2.set('')
    for i in range(len(v)):
        v[i].set(1)
    v1.set(0)
    v2.set(0)
    cmb.current(0)
    text1.delete('1.0','end')
    neighbor.delete(0, 'end')
    neighbor.insert(0, 10)

def warn():
    showwarning('Warn', '!!!')

def load_mol():
    global mols1
    filename = askopenfilename(filetypes=[('sdf files','.sdf'),
                                          ('all files','.*')])
    if filename:
        workspace_path1.set(filename)
        text1.insert('insert', 'Load ' + str(filename) + '\n')
        text1.update()
        mols1 = [mol for mol in SDMolSupplier(filename) if mol is not None]
        text1.insert('insert', 'Load finished!\n')
        text1.update()
    else:
        showinfo('Note', 'Please load a .sdf file!')

def load_lib():
    global mols2
    filename = askopenfilename(filetypes=[('sdf files','.sdf'),
                                          ('all files','.*')])
    if filename:
        workspace_path2.set(filename)
        text1.insert('insert', 'Load ' + filename + '\n')
        text1.insert('insert', 'File is large, please be patient!!!\n')
        text1.update()
        mols2 = [mol for mol in SDMolSupplier(filename) if mol is not None]
        text1.insert('insert', 'Load finished!\n')
        text1.update()
    else:
        showinfo('Note', 'Please load a .sdf file!')
        
def submit():
    text1.insert('insert', 'Starting searching...\nPlease wait!\n')
    text1.update()
    start = process_time()
    # First select descriptor or fingerprint to calculate distance.
    if v2.get() == 0: #0 means using molecular descriptor
        desc = [descriptors[idx] for idx,var in enumerate(v) if var.get() == 1]
        prop_target = cal_Descriptor(mols1, desc)
        prop_ref = cal_Descriptor(mols2, desc)
    elif v2.get() == 1: #1 means using molecular fingerprint
        prop_target = cal_FP(mols1, fingerprints[v1.get()])
        prop_ref = cal_FP(mols2, fingerprints[v1.get()])

    #Calculate distance and then sort. 
    if cmb.get() == 'tanimoto':
        dists = cdist(prop_target, prop_ref, metric='jaccard')
    else:
        dists = cdist(prop_target, prop_ref, metric=cmb.get())
    idx = argsort(dists, axis=1)
    k = int(neighbor.get())
    idx_nearest = idx[:, :k]
    f = open('result.csv', 'w')
    f.write('SMILES,Cluster Center Index\n')
    #f.write('target,' + ','.join(['mol' + str(n) for n in range(k)]) + '\n')
    #for i in range(len(mols1)):
    #    nearest_mols = [MolToSmiles(mols2[j]) for j in idx_nearest[i]]
    #    f.write(MolToSmiles(mols1[i]) + ',' + ','.join(nearest_mols) + '\n')
    #for i in range(len(mols1)):
    #   nearest_mols = [MolToSmiles(mols2[j]) for j in idx_nearest[i]]
    #   for j in nearest_mols:
    #       f.write(j + ',' + MolToSmiles(mols1[i]) + ',' + str(i + 1) + '\n')
    for i in range(len(mols1)):
        f.write(MolToSmiles(mols1[i]) + ',' + str(i + 1) + '\n')
        nearest_mols = [MolToSmiles(mols2[j]) for j in idx_nearest[i]]
        for j in nearest_mols:
            f.write(j + ',' + '' + '\n')
    
    f.close()
    
    end = process_time()
    text1.insert('insert', 'Escape time:{}s'.format(end-start) + '\n')
    text1.insert('insert', 'Search finished! Please check the result.csv in your curreny workspace!\n')
    text1.update()

def thread_it(func, *args):
    t = Thread(target=func, args=args) 
    t.setDaemon(True) 
    t.start()
    
def show():
    showinfo('About', 'A simple virtual screening tool.')
    
root = Tk()
root.title('Similarity-Based Virtual Screening')
screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()
win_width, win_height = 750, 600
x = (screen_width - win_width) / 2
y = (screen_height - win_height) / 2
root.geometry('%dx%d+%d+%d' %(win_width, win_height, x, y))
root.resizable(0,0)

menubar = Menu(root)
menubar.add_command(label = "About", command = show)
root.config(menu=menubar)
#####################################
lb1 = LabelFrame(root, text='Input Data')
lb1.grid(row=0, column=0, columnspan=4, padx=10, pady=10)

workspace_path1 = StringVar(root)
mol = Label(lb1, text='Load Molecules:')
mol.grid(row=0, column=0)
input1 = Entry(lb1, textvariable=workspace_path1, width=78)
input1.grid(row=0, column=1)

workspace_path2 = StringVar(root)
lib = Label(lb1, text='Load Library:')
lib.grid(row=1, column=0)
input2 = Entry(lb1, textvariable=workspace_path2, width=78)
input2.grid(row=1, column=1)

load1 = Button(lb1, text='OpenMol',command=load_mol)
load1.grid(row=0,column=2, padx=3)
load2 = Button(lb1, text='OpenLib',command=load_lib)
load2.grid(row=1,column=2, padx=3)

##############################################
lb2 = LabelFrame(root, text='Molecule Property')
lb2.grid(row=1, column=0, columnspan=4, padx=10)

lb21 = LabelFrame(lb2, text='Molecule Descriptors')
lb21.grid(row=0, column=0, columnspan=2, padx=3, ipadx=10)

descriptors = ['MolWt', 'NumAliphaticCarbocycles', 'NumAliphaticHeterocycles', 
             'NumAromaticHeterocycles', 'NumAromaticRings', 'NumHAcceptors', 
             'NumHDonors', 'NumHeteroatoms', 'NumRotatableBonds', 
             'RingCount', 'MolLogP', 'TPSA']
v =  []
for idx, descriptor in enumerate(descriptors):
    v.append(IntVar())
    ckb = Checkbutton(lb21, text=descriptor, variable=v[-1])
    ckb.select()
    ckb.grid(row=idx // 2, column=idx % 2, sticky='W')

lb22 = LabelFrame(lb2, text='Molecule Fingerprints')
lb22.grid(row=0, column=2, columnspan=2, padx=3, ipadx=101, ipady=9)
fingerprints = ['MACCSkeys', 'RDKitFingerprint', 'Morgan2', 'AtomPairs', 'TopoTorsion']
v1 = IntVar()
v1.set(0)
for idx, fingerprint in enumerate(fingerprints):
    rdb = Radiobutton(lb22, text=fingerprint, variable=v1, value=str(idx), command=None)
    rdb.grid(row=idx, sticky='W')

######################################
lb3 = LabelFrame(root, text='Calculate Similarity by')
lb3.grid(row=2, column=0, columnspan=4, ipadx=273, pady=10)

types = ['Descriptor', 'Fingerprint']
v2 = IntVar()
v2.set(0)
for idx, type in enumerate(types):
    rdb = Radiobutton(lb3, text=type, variable=v2, value=str(idx), command=None)
    rdb.grid(row=0, column=idx, padx=2)

###########################################
lb4 = LabelFrame(root, text='Similarity Parameters')
lb4.grid(row=3, column=0, columnspan=4, ipadx=83)

Label(lb4, text='Metrics:').grid(row=0, column=0, padx=2, pady=3)

cmb = Combobox(lb4)
cmb.grid(row=0, column=1, padx=5, pady=3)
cmb['value'] = ('cosine', 'dice', 'euclidean', 'hamming', 'minkowski', 'tanimoto')
cmb.current(0)

Label(lb4, text='Nearest Neighbors Per Molecule:').grid(row=0, column=2, pady=3)
neighbor = Entry(lb4)
neighbor.insert(0, 10)
neighbor.grid(row=0, column=3, padx=5, pady=3)

Button(root, text='Submit', command=lambda: thread_it(submit)).grid(row=4, column=0, pady=10, ipadx=15)
Button(root, text="Kill", command=None).grid(row=4, column=1, pady=10, ipadx=15)
Button(root, text="Reset", command=reset).grid(row=4, column=2, pady=10, ipadx=15)
Button(root, text='Quit', command=root.quit).grid(row=4, column=3, pady=10, ipadx=15)
#########################################
lb5 = LabelFrame(root, text='Console Log')
lb5.grid(row=5, column=0, columnspan=4, ipadx=1)

text1 = ScrolledText(lb5, width=100, height=12, background='#ffffff')
text1.grid(row=0, column=0, padx=2)
text1.update()
text1.bind('<KeyPress>', lambda e: 'break' )

root.mainloop()

