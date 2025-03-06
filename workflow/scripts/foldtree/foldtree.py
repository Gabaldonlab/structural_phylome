
import foldseek2tree
import corecut 
import argparse
import numpy as np
import pandas as pd
import os

def structblob2tree(input_folder, outfolder, tmpfolder, corefolder, outtree, 
                    overwrite = False, correction = False,
                    fastmepath = 'fastme', foldseekpath = 'foldseek', delta = 0.0001,
                    kernel = 'fident', core = False, hitthresh = .8, minthresh = .6):
    '''
    run fold tree pipeline for a folder of pdb files
    
    Parameters
    ----------
    input_folder : str
        path to folder with pdb files   
    outfolder : str
        path to output folder   
    overwrite : bool
        overwrite existing foldseek output  
    fastmepath : str    
        path to fastme executable
    quicktreepath : str 
        path to quicktree executable
    foldseekpath : str  
        path to foldseek executable 
    delta : float   
        small number to replace negative branch lengths with, default is .0001
    correction : str    
        correction method to use, either 'tajima' or 'none'
    kernel : str    
        kernel to use, either 'fident', 'lddt' or 'alntmscore'
    
    '''
    #check if the foldseek output is already there
    input_folder = os.path.join(input_folder, '')
    if os.path.exists(outfolder + '_allvall.tsv') and overwrite == False:
        print('found foldseek output, skipping foldseek')
        alnres = outfolder + '_allvall.tsv'
    else:
        alnres = foldseek2tree.runFoldseek_allvall_EZsearch(input_folder, outfolder + '_allvall.tsv', tmpfolder, foldseekpath = foldseekpath)
    
    if core == True:
        corefolder = os.path.join(corefolder, '')
        corecut.extract_core(alnres, outfolder+'.core',  hitthresh = hitthresh, minthresh = minthresh, 
                             corefolder = corefolder, structfolder = input_folder)
        if os.path.exists(outfolder + '_core_allvall.tsv') and overwrite == False:
            print('found foldseek core output, skipping foldseek')
            alnres = outfolder + '_core_allvall.tsv'
        else:
            alnres = foldseek2tree.runFoldseek_allvall_EZsearch(corefolder, outfolder + '_core_allvall.tsv', tmpfolder, foldseekpath = foldseekpath)

    res = pd.read_table(alnres, header = None)
    res[0] = res[0].map(lambda x :x.replace('.core.pdb', ''))
    res[1] = res[1].map(lambda x :x.replace('.core.pdb', ''))
    res.columns = 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,lddtfull,alntmscore'.split(',')
    
    ids = list(set(list(res['query'].unique()) + list(res['target'].unique())))
    pos = {protid : i for i,protid in enumerate(ids)}
    matrices = {kernel:np.zeros((len(pos), len(pos)))}

    #calc kernel for tm, aln score, lddt
    for idx,row in res.iterrows():
        for k in matrices:
            matrices[k][pos[row['query']] , pos[row['target']]] += row[k]
            matrices[k][pos[row['target']] , pos[row['query']]] += row[k]
    trees = {}
    for i,k in enumerate(matrices):
        matrices[k] /= 2
        matrices[k] = 1-matrices[k]
        matrices[k] = np.clip(matrices[k], 0, 1)
        
        # print(matrices[k], np.amax(matrices[k]), np.amin(matrices[k]) )
        if correction:
            if kernel == 'fident':
                factor = 19/20
            else:
                factor = 1
            matrices[k] = foldseek2tree.Tajima_dist(matrices[k], bfactor = factor)
        # np.save(outfolder + k + '_distmat.npy' , matrices[k])
        distmat_txt = foldseek2tree.distmat_to_txt(ids, matrices[k], outtree.replace(".nwk", ".txt"))
        out_tree = foldseek2tree.runFastme(fastmepath = fastmepath, clusterfile = distmat_txt)
        # out_tree = foldseek2tree.postprocess(out_tree, outfolder + '_' + k + '.nwk', delta = delta)
        out_tree = foldseek2tree.postprocess(out_tree, outtree, delta = delta)
        trees[k] = out_tree
    return alnres, trees

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='run foldtree pipeline for a folder of pdb files')
    # Files
    parser.add_argument("-i", "--in", dest="struct_dir", default="None", help="input structure folder", required=True)
    parser.add_argument("-o", "--out", dest="output_dir", default="None", help="output prefix", required=True)
    parser.add_argument("--outtree", dest="outtree", default="None", help="output tree name", required=True)
    parser.add_argument("-c", "--core", dest="core_dir", help="Tmp folder for foldseek", required=True)
    parser.add_argument("-t", "--tmp", dest="tmp_dir", default="/tmp/foldtree", help="Tmp folder for foldseek")
    parser.add_argument('--overwrite', help='overwrite existing foldseek output', action='store_true')

    parser.add_argument('--kernel', choices = ['lddt', 'fident', 'tmalign' ], default = 'fident',
    help='this is the comparison metric used to build the tree')

    # corecut module
    parser.add_argument('--corecut', help='cut the core of the proteins and realign', action='store_true')
    parser.add_argument('--hitthresh', default = .8, help='threshold for finding the boundaries of the core')
    parser.add_argument('--minthresh', default = .6, help='threshold if the core is not found')

    # Distance correction
    parser.add_argument('--correction', help='use the -ln correction for the distance matrix', action='store_true')
    parser.add_argument('--delta', default=0.0001, help='small number to replace negative branch lengths with')

    # binaries
    # parser.add_argument('--fastmepath', default='fastme', help='path to fastme binary')
    # parser.add_argument('--quicktreepath', default='quicktree', help='path to quicktree binary')
    # parser.add_argument('--foldseekpath', default='foldseek', help='path to foldseek binary')
    
    args = parser.parse_args()
    # if not all([args.positional1, args.positional2]):
    #     parser.error("Positional arguments are required.")
    # flag_dict = {}
    # flag_dict['corecut'] = True if args.corecut else False
    # flag_dict['hittresh'] = args.hittresh
    # flag_dict['minthresh'] = args.minthresh
    # flag_dict['correction'] = True if args.correction else False
    # flag_dict['overwrite'] = True if args.overwrite else False
    # flag_dict['fastmepath'] = args.fastmepath
    # flag_dict['quicktreepath'] = args.quicktreepath
    # flag_dict['foldseekpath'] = args.foldseekpath
    # flag_dict['delta'] = args.delta
    # flag_dict['kernel'] = args.kernel
    print(f"struct dir: {args.struct_dir}, output dir: {args.output_dir}")
    structblob2tree(args.struct_dir, args.output_dir, args.tmp_dir, args.core_dir, args.outtree,
                    overwrite = args.overwrite, correction = args.correction,
                    delta = args.delta, kernel = args.kernel, 
                    core = args.corecut, hitthresh = args.hitthresh, minthresh = args.minthresh)
    print('Done!')

