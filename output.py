A"""
This module is used to output and visualize the 
results
"""

import analysis as ana
import go_parser as go
import pandas as pd


def findGroups(GO , binary, mod_names, min_sup=10):
    """
    This function returns the patterns and its corresponding 
    result in a panda format
    """

    report = ana.findFSGenes(ana.__findFS(binary,min_sup),\
                             mod_names,binary)

    
    
    dfs = []
    patterns = []
    lengths = []
    
    for key in report:
        mods = [mod_names[x] for x in key]
        pattern = ','.join(mods)
        modL = len(mods)
        
        glst = report[key]
        df = GO[GO.index.isin(glst)]
        length = len(df)

        patterns+= [pattern]*length
        lengths += [modL]*length

        dfs.append(df)

    result = pd.concat(dfs)
    result['pattern']=pd.Series(patterns,index=result.index)

    result['modL']=pd.Series(lengths,index=result.index)
    return result


