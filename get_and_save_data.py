import GEOparse
import pandas as pd


def extract_data(gse, dest):
    """
    Extracts the expression matrix and the phenotypic variables of the selected set using the geoparse library and saves them in csv format. 

    Parameters:
    gse: Genetic set to be extracted
    dest: Destination folder where the data will be saved
    """
    gse = GEOparse.get_GEO(geo=gse, destdir= dest)
    # Get expression data and metadata matrices
    exprs = []
    gsmNames = []
    metadata = {}
    for gsm_name, gsm in gse.gsms.items():
        if gsm.metadata['type'][0]=='RNA':
            # Expression data
            if len(gsm.table)>0:
                tmp = gsm.table['VALUE']
                tmp.index = gsm.table['ID_REF']
                gsmNames.append(gsm_name)
                if len(exprs)==0:
                    exprs = tmp.to_frame()
                else:
                    exprs = pd.concat([exprs,tmp.to_frame()],axis=1)
            
            # Metadata
            for key,value in gsm.metadata.items():
                #print(key)
                #print(value)
                if (key=='characteristics_ch1' or key=='characteristics_ch2') and (len([i for i in value if i!=''])>1 or value[0].find(': ')!=-1):
                    #print(value)
                    tmpVal = 0
                    for tmp in value:
                        splitUp = [i.strip() for i in tmp.split(':')]
                        #print(splitUp)
                        if len(splitUp)==2:
                            if not splitUp[0] in metadata:
                                metadata[splitUp[0]] = {}
                            metadata[splitUp[0]][gsm_name] = splitUp[1]
                        else:
                            if not key in metadata:
                                metadata[key] = {}
                            metadata[key][gsm_name] = splitUp[0]
                else:
                    if not key in metadata:
                        metadata[key] = {}
                    if len(value)==1:
                        metadata[key][gsm_name] = ' '.join([j.replace(',',' ') for j in value])

    # Write expression data matrix to file
    exprs.columns = gsmNames
    with open(dest+'exprs_data.csv','w+') as outFile:
        exprs.to_csv(outFile)
    
    # Write metadata matrix to file
    with open(dest+'metadata.csv','w') as outFile:
        outFile.write('ID_REF,'+','.join(gsmNames))
        for key in metadata:
            tmp = [key]
            for gsm_name in gsmNames:
                if gsm_name in metadata[key]:
                    tmp.append(metadata[key][gsm_name])
                else:
                    tmp.append('NA')
            outFile.write('\n'+','.join(tmp))


def get_data(dest):
    """
    Accesses the selected folder to read the data stored in it.

    Parameters:
    dest: Folder where we extract your data

    Returns:
    expres: Expression matrix
    metadata: Phenotypic variables
    """
    expres = pd.read_csv(dest+'exprs_data.csv')

    metadata_ = pd.read_csv(dest+'metadata.csv')

    metadata = metadata_.set_index('ID_REF').T

    return expres, metadata



