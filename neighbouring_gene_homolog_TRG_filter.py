
## python neighbouring_gene_homolog_TRG_filter.py
import os



def for_difference(function_key):
    Output = None
    
    left_nz = function_key[0][-2]
    left_end = function_key[0][-3]
    
    right_nz = function_key[1][-2]
    right_start = function_key[1][-4]
    
    if right_nz == left_nz:
        calculate_diff = int(right_start) - int(left_end) 
        

        Output = calculate_diff


    return(Output)




GENUS_TAXID = '1386'

INDIR = 'Neighbouring_gene_homolog_TRG'
OUTDIR = 'Neighbouring_gene_homolog_TRG_Filtered'


IN_GENUS_PATH = os.path.join(INDIR, GENUS_TAXID)
if not os.path.exists(IN_GENUS_PATH):
    print('Input directory {} does not exist.'.format(IN_GENUS_PATH))
    exit()

OUT_GENUS_PATH = os.path.join(OUTDIR, GENUS_TAXID)
if not os.path.exists(OUT_GENUS_PATH):
    os.makedirs(OUT_GENUS_PATH)


for species_taxid in os.listdir(IN_GENUS_PATH):
    in_species_path = os.path.join(IN_GENUS_PATH, species_taxid)
    
    out_species_path = os.path.join(OUT_GENUS_PATH, species_taxid)
    if not os.path.exists(out_species_path):
        os.makedirs(out_species_path)
    
    for files in os.listdir(in_species_path):
        in_files_path = os.path.join(in_species_path,files)
        print(in_files_path)

    
        out_files_path = os.path.join(out_species_path, files)
        
        out_file = open(out_files_path, 'w')
        out_file.write('Q.strain\tQ.accession\tleftgene\tstart\tend\trightgene\tstart\tend\tSub.strain\tSub.accession\tleft_homolog_id\tstart\tend\tright_homolog_id\tstart\tend\tdifference\n')

        Pairing_dic = {}
        
        fh = open(in_files_path)
        for line in fh:
            if line.startswith('('):
                continue
            sl = line.rstrip( ).split("\t")
            

            if sl[7] == 'None':
                GCF_key = line
            else:
                GCF_key = sl[5] +'\t'+ sl[6]+'\t'+ sl[10] + '\t'+ sl[11] 
                
                
            if GCF_key not in Pairing_dic:
                Pairing_dic[GCF_key] = []
            Pairing_dic[GCF_key].append(sl)

        fh.close()  

                                 
        for k, v  in Pairing_dic.items():
            if not k.startswith('left') and not k.startswith('right') and len(Pairing_dic[k])!=1:
                
                diff_return = for_difference(Pairing_dic[k])
                
                
                out_file.write(str(v[0][6])+ '\t'+str(v[0][5]) +'\t'+ str(v[0][2])+'\t'+str('\t'.join(v[0][3:5]))+'\t'+str(v[1][2])+'\t'+str('\t'.join(v[1][3:5]))+'\t'+ str(v[0][11])+'\t'+str(v[0][10])+'\t'+ str(v[0][7])+'\t'+str('\t'.join(v[0][8:10]))+'\t'+str(v[1][7])+'\t'+str('\t'.join(v[1][8:10]))+'\t'+str(diff_return)+'\n')
                
                
        out_file.close()
                   
            
