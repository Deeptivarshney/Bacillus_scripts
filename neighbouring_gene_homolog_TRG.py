import os
from Bio import SeqIO
import json
import argparse
import copy

#python3 neighbouring_gene_homolog_TRG.py  --indir1 Bacillus_Database_GFF_Sorted_Contig_info_add --indir2 Step3_homolog_coordinates_info_output --taxid 1386 --txt /media/data/TRG_project/NCBI_bacteria/Bacillus/1386.Species_TRGs_blastagainst_bacillus-genus.idea2result.txt --out Neighbouring_gene_homolog_TRG

parser = argparse.ArgumentParser(description='Process some integers.')

parser.add_argument('--txt', '-txt', '-t', dest='txt', metavar="FILE",
                    help='TXT file for TRGs ids')

parser.add_argument('--indir1', '--i1', '-i1', metavar="DIRECTORY", type=str, default="Bacillus_Database_GFF_Sorted_Contig_info_add",
                    help='input directory with downloaded genus subdirectories [DEFAULT: %(default)s]')


parser.add_argument('--indir2', '--i2', '-i2', metavar="DIRECTORY", type=str, default="Step3_homolog_coordinates_info_output",
                    help='input directory with genus subdirectories for homologs for each species [DEFAULT: %(default)s]')

parser.add_argument('--taxid', '-taxid', dest='taxid', type=int,  default=1386,
                    help='Taxonomy ID of query group of organisms [DEFAULT: %(default)s]')
parser.add_argument('--outdir', '--o', '-o', metavar="DIRECTORY", type=str, default="Neighbouring_gene_homolog_TRG",
                    help='output directory to save results [DEFAULT: %(default)s]')


args = parser.parse_args()


def function_for_left(f_species, f_otherstrain,f_species_to_strain, index_dic,left,f_homolog, f_dic_assembly):
    
    flag = None
    if left in homolog.keys():

        for otherSpeciesGeneID in f_homolog[left]:
            for otheraccession in index_dic[f_otherstrain]:
                if otherSpeciesGeneID in index_dic[f_otherstrain][otheraccession]:
                    #flag = left+'\t'+otherSpeciesGeneID+'\t'+'\t'.join(f_dic_assembly[f_otherstrain][otheraccession][otherSpeciesGeneID])+'\t'+ otheraccession+'\t'+f_otherstrain 
                    flag = otherSpeciesGeneID+'\t'+'\t'.join(f_dic_assembly[f_otherstrain][otheraccession][otherSpeciesGeneID])+'\t'+ otheraccession+'\t'+f_otherstrain 

                    if flag != None:
                    #flag = left,f_dic_assembly[f_otherstrain][otheraccession][left],otherSpeciesGeneID,f_dic_assembly[f_otherstrain][otheraccession][otherSpeciesGeneID], otheraccession,f_otherstrain
                        break
            if flag != None:
                break
    else:
        flag = None
    #print(flag)
    return(flag)



IN_GENUS_PATH = os.path.join(args.indir1, str(args.taxid))
if not os.path.exists(IN_GENUS_PATH):
    print('Input directory {} does not exist.'.format(IN_GENUS_PATH))
    exit()

OUT_GENUS_PATH = os.path.join(args.outdir, str(args.taxid))
if not os.path.exists(OUT_GENUS_PATH):
    os.makedirs(OUT_GENUS_PATH)



dic_assembly = {}
Index_Taxid = {}
Species_to_strain = {}

for species_taxid in os.listdir(IN_GENUS_PATH):
    in_species = os.path.join(IN_GENUS_PATH, species_taxid)


    
    if species_taxid not in Species_to_strain:
        Species_to_strain[species_taxid] = []

    for gffassembly_dir in os.listdir(in_species):
        gffassembly_path = os.path.join(in_species,gffassembly_dir)
        Species_to_strain[species_taxid].append(gffassembly_dir)

        gffassembly_name = '{}_genomic_gff.txt'.format(gffassembly_dir)
        gff_assembly_name_path = os.path.join(gffassembly_path, gffassembly_name)
        
        
        gffopen = open(gff_assembly_name_path)
        if gffassembly_dir not in dic_assembly and gffassembly_dir not in Index_Taxid:
            dic_assembly[gffassembly_dir] = {}
            
            Index_Taxid[gffassembly_dir] = {}
            
            for evline in gffopen:
                splitevline = evline.split()

                if splitevline[1] not in dic_assembly[gffassembly_dir]:
                    dic_assembly[gffassembly_dir][splitevline[1]] = {}

                if splitevline[0] not in dic_assembly[gffassembly_dir][splitevline[1]]:
                    dic_assembly[gffassembly_dir][splitevline[1]][splitevline[0]] = (splitevline[3],splitevline[4])
                    
                if splitevline[1] not in Index_Taxid[gffassembly_dir]:
                    Index_Taxid[gffassembly_dir][splitevline[1]]= []
                Index_Taxid[gffassembly_dir][splitevline[1]].append(splitevline[0])

        
        gffopen.close()

# for key, value in Species_to_strain.items():
#     print(key,value)

IN_HOMOLOG_PATH = os.path.join(args.indir2, str(args.taxid))

#oh = open('homolog_left.3.txt','a')
trgfile = open(args.txt)
fh = trgfile.readlines()[1:]
for line in fh:
    print('Query species line '+line)
    sl = line.rstrip( ).split("\t")
    species = sl[0]
    out_species = os.path.join(OUT_GENUS_PATH, species)
    if not os.path.exists(out_species):
        os.makedirs(out_species)
    

    ## Make dictionary for homologs ids for working query species 

    in_species_path = os.path.join(IN_HOMOLOG_PATH, species)

    homolog = {}
    
    for query_homo_info in os.listdir(in_species_path):
        query_filename = '{}'.format(query_homo_info)
        query_homofilenamepath = os.path.join(in_species_path,query_filename)
        #print(query_homofilenamepath)
        
        homologresultfile = open(query_homofilenamepath)
        for homologline in homologresultfile:
            hmsplit = homologline.split()
            queryid = hmsplit[0]
            subid = hmsplit[1]
                        
            if queryid not in homolog:
                homolog[queryid] = set()
            homolog[queryid].add(subid)      ## Done -  dictionary for homologs ids for working query species 



# ## continue: work with TRG_txt file      
    
    TRG_ids_column = sl[4]
    splitTRG_id = TRG_ids_column.split('|')
    
    for singleTRG in splitTRG_id:
        out_filename = '{}.TRG.txt'.format(singleTRG)
        out_filenamepath = os.path.join(out_species, out_filename)
        oh = open(out_filenamepath,'w')
        for Strains in Species_to_strain[species]:
            for q_accession in Index_Taxid[Strains]:
                if singleTRG in Index_Taxid[Strains][q_accession]:
                    indexofsingleTRG = Index_Taxid[Strains][q_accession].index(singleTRG)
                    #left = Index_Taxid[Strains][q_accession][indexofsingleTRG-1]
                    ##right = Index_Taxid[Strains][q_accession][indexofsingleTRG+1]
                    #print(left,indexofsingleTRG, singleTRG, right,indexofsingleTRG+1, Strains)
                    #return_hit_for_left = None 
                    #for index in reversed(range(0, indexofsingleTRG)):
                    if indexofsingleTRG > 0:
                        index = indexofsingleTRG-1
                        left = Index_Taxid[Strains][q_accession][index]
                        indexSwap=index
                        leftAccnPos=left,q_accession,Strains,dic_assembly[Strains][q_accession][left]
                        #print (left, dic_assembly[Strains][q_accession][left],Strains, q_accession)
                        leftSwap=left
                        for outSpecies in Species_to_strain.keys():
                            if outSpecies!=species:
                                for outerstrains in Species_to_strain[outSpecies]:
                                    #print ('Left',singleTRG)
                                    return_hit_for_left = None
                                    while return_hit_for_left==None and index >0:
                                        return_hit_for_left = function_for_left(outSpecies,outerstrains,Species_to_strain, Index_Taxid,left,homolog,dic_assembly)
                                        #leftAccnPos='left'+'\t'+str(index)+'\t'+left+'\t'+q_accession+'\t'+Strains+'\t'+'\t'.join(dic_assembly[Strains][q_accession][left])
                                        leftAccnPos='left'+'\t'+str(index)+'\t'+left+'\t'+'\t'.join(dic_assembly[Strains][q_accession][left])+'\t'+q_accession+'\t'+Strains
                                        index= index-1
                                        left = Index_Taxid[Strains][q_accession][index]
                                            
                                    oh.write(str(leftAccnPos)+'\t'+str(return_hit_for_left)+'\n')
                                    #print(leftAccnPos,return_hit_for_left)
                                    #dic_assembly[Strains][q_accession][left]
                                    left=leftSwap
                                    index=indexSwap
                    if indexofsingleTRG < len(Index_Taxid[Strains][q_accession]) and len(Index_Taxid[Strains][q_accession])!=1:
                        index_right = indexofsingleTRG+1
                        right = Index_Taxid[Strains][q_accession][index_right]
                        right_indexSwap=index_right
                        rightAccnPos=right,q_accession,Strains,dic_assembly[Strains][q_accession][right]
                        rightSwap=right
                        
                        for right_outSpecies in Species_to_strain.keys():
                            if right_outSpecies!=species:
                                for right_outerstrains in Species_to_strain[right_outSpecies]:
                                    return_hit_for_right = None
                                    while return_hit_for_right==None and index_right <=len(Index_Taxid[Strains][q_accession]):
                                        return_hit_for_right = function_for_left(right_outSpecies,right_outerstrains,Species_to_strain, Index_Taxid,right,homolog,dic_assembly)
                                        #rightAccnPos='right' +'\t'+str(index_right)+'\t'+right+'\t'+q_accession+'\t'+Strains+'\t'+'\t'.join(dic_assembly[Strains][q_accession][right])
                                        rightAccnPos='right'+'\t'+str(index_right)+'\t'+right+'\t'+'\t'.join(dic_assembly[Strains][q_accession][right])+'\t'+q_accession+'\t'+Strains
                                        index_right = index_right+1
                                        if index_right>len(Index_Taxid[Strains][q_accession]):
                                            break
                                        try:
                                            right = Index_Taxid[Strains][q_accession][index_right]
                                        except:
                                            print (rightAccnPos,right,index_right)
                                            break


                                    
                                    oh.write(str(rightAccnPos)+'\t'+str(return_hit_for_right)+'\n')
                                    #print(rightAccnPos,return_hit_for_right)
                                    right=rightSwap
                                    index_right=right_indexSwap
                                    


                        #right gene
    oh.close()
    

        
