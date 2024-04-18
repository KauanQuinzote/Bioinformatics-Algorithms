from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#arch paths
fasta_file = "C:\\Users\\kauan\\OneDrive\\Documentos\\Faculdade\\3 Semestre\\Bioinformatics-Algorithms\\Potential Encoded Dengues RNA Proteins\\IN Dengue 21 segmentos.fasta"
output_fasta = 'C:\\Users\\kauan\\OneDrive\\Documentos\\Faculdade\\3 Semestre\\Bioinformatics-Algorithms\\Potential Encoded Dengues RNA Proteins\\Potential Encoded Proteins.fasta'


def reading_dengue_fasta():

    start_codon = "M" #the start codon parameter 
    end_codon = "*"  #end codon parameter
    records = []
        
    for seq_record in SeqIO.parse(fasta_file, "fasta"): #for each sequence at fasta file
       
        seq_id = seq_record.id
        sequence = str(seq_record.seq)    
        
        for frame in range(3):#generate each frame
            
            aft_star = None 
            protein_num = 1
            i = 0 #the stop index to return to read the sequence
            
            dna_seq = Seq(sequence[frame:])
            protein_seq = dna_seq.translate()
        
            while i < len(protein_seq): #while the sequence is not finished
                
                if protein_seq[i] == start_codon: #ensures if char is a start codon
                
                    if aft_star is None: #then, change the index to continue reading the seq
                        aft_star = i
                    
                elif protein_seq[i] == end_codon: #ensures if char is a end condon
                    
                    if aft_star is not None: #if the index is to continue reading is true
                
                        amino = protein_seq[aft_star:i + 1] #split the amino from sequence
                        record = SeqRecord(amino, id=f"{seq_id}, Frame {frame + 1}, proteína {protein_num}\n\n", description="") #make a SeqRecord obj from the splited sequence
                        
                        records.append(record) #append the list
                        protein_num += 1

                        aft_star = None
                i += 1 #append caractere from de sequence
                
        SeqIO.write(records, output_fasta, "fasta")#write all at doc


reading_dengue_fasta() #call the script