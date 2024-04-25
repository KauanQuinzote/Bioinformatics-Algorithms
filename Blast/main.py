from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

input_file = "C:\\Users\\kauan\\OneDrive\\Documentos\\Faculdade\\3 Semestre\\Bioinformatics-Algorithms\\Blast\\sequence.gp"
output_file = "C:\\Users\\kauan\\OneDrive\\Documentos\\Faculdade\\3 Semestre\\Bioinformatics-Algorithms\\Blast\\regions.fasta"
amino_file =  R"C:\Users\kauan\OneDrive\Documentos\Faculdade\3 Semestre\Bioinformatics-Algorithms\Blast\sequence.fasta"

x = 1

def extrat_algarithms(string): #this will take the algarisms at Region's line
    start_seq = ''
    end_seq = ''
    x = int
    
    for char in string:
        
        if char.isdigit():     
            start_seq +=char 
        if char == '.':
            x = string.index(char)
            break
        
    for char in string[x:]:
        if char.isdigit():  
            end_seq +=char
    
    return start_seq, end_seq

def making_fasta(start_seq, end_seq): #this will make the fasta with a sliced sequence from the polyprotein that you chase
    
    start_seq = int(start_seq)
    end_seq = int(end_seq)
    
    record = SeqIO.read(amino_file, "fasta")
    record = record.seq
    
    slice_record = ''
    
    for index, char in enumerate(record):
        if index >= start_seq and index <= end_seq:
            slice_record += char
    
    with open(output_file, "a") as f: #write the sliced sequence
        SeqIO.write(SeqRecord(slice_record, id=f'sliced_{x}', description='Sliced from Dengue_AAN60367.2'), f, 'fasta')

def main():
    global x
    
    with open (input_file, 'r') as f: #open a fasta file that contains your polyprotein sequence
        content = f.readlines() #this will read
    
    for line in content: #for here, i wanna slice just the regions
        if 'Region' in line:       
            start_seq, end_seq = extrat_algarithms(line)
            making_fasta(start_seq, end_seq)
            x+=1
main()