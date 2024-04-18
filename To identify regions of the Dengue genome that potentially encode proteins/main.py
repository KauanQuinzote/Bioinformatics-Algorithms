from Bio import SeqIO

input_file = 'C:\\Users\\kauan\\OneDrive\\Documentos\\Faculdade\\3 Semestre\\Bioinformatics-Algorithms\\To identify regions of the Dengue genome that potentially encode proteins\\AF226685.2.fna'
output_file = 'C:\\Users\\kauan\\OneDrive\\Documentos\\Faculdade\\3 Semestre\\Bioinformatics-Algorithms\\To identify regions of the Dengue genome that potentially encode proteins\\potential_proteins_sequence.fasta'

def pickup_proteins():
    
    frame = 0
    amino = []
    
    while frame < 3:
    
        for sequence in SeqIO.parse(input_file, "fasta"):
                
            gen = sequence.seq.translate() #selecione diretamente a sequencia de DNA e traduza
                     
            for x in gen:
                if x == 'M':
                    
                    while x != '*':
                        amino.append(x)

            print(f'{amino}\n')
                
        frame+=1

pickup_proteins()