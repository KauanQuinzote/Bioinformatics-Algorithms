from Bio import SeqIO

input_file = 'C:\\Users\\kauan\\OneDrive\\Documentos\\Faculdade\\3 Semestre\\Bioinformatics-Algorithms\\To identify regions of the Dengue genome that potentially encode proteins\\AF226685.2.fna'
output_file = 'C:\\Users\\kauan\\OneDrive\\Documentos\\Faculdade\\3 Semestre\\Bioinformatics-Algorithms\\To identify regions of the Dengue genome that potentially encode proteins\\potential_proteins_sequence.fasta'


def pickup_sequences(frame_1, frame_2, frame_3):
    
    sequences_frame1 = frame_1.split('*')
    sequences_frame2 = frame_2.split('*')
    sequences_frame3 = frame_3.split('*')
    
    with open(output_file, "a") as f: # abre arquivo txt destino e atribui o nome "f"
        x = y = z = 1
                  
        for sequence in sequences_frame1:
            f.write(f"AF226685.2, Frame 1, proteína {x}\n\n\t{sequence}\n\n")
            x+=1
            
        for sequence in sequences_frame2:
            f.write(f"AF226685.2, Frame 2, proteína {y}\n\n\t{sequence}\n")
            y+=1
            
        for sequence in sequences_frame3:
            f.write(f"AF226685.2, Frame 2, proteína {z}\n\n\t{sequence}\n")
            z+=1
            
        f.close()

def pickup_proteins():
    
    for sequence in SeqIO.parse(input_file, "fasta"):
        
        print(sequence)
        DNA = sequence.seq #selecione diretamente a sequencia de DNA e esqueça o resto
        RNA = DNA.transcribe() #transcreva o DNA
        
        frame_1 = RNA.translate()
        frame_2 = RNA[1:].translate()
        frame_3 = RNA[2:].translate()
    
    pickup_sequences(frame_1, frame_2, frame_3)

pickup_proteins()