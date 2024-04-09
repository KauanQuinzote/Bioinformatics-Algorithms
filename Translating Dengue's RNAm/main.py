from Bio import SeqIO
from Bio.Seq import Seq

input_file = 'C:\\Users\\kauan\\OneDrive\\Documentos\\Faculdade\\3 Semestre\\Bioinformatics-Algorithms\\Tradução de RNA\\IN Dengue 21 segmentos corrigido.fasta'
output_RNA = 'C:\\Users\\kauan\\OneDrive\\Documentos\\Faculdade\\3 Semestre\\Bioinformatics-Algorithms\\Tradução de RNA\\Dengue_RNA_Transcrito.fasta'
protein_dir = 'C:\\Users\\kauan\\OneDrive\\Documentos\\Faculdade\\3 Semestre\\Bioinformatics-Algorithms\\Tradução de RNA\\'

def count_string_backflip(sequence, substring, frame_id):
    
    reverse_sequence = sequence[::-1]
    # Dependendo de como você está contando os codons em cada frame, você pode precisar ajustar esta função
    start_index = frame_id - 1  # Ajuste o índice de início com base no frame
    count = 0
    index = reverse_sequence.find(substring, start_index)  # Encontre a primeira ocorrência da substring no frame
    while index != -1:
        count += 1
        index = reverse_sequence.find(substring, index + 1)  # Encontre a próxima ocorrência da substring
    return count

def backfliptranscript(RNA, seq_id):
    
    output_frames = 'C:\\Users\\kauan\\OneDrive\\Documentos\\Faculdade\\3 Semestre\\Bioinformatics-Algorithms\\Tradução de RNA\\'
    frame_id = 4
    
    aminos = { #dicionário que o Shida pediu com os aminoácidos
        "AUG": "Metionina",
        "UUU": "Fenilalanina",
        "UUC": "Fenilalanina",
        "UUA": "Leucina",
        "UUG": "Leucina",
        "CUU": "Leucina",
        "CUC": "Leucina",
        "CUA": "Leucina",
        "CUG": "Leucina",
        "AUU": "Isoleucina",
        "AUC": "Isoleucina",
        "AUA": "Isoleucina",
        "GUU": "Valina",
        "GUC": "Valina",
        "GUA": "Valina",
        "GUG": "Valina",
        "UGU": "Cisteína",
        "UGC": "Cisteína",
        "UGG": "Triptofano",
        "UAA": "Codão de terminação",
        "UAG": "Codão de terminação",
        "UGA": "Codão de terminação"
    }
    
    while frame_id < 7:
        
        amount_aminos = {key: 0 for key in aminos.keys()}# Inicialize o contador para cada codon

        with open(f'{output_frames}Frame_{frame_id}.txt', 'a') as f:
            f.write(f"Sequencia : {seq_id}\n")
            
            for key in aminos.keys():
                
                amount_aminos[key] = count_string_backflip(RNA, key, frame_id)
                
                f.write(f'\t{aminos[key]} : {amount_aminos[key]}\n')
            f.write("\n")    
                
        frame_id +=1       

def count_string(sequence, substring, frame_id):
    # Dependendo de como você está contando os codons em cada frame, você pode precisar ajustar esta função
    start_index = frame_id - 1  # Ajuste o índice de início com base no frame
    count = 0
    index = sequence.find(substring, start_index)  # Encontre a primeira ocorrência da substring no frame
    while index != -1:
        count += 1
        index = sequence.find(substring, index + 1)  # Encontre a próxima ocorrência da substring
    return count

def translate(RNA, seq_id):
    frame_id = 1;
    output_frames = 'C:\\Users\\kauan\\OneDrive\\Documentos\\Faculdade\\3 Semestre\\Bioinformatics-Algorithms\\Tradução de RNA\\'
    
    while frame_id < 4:
        
        translated = RNA[frame_id-1:].translate()
        
        with open(f'{output_frames}Translated_{frame_id}.txt', 'a') as f:
            f.write(f"Sequencia : {seq_id}\n")
            f.write(f'\t{translated}\n')
            
        frame_id +=1
    
def translate_proteins(RNA, seq_id):
        
    output_frames = 'C:\\Users\\kauan\\OneDrive\\Documentos\\Faculdade\\3 Semestre\\Bioinformatics-Algorithms\\Tradução de RNA\\'
    frame_id = 1
    
    aminos = { #dicionário que o Shida pediu com os aminoácidos
        "AUG": "Metionina",
        "UUU": "Fenilalanina",
        "UUC": "Fenilalanina",
        "UUA": "Leucina",
        "UUG": "Leucina",
        "CUU": "Leucina",
        "CUC": "Leucina",
        "CUA": "Leucina",
        "CUG": "Leucina",
        "AUU": "Isoleucina",
        "AUC": "Isoleucina",
        "AUA": "Isoleucina",
        "GUU": "Valina",
        "GUC": "Valina",
        "GUA": "Valina",
        "GUG": "Valina",
        "UGU": "Cisteína",
        "UGC": "Cisteína",
        "UGG": "Triptofano",
        "UAA": "Codão de terminação",
        "UAG": "Codão de terminação",
        "UGA": "Codão de terminação"
    }
    
    while frame_id < 4:
        
        amount_aminos = {key: 0 for key in aminos.keys()}# Inicialize o contador para cada codon

        with open(f'{output_frames}Frame_{frame_id}.txt', 'a') as f:
            f.write(f"Sequencia : {seq_id}\n")
            
            for key in aminos.keys():
                
                amount_aminos[key] = count_string(RNA, key, frame_id)
                
                f.write(f'\t{aminos[key]} : {amount_aminos[key]}\n')
            f.write("\n")    
                
        frame_id +=1
        
    backfliptranscript(RNA, seq_id)
    translate(RNA, seq_id)  

def transcript():
    
    seq_id = 1
       
    with open(output_RNA, "w") as f: # abre arquivo txt destino e atribui o nome "f"
        for sequence in SeqIO.parse(input_file, "fasta"): #para cada ID no arquivo fasta
            DNA = sequence.seq #selecione diretamente a sequencia de DNA e esqueça o resto
            RNA = DNA.transcribe() #transcreva o DNA

            f.write(f'Sequência de RNA: {seq_id}\n{RNA}\n\n') #escreva no arquivo RNA Degue.fasta isso
            translate_proteins(RNA, seq_id)
            seq_id +=1                

transcript()
