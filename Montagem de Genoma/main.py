from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def encontrar_subsequencia_comum(sequencia1, sequencia2):
    #essa função encontra a maior subsequencia comum entre duas sequencias recebidas
    len1, len2 = len(sequencia1), len(sequencia2)
    max_len = 0
    fim = 0
    for i in range(len1):
        for j in range(len2):
            l = 0
            while (i + l < len1 and j + l < len2 and sequencia1[i + l] == sequencia2[j + l]):
                l += 1
            if l > max_len:
                max_len = l
                fim = i + l
    return sequencia1[fim - max_len:fim]

def concatenar_sem_comum(sequencia1, sequencia2):
    #essa função concatena as 2 sequencias após achar a maior subsequencia comum entre elas
    subsequencia_comum = encontrar_subsequencia_comum(sequencia1, sequencia2)
    if not subsequencia_comum:
        return 'Não compatível'
    # Remover a subsequência comum apenas uma vez
    pos1 = sequencia1.find(subsequencia_comum)
    pos2 = sequencia2.find(subsequencia_comum)
    parte1 = sequencia1[:pos1]
    parte2 = sequencia2[pos2 + len(subsequencia_comum):]
    return parte1 + subsequencia_comum + parte2

def OpenFasta():
    #abrindo o arquivo
    lst = []
    
    #pegando cada sequencia no reads4.fasta    
    for seq in SeqIO.parse("reads4.fasta", "fasta"):
        read = seq.seq        
        lst.append(read)
        
    #combinando cada sequencia dessa com sequencia sucessora
    for i in range(3):    
        word = concatenar_sem_comum(lst[i], lst[i+1])
        print(word)
        lst[i+1] = word
        
    #transformando objetos Seq em Seqrecord
    seq_record = SeqRecord(lst[3])
    seq_record.id = "No ID"
    
    #gravando no fasta
    SeqIO.write(seq_record, "contig.fasta", "fasta")
    
OpenFasta()    