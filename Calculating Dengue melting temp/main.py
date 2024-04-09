from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import matplotlib.pyplot as plt
import pandas as pd

arc_fasta_pc = 'C:\\Users\\kauan\\OneDrive\\Documentos\\Faculdade\\3 Semestre\\AlgBio\\Calculando temperatura de anelamento\\IN Dengue 21 segmentos.fasta'  # endereço do arq fasta
output_file = "C:\\Users\\kauan\\OneDrive\\Documentos\\Faculdade\\3 Semestre\\AlgBio\\Calculando temperatura de anelamento\\sequencia_fasta.txt"  # endereço do txt final
fasta_txt = "C:\\Users\\kauan\\OneDrive\\Documentos\\Faculdade\\3 Semestre\\AlgBio\\Calculando temperatura de anelamento\\fasta_txt.txt"

def Construct():
    
    datas = {'Sequence': [], 'GC percent' : [], 'Melting Temperature': []}
    
    seq_id = 1
    
    with open(output_file, "w") as f:  # abre arquivo txt destino e atribui o nome "f"
        for sequence in SeqIO.parse(arc_fasta_pc, "fasta"):
            
            Writing_bases(f, seq_id, sequence)
            
            GC = Amount_gc(f, sequence)
            
            melting_temp = Melting_temperature(f, GC, sequence)
            
            Table_append(datas, GC*100, melting_temp, seq_id)
            
            seq_id +=1
            
        Table_generator(datas)

def Writing_bases(f, seq_id, sequence):
    
    base_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    
    for base in sequence.seq:
        if base in base_counts:
            base_counts[base] += 1
            
    total = base_counts.values()
            
    # Escreve as contagens no arquivo de saída
    f.write(f"========================Bases da sequência {seq_id} ========================\n")
    f.write(f"A: {base_counts['A']}\n")
    f.write(f"T: {base_counts['T']}\n")
    f.write(f"C: {base_counts['C']}\n")
    f.write(f"G: {base_counts['G']}\n")
    f.write(f"TOTAL BASES: {sum(total)}\n\n")

def Amount_gc(f, sequence):
    
    GC = gc_fraction(sequence, ambiguous="ignore")
    
    f.write(f"%GC:\n")
    f.write(f"%{GC*100}\n")
    
    return GC

def Melting_temperature(f, GC, sequence):
    
    t = 64.9 + (0.41 * GC) - (500/len(sequence))
    
    return t

def To_txt():
     with open(fasta_txt, "w") as f:  # abre arquivo txt destino e atribui o nome "f"
        for sequence in SeqIO.parse(arc_fasta_pc, "fasta"):
            
            f.write("\n{}\n".format(sequence))

def Table_append(datas, GC, melting_temp, seq_id):
    
    datas['Sequence'].append(seq_id)
    datas['GC percent'].append(GC)
    datas['Melting Temperature'].append(melting_temp)
    
def Table_generator(datas):
    
    df = pd.DataFrame(datas)
    
    plt.plot(df['Sequence'], df['GC percent'], marker='o')
    plt.xlabel('GC percent')
    plt.ylabel('Melting Temeprature')
    plt.title('GC% X Melting Temperature (C°)')
    plt.grid(True)  # Adiciona a grade ao gráfico
    plt.xticks(range(1, 22))  # Especifica os valores no eixo x de 1 a 10
    plt.tight_layout()  # Ajustar o layout para melhor visualização
    plt.savefig('grafico.png')
    df.to_excel('tabela_GC.xlsx', index=False)

Construct()
To_txt()