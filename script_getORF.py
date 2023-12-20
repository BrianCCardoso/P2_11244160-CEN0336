#!/usr/bin/env python3
#importacões necessárias no programa
import re
import sys

#declarando varáveis
dic_genes = {}

class NotFASTAError(Exception):
    pass

def open_reading_frame(dic_genes):
    dic_genes_openrc = {}
    for id in dic_genes.keys():
        dic_genes_openrc[id] = {'1': '','2' : '','3' : '','4' : '','5' : '','6' : ''}
        tratando_sequencia = dic_genes[id]
        tratando_sequencia_reverso = tratando_sequencia[::-1].replace('G','c').replace('C','g').replace('T','a').replace('A','t').upper()
        for i in range(0,3):
            dic_genes_openrc[id][str(i+1)] = tratando_sequencia[i:]
            dic_genes_openrc[id][str(i+4)] = tratando_sequencia_reverso[i:]
            i+=1
    return dic_genes_openrc

def max_reading_frame(dic_openrc):
    dic_max_ORF = {}
    for id in dic_openrc.keys():
        dic_max_ORF[id] = {'Sequência' : '', 'Começo' : '', 'Fim' : '', 'Frame de leitura' : '', 'Tamanho' : 0}
        for frame in dic_openrc[id]:
            sequencia_compilada = ''.join(dic_openrc[id][frame])
            for match in re.finditer(r'ATG.*?(TAA|TGA|TAG)',sequencia_compilada):
                orf = match.group()
                if len(orf) > len(dic_max_ORF[id]['Sequência']) and len(orf) % 3 == 0:
                    dic_max_ORF[id]['Sequência'] = orf
                    dic_max_ORF[id]['Começo'] = match.start() + 1
                    dic_max_ORF[id]['Fim'] = match.end() +1
                    dic_max_ORF[id]['Frame'] = int(frame)
                    dic_max_ORF[id]['Tamanho'] = len(orf)
    return dic_max_ORF

def tranlation(dic_max):
    translation_table = {
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
    'AAT':'N', 'AAC':'N',
    'GAT':'D', 'GAC':'D',
    'TGT':'C', 'TGC':'C',
    'CAA':'Q', 'CAG':'Q',
    'GAA':'E', 'GAG':'E',
    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
    'CAT':'H', 'CAC':'H',
    'ATT':'I', 'ATC':'I', 'ATA':'I',
    'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    'AAA':'K', 'AAG':'K',
    'ATG':'M',
    'TTT':'F', 'TTC':'F',
    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'TGG':'W',
    'TAT':'Y', 'TAC':'Y',
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    'TAA':'*', 'TGA':'*', 'TAG':'*'
}
    dic_proteina = {}
    for id in dic_max.keys():
        proteina = ''
        seq_orf = dic_max[id]['Sequência']
        codons = re.findall(r'(.{3})',seq_orf)
        for codon in codons:
            if codon in translation_table:
                proteina += translation_table[codon]
            else:
                print(f'Esse codon {codon} não existe')
                sys.exit()
        dic_proteina[id] = proteina
    return dic_proteina


try:
    arquivo = 'C:\\Users\\bccle\\OneDrive\\Documentos\\P2_11244160-CEN0336\\arquivo_fasta.fa'
    if not arquivo.endswith('.fa'):
        raise NotFASTAError('Não é um arquivo fasta')
    with open (arquivo,'r') as arquivo_lido:
        for linha in arquivo_lido:
            if linha.startswith('>'):
                linha = linha.rstrip()
                gene_id = re.search(r'>(\S+)(\s+)?',linha)
                gene_id = gene_id.group(1)
                dic_genes[gene_id] = ''
                sequencia = ''
            else:
                sequencia += linha.rstrip().upper()
                dic_genes[gene_id] = sequencia
        dic_max_orf = max_reading_frame(open_reading_frame(dic_genes))
        dic_proteinas = tranlation(dic_max_orf)
        with open('ORF.fna','w') as saida_fna:
            for id in dic_max_orf.keys():
                saida_fna.write(f">{id}_frame_{dic_max_orf[id]['Frame']}_{dic_max_orf[id]['Começo']}_{dic_max_orf[id]['Fim']}\n")
                for i in range(0, len(dic_max_orf[id]['Sequência']),60):
                    saida_fna.write(dic_max_orf[id]['Sequência'][i:i+60]+'\n')
        with open('ORF.faa','w') as saida_faa:
            for id in dic_max_orf.keys():
                saida_faa.write(f">{id}_frame_{dic_max_orf[id]['Frame']}_{dic_max_orf[id]['Começo']}_{dic_max_orf[id]['Fim']}\n")
                for i in range(0, len(dic_proteinas[id]),60):
                    saida_faa.write(dic_proteinas[id][i:i+60]+'\n')


except IndexError:
    print('Por favor insira o nome de um arquivo')
except IOError:
    print(f'O arquivo {arquivo} não existe.')
except NotFASTAError:
    print(f"O arquivo '{arquivo}' não é do formato FASTA")

