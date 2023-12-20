#!/usr/bin/env python3

#define as variaveis como 0
total = 0
contador_notas = 0

while contador_notas < 10:#loop para armazenar as 10 notas no período
    try:
        nota = float(input("Digite uma nota entre 0 e 10: "))#transforma a nota em tipo float caso seja possivel
        if nota<=10 and nota>=0:#define um limite para as notas
            total += nota
            contador_notas += 1
        else:
            print("Nota inválida. Tente novamente.")#caso o limite não seja respeitado
            continue
    except ValueError:
        print('Digite apenas número inteiros ou flutuantes!')#mensagem de erro caso não seja digitado números
media_disciplinas = total/10#calcula a média
print(f'A média da disciplina é {media_disciplinas}')#imprime a média