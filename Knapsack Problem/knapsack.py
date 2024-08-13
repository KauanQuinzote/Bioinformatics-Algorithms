
def main():
    # Lista de itens, cada um com um valor e peso
    itens = [
        {'value': 5, 'weight': 5},
        {'value': 4, 'weight': 6},
        {'value': 7, 'weight': 8},
        {'value': 7, 'weight': 4}
    ]

    def knapsack(items, capacity):
        # Extrai listas de pesos e valores dos dicionários de itens
        weights = [item['weight'] for item in items]
        values = [item['value'] for item in items]
        
        n = len(weights)  # Número de itens
        
        # Cria a tabela DP com dimensões (n+1) x (capacity+1)
        # Inicializa todas as entradas com 0
        table = [[0] * (capacity + 1) for _ in range(n + 1)]

        # Preenche a tabela DP
        for i in range(1, n + 1):
            for j in range(1, capacity + 1):
                if weights[i - 1] <= j:
                    # Inclui o item i ou não inclui
                    # Compara o valor de incluir o item i com o valor de não incluir o item i
                    table[i][j] = max(table[i - 1][j], values[i - 1] + table[i - 1][j - weights[i - 1]])
                else:
                    # Não pode incluir o item i, então o valor é o mesmo do subproblema anterior
                    table[i][j] = table[i - 1][j]

        # A solução final está na última célula da tabela, que representa
        # a capacidade máxima com todos os itens considerados
        return table[n][capacity]

    # Define a capacidade máxima da mochila
    capacity = 13
    
    # Chama a função knapsack e imprime o valor máximo que pode ser obtido
    answer = knapsack(itens, capacity)
    
    with open('ans.txt', 'w') as ans:
        ans.write("Resposta: {}".format(answer))
        

# Chama a função principal para executar o código
main()
