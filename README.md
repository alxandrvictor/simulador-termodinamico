# Este projeto consiste em um simulador termodinâmico desenvolvido em **Fortran 90/95** para o cálculo e dimensionamento de sistemas de refrigeração por compressão de vapor. 

Trabalho final desenvolvido para o curso **"Métodos Numéricos: Introdução Do Cálculo à Simulação Numérica"**, ministrado pelo **Prof. Dr. Rafael Gabler Gontijo** (L2C / Canal Ciência e Brisa).

## Tecnologias e Métodos Numéricos

O código resolve o balanço de massa e energia do ciclo utilizando:
* **Newton-Raphson:** Solver principal para sistemas de equações não-lineares.
* **Eliminação de Gauss:** Resolução do sistema linear a cada iteração.
* **Diferenças Finitas:** Aproximação da Matriz Jacobiana de sensibilidade.
* **Simulação da Compressão:** Modelagem do estado de superaquecimento e eficiência isentrópica.
* **Monte Carlo:** Geração de perturbações estocásticas na carga térmica.

## Configuração da Simulação

O simulador é versátil e permite a análise de diferentes cenários. Para obter resultados específicos, o usuário deve configurar o código antes da compilação:

### 1. Escolha do Caso (Fluido)
No início do código, é necessário selecionar o perfil de operação desejado:
* **Perfil Doméstico:** Utiliza o fluido **R134a**, típico de geladeiras residenciais.
* **Perfil Comercial:** Utiliza o fluido **R404A**, voltado para balcões frigoríficos e câmaras frias.

### 2. Alteração de Variáveis de Entrada
Para verificar como o sistema reage a diferentes condições (análise de sensibilidade), o usuário pode modificar as variáveis de contorno diretamente no arquivo fonte:
* **Q (Carga Térmica):** Altera a demanda de resfriamento inicial.
* **Te (Temperatura de Evaporação):** Define a temperatura no trocador de calor interno.
* **Parâmetros de Monte Carlo:** Ajuste da amplitude do desvio padrão para simular a intensidade do "susto" térmico no sistema.

**Nota:** O valor de **Tc (Temperatura de Condensação)** e a vazão mássica são calculados iterativamente pelo solver, adaptando-se às variações estocásticas impostas durante as 20 iterações.

> **Nota:** A mudança nestes parâmetros impactará diretamente o cálculo da vazão mássica ($x_4$) e o deslocamento final do compressor.

## Como Executar

1. **Configuração:** Abra o arquivo `.f90` e ajuste as variáveis de entrada e o fluido desejado.
2. **Compilação:** Certifique-se de ter o `gfortran` instalado.
   ```bash
   gfortran -o simulador main.f90

## Análise de Resultados

Ao rodar a simulação, o console exibirá um log detalhado contendo:

* **Monitoramento em Tempo Real:** Série temporal de 20 iterações demonstrando a convergência do solver frente às variações de carga.
* **Estabilidade Numérica:** Validação do uso do fator de relaxação e da sub-rotina de Eliminação de Gauss frente a perturbações.
* **Relatório de Eficiência:** Cálculo do COP e Eficiência Exergética calibrados para padrões reais (35% a 45%).
* **Dimensionamento de Hardware:** Cálculo do deslocamento volumétrico médio (cm^3/rev) para seleção técnica de compressores.


## Licença

Este projeto está licenciado sob a Licença MIT - consulte o arquivo `LICENSE` para mais detalhes.

## Créditos e Referências

Projeto desenvolvido como requisito para a conclusão do curso de **Métodos Numéricos**, sob orientação do **Prof. Dr. Rafael Gabler Gontijo**.

* **Canal Ciência e Brisa:** [YouTube](https://www.youtube.com/@cienciaebrisa)
* **L2C Treinamentos:** [Portal L2C](https://l2ctreinamentos.com.br/)
* **Instrutor:** Dr. Rafael Gabler Gontijo

---
*Desenvolvido como portfólio de Engenharia e Simulação Numérica.*
