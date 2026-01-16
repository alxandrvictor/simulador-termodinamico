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
* **Perfil Doméstico:** Utilizando o fluido **R134a**, típico de geladeiras residenciais.
* **Perfil Comercial:** Utilizando o fluido **R404A**, voltado para balcões frigoríficos e câmaras frias.
* **Perfil Comercial:** Utilizando o fluido **R134a**, também voltado para balcões frigoríficos e câmaras frias, mas com propriedades termodinâmicas diferentes.
* 
### 2. Alteração de Variáveis no Código-Fonte
Para modificar as condições de contorno de engenharia, localize o bloco `! --- 2. CONFIGURAÇÃO INICIAL ---` no arquivo `simulador_termodinamico.f90`:

| Variável | Função | Localização |
| :--- | :--- | :--- |
| `Te` | Temperatura de Evaporação ($^\circ C$) | Dentro do `select case` do perfil |
| `Q_base` | Carga Térmica nominal ($kW$) | Define a capacidade do sistema |
| `n_is` | Eficiência Isentrópica | Modelagem real do compressor |
| `tempo_total` | Número de passos da simulação | No loop principal de tempo |

> **Nota:** Para ajustar a intensidade da simulação de Monte Carlo, modifique os multiplicadores de `Q_atual` no bloco `! --- 5. LOOP PRINCIPAL ---`.

---

## 📊 Análise de Resultados

O simulador entrega uma saída detalhada que permite validar o dimensionamento do sistema:

1.  **Convergência Numérica:** Monitoramento do erro residual e estabilidade da Jacobiana.
2.  **Certificação Energética:** Classificação do sistema (Classe A a D) baseada na **Eficiência de 2ª Lei (Exergética)**.
3.  **Dimensionamento de Hardware:** Cálculo do deslocamento volumétrico médio ($cm^3/rev$) para seleção de compressores reais (Ex: Embraco, Bitzer).
4.  **Análise Econômica:** Estimativa de custo mensal baseada em tarifas de energia configuráveis.

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
