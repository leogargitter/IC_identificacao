# IC_identificacao
Arquivos relacionados a minha IC de 2020-2021


\section*{Resumo}

Amplificadores de potência (PAs - \textit{Power Amplifiers}) são dispositivos que tem como função a amplificar um sinal com o objetivo de obter a potência desejada, esses são os dispositivos com maior consumo de potência dentro dos sistemas de comunicação sem fio além de operarem de forma não linear. Para contornar essas desvantagens dos PAs existem técnicas como a pré-distorção digital (DPD - \textit{Digital Pré Distortion}) e uma forma de se definir o bloco de pré distorção são modelos polinomiais como o baseado em polinômios com memória (MP - \textit{Memory Polynomial}). O foco desse trabalho será na fase de identificação dos modelos comportamentais, ou seja, a etapa onde são definidos os elementos necessários para a implementação do modelo. Em trabalhos anteriores foram trabalhadas as implementações dos modelos por tabelas de busca e usando matrizes reais, portanto o objetivo desse trabalho é estudar métodos de identificação que definam os valores a serem inseridos nas tabelas de busca ou matrizes reais sem o conhecimento prévio dos coeficientes.

\section*{Palavras-chave}

modelagem comportamental; amplificadores de potência; identificação de modelos


\section*{Pesquisa Bibliográfica}

Amplificadores de potência são dispositivos que tem como função a amplificação da potência de um determinado sinal de entrada de maneira que o sinal de saída tenha a potência necessária, esses dispositivos são fundamentais para os sistemas de comunicação sem fio. Apesar da grande importância nas comunicações sem fio, os PAs são os dispositivos com maior consumo de potência além de operarem de maneira não linear, podendo gerar distorções no sinal amplificado. Uma maneira de contornar o problema da não linearidade é a técnica pré-distorção digital que consiste em passar o sinal de entrada por um bloco de pré distorção com uma não linearidade complementar à do PA, de forma que essa pré distorção em cascata com o PA resulte em um ganho linear, outra vantagem da linearização é o aumento da eficiência do PA.

Nesse trabalho, para a aplicação das técnicas de identificação, será considerado o modelo Memory Polynomial descrito pela equação a seguir:

\begin{equation} \label{MP}
\tilde{y}(n) =  \sum_{p=1}^P \sum_{m=0}^M \tilde{b}_{mp} \: \tilde{x}(n-m) |\tilde{x}(n-m)|^{p-1}
\end{equation}, onde $P$ representa a ordem polinomial, $M$ é a ordem do efeito de memória, $\tilde{x}$ representa o sinal de entrada do PA, $\tilde{b}_{mp}$ são os coeficientes e $\tilde{y}$ será o sinal de saída estimado. Esse modelo é baseado na Série de Volterra no tempo discreto que tem como característica ser linear nos parâmetros. Isso significa que os coeficientes do modelo podem ser obtidos através de um sistema de identificação linear como o Método dos Mínimos Quadrados (MMQ). Uma característica importante para a reprodução desse modelo é que as multiplicações devem ser feitas sempre entre elementos em um mesmo instante de tempo.

Esse modelo também pode ser escrito da seguinte forma:

\begin{equation} \label{MP_fpol}
    \tilde{y}(n) =  \sum_{m=0}^M \tilde{x}(n-m) f^{pol}_m(|\tilde{x}(n-m)|) 
\end{equation}, onde $f^{pol}_m$ podem ser chamadas de funções polinomiais unidimensionais. A abordagem de implementação usando tabelas de busca (LUTs - \textit{Look Up Tables}) substituí essas funções polinomiais por tabelas de busca, dessa forma reduzindo a quantidade de operações necessárias para a implementação do modelo.

Além de métodos como o MMQ em um artigo usado como base para este trabalho foi proposta outra técnica de identificação. A técnica proposta tem como foco a implementação do modelo por tabelas de busca, encontrando os valores para o preenchimento da tabela sem conhecimento dos coeficientes $\tilde{b}_{mp}$. Para a aplicação desse método de identificação primeiramente devemos montar a matriz $X_{LUT}$:


\begin{equation}

    X_{LUT} = 
    \begin{bsmallmatrix}
        1blocodeQcolunasdalinha1 & 2blocodeQcolunasdalinha1 & \dots & (M+1)-ésimoblocodeQcolunasdalinha1\\
         & \vdots & & \\
        1blocodeQcolunasdalinha(n-M) & 2blocodeQcolunasdalinha(n-M) & \dots & (M+1)-ésimoblocodeQcolunasdalinha(n-M)\\
    \end{bsmallmatrix}  
    
    
\end{equation}

A matriz acima foi representada em blocos de Q colunas e 1 linha, o valor de Q é definido previamente e corresponde ao número de entradas das tabelas de busca. Essa matriz também pode ser escrita da seguinte forma:

\begin{equation} \label{x_lut}
    X_{LUT} = 
    \begin{bmatrix}
    \tilde{x}(1+M) f_0^{LUT}(|\tilde{x}(1+M)|) & \tilde{x}(M) f_1^{LUT}(|\tilde{x}(M)|) & \dots &  \tilde{x}(1) f_M^{LUT}(|\tilde{x}(1)|)\\
     & \vdots & & \\
    \tilde{x}(n) f_0^{LUT}(|\tilde{x}(n)|) & \tilde{x}(n-1) f_1^{LUT}(|\tilde{x}(n-1)|) & \dots &  \tilde{x}(n-M) f_M^{LUT}(|\tilde{x}(n-M)|)
    \end{bmatrix}
\end{equation}

Portanto para a montagem da matriz precisamos considerar saída interpolada da LUT:

\begin{equation}
    f^{pol}_m (|\tilde{x}(n-m)|)\sim f^{LUT}_m (|\tilde{x}(n-m)|) = \tilde{s}_{m(q-1)} + \left[ \frac{\tilde{s}_{m(q)} - \tilde{s}_{m(q-1)}}{e_{m(q)} - e_{m(q-1)}} \right] [|\tilde{x}(n-m)| - e_{m(q-1)}]
\end{equation}

Para estimar essa saída a condição $e_{m(q-1)}<|\tilde{x}(n-m)|<e_{mq}$ com $q$ variando entre 2 e Q deve ser satisfeita. A partir da equação anterior verifica-se que a saída interpolada de uma LUT depende das informações de duas posições consecutivas da LUT. Assim para preencher cada bloco da matriz $X_{LUT}$ consideramos as seguintes equações: 

\begin{equation}
    Valor_{(q-1)_{esimacoluna}} = \tilde{x}(n-m) \left[ 1 - \frac{[|\tilde{x}(n-m)| - e_{m(q-1)}]}{e_{m(q)} - e_{m(q-1)}} \right]
\end{equation}

\begin{equation}
    Valor_{(q)_{esimacoluna}} = \tilde{x}(n-m) \left[ \frac{[|\tilde{x}(n-m)| - e_{m(q-1)}]}{e_{m(q)} - e_{m(q-1)}} \right]
\end{equation}

Com a matriz $X_{LUT}$ montada é apenas aplicar o MMQ usando essa matriz para se obter as saídas das LUTs:

\begin{equation}
    s = (X_{LUT}^* X_{LUT})^{-1} (X_{LUT}^* Y) =
    \begin{bmatrix}
    \tilde{s}_{01} \\
    \vdots \\
    \tilde{s}_{0Q} \\
    \vdots \\
    \tilde{s}_{M1} \\
    \vdots \\
    \tilde{s}_{MQ}
    \end{bmatrix}
\end{equation}

\section*{Metodologia}

Esse trabalho tem o intuito de investigar com mais detalhes a etapa de identificação de modelos comportamentais para amplificadores de potência. Para fins práticos será considerado o modelo Memory Polynomial. No desenvolvimento do código que representa o modelo é utilizada a linguagem de programação de uso geral Python 3, em conjunto com as bibliotecas Numpy e SciPy. 

Os dados fornecidos foram divididos em dois conjuntos, um com o propósito da extração dos coeficientes e outro para a validação do modelo. Esses dados foram obtidos a partir de um PA baseado em um amplificador do tipo LDMOS excitado por uma portadora com frequência de 2GHz, modulada por um sinal WCDMA de largura de banda igual a 3,84 MHz e com frequência de amostragem de 30,72 MHz.

A métrica utilizada para avaliar os resultados do modelo é o NMSE (\textit{Normalized Mean Square Error} - Erro Quadrático Médio Normalizado), que é obtido a partir da equação \eqref{NMSE}:

\begin{equation} \label{NMSE}
NMSE = 10 log_{10} \frac{\sum^N_{n=1} |e(n)|^2}{\sum^N_{n=1} |y_{real}(n)|^2}
\end{equation}, onde $e(n)$ é o erro entre a saída esperada e a saída estimada pelo modelo em um instante $n$, $y_{real}(n)$ é o valor da amostra do sinal de saída no instante $n$ e $N$ representa o número total de amostras. Quanto menor for o valor calculado do NMSE melhores são os resultados. Também será considerado o número de operações necessárias para a implementação de cada modelo.


\section*{Resultados Parciais Alcançados}