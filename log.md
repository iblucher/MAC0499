# LOG TCC

## MAIO
### TODO
- [x] Ler todos os dados FITS: ambas as estrelas reais, dados do Sol e estrelas padrão
- [x] Fazer a divisão das 4 estrelas e plotar espectro resultante
- [x] Ler artigo "Using a model for telluric absorption in full-spectrum fits" 
- [x] Ler paper do TAPAS
- [x] Ler paper do MolecFit

### DÚVIDAS
- Fluxo do Sol e de Arcturus com valores negativos. Como lidar com isso? (vi na internet sobre zero-clip nos valores negativos, mas melhor confirmar)

### OBSERVAÇÕES
- Dados UVES em PrimaryHDU, com os dados de fluxo como "data", e dados que auxiliam montagem do comprimento de onda no "header" (CRVAL1 e CDELT1)
- Dados ardata.fits em BinTableHDU, 4 colunas com comprimento de onda (já o mesmo intervalo pra todas as estrelas) e os fluxos do Sol, Arcturus e a estrela telúrica. 
- Divisão da estrela HD110379 gerou artifícios (pico de fluxo > 1.0) e pequenas mudanças visuais
- Divisão da estrela HD186791 gerou artifícios (pico de fluxo > 1.0) e mais difícil ver mudanças visuais
- Divisão do Sol gerou artifícios (pico de fluxo > 1.0), mas é muito difícil de observar devido à quantidade de linhas no espectro
- Divisão de Arcturus gerou artifícios (picos de fluxo >> 1.0). Divisão gerou fluxo muito grande, principalmente perto do comprimento de onda 9000Å. Pode ter por causa do zero-clip dos valores negativos do fluxo.
- Todas as leituras até o momento indicam que o modelo usado para a transmissão atmosférica é o LBLRTM. Vale a pena entender como funciona o código? Ou será que vamos usar o próprio TAPAS para fazer várias simulações e tirar medidas em cima disso? Definir melhor o experimento a ser implementado. 
- Todos os modelos de transmissão atmosférica vistos até agora usam como base o LBLRTM e a base de dados de linhas espectrais HITRAN para modelar a line-spread function e espectro da atmosfera.

### DONE
- Funções de manipulação dos arquivos estelares e de plot do espectro
- Divisão dos espectros UVES
- Divisão dos espectros ardata
- Leituras dos papers