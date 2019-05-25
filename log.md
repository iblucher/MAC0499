# LOG TCC

## MAIO
### TODO
- [x] Ler todos os dados FITS: ambas as estrelas reais, dados do Sol e estrelas padrão
- [x] Fazer a divisão das 4 estrelas e plotar espectro resultante
- Ler paper do TAPAS e entender como a transmissão atmosférica é modelada
- Ler paper do MolecFit e entender o que faz e como faz (ver código no github)

### DÚVIDAS
- Fluxo do Sol e de Arcturus com valores negativos. Como lidar com isso? (vi na internet sobre zero-clip nos valores negativos, mas melhor confirmar)

### OBSERVAÇÕES
- Dados UVES em PrimaryHDU, com os dados de fluxo como "data", e dados que auxiliam montagem do comprimento de onda no "header" (CRVAL1 e CDELT1)
- Dados ardata.fits em BinTableHDU, 4 colunas com comprimento de onda (já o mesmo intervalo pra todas as estrelas) e os fluxos do Sol, Arcturus e a estrela telúrica. 
- Divisão da estrela HD110379 gerou artifícios (pico de fluxo > 1.0) e pequenas mudanças visuais
- Divisão da estrela HD186791 gerou artifícios (pico de fluxo > 1.0) e mais difícil ver mudanças visuais
- Divisão do Sol gerou artifícios (pico de fluxo > 1.0), mas é muito difícil de observar devido à quantidade de linhas no espectro
- Divisão de Arcturus gerou artifícios (picos de fluxo >> 1.0). Divisão gerou fluxo muito grande, principalmente perto do comprimento de onda 9000Å. Pode ter por causa do zero-clip dos valores negativos do fluxo.

### DONE
- Funções de manipulação dos arquivos estelares e de plot do espectro
- Divisão dos espectros UVES
- Divisão dos espectros ardata