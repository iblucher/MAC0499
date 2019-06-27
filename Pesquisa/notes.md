# PESQUISA TCC

## REUNIÃO 10/06

### TAPAS
- "TAPAS computes the atmospheric transmission in the line-of-sight to the target indicated by the user"
-  O usuário precisa fornecer o horário, a localização no solo e as coordenadas equatoriais do alvo ou o ângulo zênite da linha de visada. O perfil atmosférico daquele lugar naquele momento é recuperado (temperatura, pressão, umidade, quantidade de ozônio) da base de dados do ETHER.  
-  A transmissão atmosférica é calculada por um software de LBLRTM para vários gases como: : O2, H2O, O3, CO2.
- "Importantly, our advice for an optimal use of
the TAPAS tool, for the second and third purposes, is the downloading
of the H2O and O2 absorptions spectra separately, which allows refined adjustments and taking into account differential
behavior of the two species at the time of the observations."
- "In order to compute the atmospheric transmission,
TAPAS makes use of LBLRTM (Line-By-Line Radiative
Transfer Model) (Clough and Iacono, 1995) code and the
2012 HITRAN spectroscopic data base HITRAN (this is a
high-resolution transmission molecular absorption database,
www.cfa.harvard.edu/hitran/, (Rothman et al., 2011) for the
HITRAN 2008 )."
    - Existe uma API disponível para o HITRAN implementada em Python. Se usarmos o LBLRTM podemos fazer uso dela.
- "One thing that is not accounted
for in the computation is the wavelength shift due to Doppler
effect induced by the atmospheric wind" -- essa correção foi introduzida depois do uso do TAPAS no artigo (Fig. 3)
- "The mathematically correct method
would be to deconvolve the observed spectrum from the ILSF
with a considerable wavelength oversampling, and then divide
by the transmittance at the highest possible resolution. In practice,
however, even with high SNR spectra, the data noise would
not allow this mathematically correct method. This shortcoming
of the division method is obviously more important with strong
absorption features than with small ones." -- não entendi porque seria necessário esse "deconvolve" para o método que é o melhor matematicamente.
- Depois de introduzir o método usado para a computação transmissão da atmosfera, o paper faz as divisões de espectros estelares pelo resultado do TAPAS separadamente para diferentes moléculas (e suas sub-faixas no espectro) e observa os resultados nelas
- Usos para o TAPAS no paper: correct the observed spectrum from telluric absorption.
For H2O, it may need some manipulation outside of the TAPAS
environment,or several calls to TAPAS to get a good fit for all
the lines (pode ser interessante se usarmos o TAPAS assim).
- TAPAS should be most useful on cool stars, where stellar
lines are narrow, like O2 and H2O lines. Since these stars are
heavily used in the search for exo-planets, it is particularly interesting
to use the atmospheric lines for a wavelength standard
(Figueira et al., 2010). Also, if well corrected with TAPAS,
some telluric contaminated regions which are presently discarded
from the exo-planet search through radial velocity variations
could be added to the analysis for a better retrieval of Vr
changes. -- isso pode ser interessante de mencionar no TCC 

### LBLRTM
- Implementação original em Fortran
- Existem wrappers em Python para o LBLRTM 

### MolecFit
- "Molecfit combines a publicly available radiative transfer code, a molecular line database, atmospheric profiles, and
various kernels to model the instrument line spread function. The atmospheric profiles are created by merging a standard atmospheric
profile representative of a given observatory’s climate, of local meteorological data, and of dynamically retrieved altitude profiles for
temperature, pressure, and humidity."
- "The availability of such a general tool for telluric absorption correction may improve future observational and analysing strategies, as
well as empower users of archival data."
- "Molecfit is a versatile tool for modelling and correcting telluric
absorption lines. In its most common use, molecfit fits a
spectrum of the transmission of the Earth’s atmosphere over narrow,
user-selected ranges (inclusion regions) of the spectrum of a
science target." -- usa espectro input do usuário para saber as regiões sobre as quais deve modelar a transmissão atmosférica e devolve tanto o espectro da atmosfera quanto o espectro corrigdo (input / modelo)
- Possui uma interface gráfica
- Passo as passo:
    - 1) Retrieval of atmospheric profiles used as input for LBLRTM
    - 2) LBLRTM que simula emissão e transmissão atmosférica
    - 3) Base de dados de espectroscopia molecular
    - 4) Calculadora de radiação de corpo cinza para contabilizar o instrumento e o telescópio
    - 5) Escolha de line spread functions
    - 6) Fitting algorithm com vários ajustes e calibrações

### "Using a model for telluric absorption in full-spectrum fits"
- https://arxiv.org/abs/1312.2450
- "We present a new method, where we use a model
for the transmission of the Earth’s atmosphere in a full-spectrum fit, which determines
the parameters for the stellar and Earth’s atmosphere simultaneously"
- "A/B type stars show strong hydrogen absorption
features that need to be taken into account when removing the tellurics. The latter item can be
addressed by interpolating linearly over the hydrogen lines, but this only works well if one is not
particularly interested in these regions of the spectrum"
- "The need of a standard star can be avoided completely when using theoretical models for
the atmospheric transmission, which for this work have been computed using LBLRTM (LineBy-Line
Radiative Transfer Model)"
- Usam um modelo para a transmissão da atmosfera e levam em conta moléculas mas proeminentes (como H2O). Seria legal entendermos como funciona a computação de G(x) e se é útil para o trabalho.

## REUNIÃO 27/06

### Organização monografia

- Introdução
    - Na introdução contextualize o tema, apresente as motivações, os objetivos, enfatize a proposta/abordagem/elemento central trabalhado, liste as principais contribuições ou resultados (se houverem), e descreva a organização do texto.
    - Tema: filtragem de ruídos telúricos em sinais astronômico
    - Contextualização
    - Estrutura do trabalho
- Desenvolvimento
    - Na parte de desenvolvimento, devem ser apresentados todos os elementos necessários para o entendimento do tema estudado. Por exemplo, definições, conceitos fundamentais, técnicas e métodos, algoritmos, entre outros. Esses elementos podem ser apresentados por meio da descrição explícita, ou por meio de citações bibliográficas. Quando for o caso, devem ser apresentados também os experimentos, os resultados e as discussões pertinentes.
    - Fundamentação teórica
        - Astronomica
        - Compuatacional
    - Mais desenvolvimento
    - Implementação
        - Ferramentas
    - Resultados

- Conclusão
    - Na conclusão, deve-se "fechar" o trabalho. Para tanto, pode-se listar as conclusões do seu trabalho, levando-se em conta as motivações e objetivos listados na introdução. Pode-se também incluir possíveis desdobramentos do seu trabalho.
    - Retomar assunto e tema e reafirmar importância da pesquisa
    - Repostas:
        - aos objetivos
        - às hipóteses
        - aos problemas da pesquisa
