#!/bin/bash
# Shell script gerador da função de distribuição gaussiana
# Author: Petrucio Ricardo Tavares de Medeiros
# Date: February, 2016.
# Attention: The installation of Opencv need OpenCV, GNUPLOT and foveatedFeatures.

# Compilando os arquivos
make

# Execução do arquivo
make generatorExecutavel

# Gerando o plot
gnuplot generator.gnu
