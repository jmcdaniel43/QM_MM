#!/bin/bash

./gdma < psi4_dma_datafile > temp.dma
./parse_psiDMA.pl temp.dma > temp_format.dma
./mpfit temp_format.dma
