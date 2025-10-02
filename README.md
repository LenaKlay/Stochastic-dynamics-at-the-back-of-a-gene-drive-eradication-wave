# Stochastic dynamics at the back of a gene drive eradication wave

## Introduction

Artificial gene drive is a genetic engineering technology that could be used for the control of natural populations. Gene drive alleles bias their own transmission and can therefore spread in a population within a relatively small number of generations, even if they are deleterious. Understanding the potential outcomes of this technology, including the modification and/or the eradication of a natural population, is essential before real-world applications are considered. Here is the code supporting the analyses presented in the paper “Stochastic dynamics at the back of a gene drive eradication wave”.

## Authors

This code was written by Léna Kläy. The project was carried out in collaboration with Vincent Calvez, Florence Débarre and Léo Girardin.

## Abstract of the article

Gene drive alleles bias their own inheritance to offspring. They can fix in a wild-type population in spite of a fitness cost, and even lead to the eradication of the target population if the fitness cost is high. However, this outcome may be prevented or delayed if areas previously cleared by the drive are recolonised by wild-type individuals. Here, we investigate the conditions under which these stochastic wild-type recolonisation events are likely and when they are unlikely to occur in one spatial dimension. More precisely, we examine the conditions ensuring that the last individual carrying a wild-type allele is surrounded by a large enough number of drive homozygous individuals, resulting in a very low chance of wild-type recolonisation. To do so, we make a deterministic approximation of the distribution of drive alleles within the wave, and we split the distribution of wild-type alleles into a deterministic part and a stochastic part. Our analytical and numerical results suggest that the probability of wild-type recolonisation events increases with lower fitness of drive individuals and with smaller local carrying capacity. Numerical simulations show that these results extend to two spatial dimensions. The role of the migration rate however, is less clear but has a lower impact. We further demonstrate that, in the event of wild-type recolonization, the probability of subsequent drive reinvasion decreases with smaller values of the intrinsic growth rate of the population. Overall, our study paves the way for further analysis of wild-type recolonisation at the back of eradication traveling waves.

## Contents

This repository contains various folders:

1) `Functions` contains the code to run the simulations (.py), as well as a `README.rmd` file detailing each function,

2) `Migale` contains the code to run the heaviest simulations on the cluster Migale (INRAE, doi: 10.15454/1.5572390655343293E12) as well as some outputs,

3) `Outputs` stores the results of the simulations.

