# Learning optical flow for fast MRI reconstruction

[![license](https://img.shields.io/github/license/tschmoderer/iathena.svg?color=blue)](https://github.com/tschmoderer/iathena/blob/master/LICENSE)

This repository contains the MATLAB code for our algorithm **MRIR-DLMC** used to reconstruct undersampled MRI data.

**Article** : [![doi badge](https://img.shields.io/badge/doi%3A-10.1088%2F1361--6420%2Fac164a-blue)](https://dx.doi.org/10.1088/1361-6420/ac164a) [![arxiv badge](https://img.shields.io/badge/arxiv%3A-2004.10464-brightgreen)](https://arxiv.org/abs/2004.10464)

> Reconstructing high-quality magnetic resonance images (MRI) from undersampled raw data is of great interest from both technical and clinical point of views. To this date, however, it is still a mathematically and computationally challenging problem due to its severe ill-posedness, resulting from the highly undersampled data leading to significantly missing information. Whilst a number of techniques have been presented to improve image reconstruction, they only account for spatio-temporal regularisation, which shows its limitations in several relevant scenarios including dynamic data. In this work, we propose a new mathematical model for the reconstruction of high-quality medical MRI from few measurements. Our proposed approach combines - *in a multi-task and hybrid model* - the traditional compressed sensing formulation for the reconstruction of dynamic MRI with motion compensation by learning an optical flow approximation. More precisely, we propose to encode the dynamics in the form of an optical flow model that is sparsely represented over a learned dictionary. This has the advantage that ground truth data is not required in the training of the optical flow term. Furthermore, we present an efficient optimisation scheme to tackle the non-convex problem based on an alternating splitting method. We demonstrate the potentials of our approach through an extensive set of numerical results using different datasets and acceleration factors. Our combined approach reaches and outperforms several state of the art techniques. Finally, we show the ability of our technique to transfer phantom based knowledge to real datasets.     

**Cite** : 

```latex
@article{schmodererDLMCR_2021,
	doi = {10.1088/1361-6420/ac164a},
	url = {https://dx.doi.org/10.1088/1361-6420/ac164a},
	year = {2021},
	month = {aug},
	publisher = {IOP Publishing},
	volume = {37},
	number = {9},
	pages = {095007},
	author = {Timothée Schmoderer and Angelica I Aviles-Rivero and Veronica Corona and Noémie Debroux and Carola-Bibiane Schönlieb},
	title = {Learning optical flow for fast MRI reconstruction},
	journal = {Inverse Problems}
}
```
