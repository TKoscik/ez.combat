# ComBat Harmonization for multi-site/batch data with Dataframe Objects
## Version 1.0.0

****

**Maintainer**: Timothy R. Koscik, <timothy-koscik@uiowa.edu>
**License**: Artistic License 2.0

This is a modification of the ComBat function code from the sva package that can be found at:
**Standard Version**: https://bioconductor.org/packages/release/bioc/html/sva.html
**Prior Modified Version**: https://github.com/Jfortin1/ComBatHarmonization

**Current Modifications**:
The standard version has been modified to add an interface to standard R dataframes as input.
- In the standard version, input data must be a p x n matrix, with p *features in rows* and n *participants in columns*, while in the modified version, input data can be n x p dataframes, with p *features in columns* and n *participants in rows*.
- The modified version adds procedures to allow batch variables to be within the same dataframe as the data as well as methods to select and exclude subsets of features (columns) in dataframes.
- Various formatting changes to the standard version source code have been mode.

****

**References**: If you are using ComBat for the harmonization of multi-site imaging data, please cite the following papers:

|       | Citation     | Paper Link
| -------------  | -------------  | -------------  |
| Original ComBat paper for gene expression array | W. Evan Johnson and Cheng Li, **Adjusting batch effects in microarray expression data using empirical Bayes methods**. Biostatistics, 8(1):118-127, 2007. | [Link](https://academic.oup.com/biostatistics/article/8/1/118/252073/Adjusting-batch-effects-in-microarray-expression) |
| ComBat for multi-site DTI data | Jean-Philippe Fortin, Drew Parker, Birkan Tunc, Takanori Watanabe, Mark A Elliott, Kosha Ruparel, David R Roalf, Theodore D Satterthwaite, Ruben C Gur, Raquel E Gur, Robert T Schultz, Ragini Verma, Russell T Shinohara. **Harmonization Of Multi-Site Diffusion Tensor Imaging Data**. NeuroImage, 161, 149-170, 2017 | [Link](https://www.sciencedirect.com/science/article/pii/S1053811917306948?via%3Dihub#!) |
| ComBat for multi-site cortical thickness measurements | Jean-Philippe Fortin, Nicholas Cullen, Yvette I. Sheline, Warren D. Taylor, Irem Aselcioglu, Philip A. Cook, Phil Adams, Crystal Cooper, Maurizio Fava, Patrick J. McGrath, Melvin McInnis, Mary L. Phillips, Madhukar H. Trivedi, Myrna M. Weissman, Russell T. Shinohara. **Harmonization of cortical thickness measurements across scanners and sites**. NeuroImage, 167, 104-120, 2018  |[Link](https://www.sciencedirect.com/science/article/pii/S105381191730931X)|
