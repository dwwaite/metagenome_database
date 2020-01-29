# The *metagenome_database* project

This is an attempt to create a unified system for storing metagenomic data. The idea will be to take the results of the following bioinformatic processes:

1. Metagenome assembly
1. Gene and non-coding sequence prediction
1. Gene annotation
1. Metagenome binning (creating MAGs)
1. MAG quality assessment
1. MAG taxonomy assignment

And store all data in a single file, with timestamps and versioning data (for example, if gene annotations are changed through time).

Once data is stored, users will be able to extract slices of their data for analysis. The hope for this system is that it will simplify the workflow for students and collaborators of the Handley Lab, so that we avoid the issue where people have different versions of files, scattered across different places.

### Under the hood

The system exists as a series of `python3` scripts, using `SQLite` to create a relational database with the following ERD:

![](https://github.com/dwwaite/metagenome_database/blob/master/docs/figures/handleydb_erd.png)

*(Apologies for the fact that this is clearly a screenshot from `Microsoft Access`. I'll make a nicer one soon!)*
