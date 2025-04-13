Reviewer: 2
3) on page 6, it is mentioned that 11,253, 377, 128, 95, etc. 
features were generated through MNA database corresponding to Fig 2A, 2B, 
so it would be nice to know 
how many features are formed by experimental data and how many features are formed by theoretical value data in the generated features, 
and how many features are categorized into Level A, B, and C as mentioned in the manuscript.

A: 
E_IS_Counts.py 统计CCML 和非CCMSL的数量，其实统计的是 > 0.7的数量， ≥0.7的数量会稍多

ClassCount.py 统计 Level A、B1,B2、C1、C2的数量


Section 2.2 区分一下experimental 和 in-silico的注释 
To clarify the distribution of features their categorization, we provide the following tables. All code and data used are available on GitHub
Table 1. Generated features from experimental and in-silico libraries
Hits	Experimental	In-silico
Top1	16	9
Top5	51	44
Top20	125	252
Table 2. Feature Categorization into Annotation Levels
Level	A	B1	B2	C1	C2
Counts	0	62	32	5	0
