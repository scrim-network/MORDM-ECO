#Many-objective robust decision making for managing an ecosystem with a deeply uncertain threshold response

README files last updated by Kelsey Ruckert, May 26, 2016.

##Citation
This code is intended to accompany the results of the paper listed below. _Please cite this source when using this code._
>Singh R., P.M. Reed, and K. Keller (2015). Many-objective robust decision making for managing an ecosystem with a deeply uncertain threshold response. _Ecology and Society_ **20**(3): 12. doi: [10.5751/ES-07687-200312.](http://www.ecologyandsociety.org/vol20/iss3/art12/)

##Overview
This code introduces a framework that allows decision makers or users to pose multiple objectives, explore trade-offs between strategies that are robust to to deep uncertainties. The framework, referred to as many-objective robust decision making (MORDM), employs multiobjective evolutionary search to identify trade-offs between strategies, re-evaluates their performance under deep uncertainty, and uses interactive visual analytics to support the selection of robust management strategies.

The code demonstrates MORDM on a stylized decision problem posed by the management of a lake in which surpassing a pollution threshold causes eutrophication. The MORDM framework enables the discovery of strategies that balance multiple preferences and perform well under deep uncertainty. This decision analytic framework allows the decision makers to select strategies with a better understanding of their expected trade-offs (traditional uncertainty) as well as their robustness (deep uncertainty).

##Notes
This code requires:

- C++ _(multi-objective optimization)_
- R _(visualizing output)_

The files on this repository include only the main data, files with output policy vectors, time evolution of phosphorus states, and other small files. Github doesn't support files exceeding 100 MB.

A number of code files need the parallel version of BORG to execute while other files are independent of BORG as they just post-process the output from BORG.  In addition, users may have to make minor edits such as changing the file paths.

Note: Not all of the visualization files are included because one or two need Aerovis outputs. Although, most of the main result figures are included.

##Contact
Riddhi Singh: <riddhi@iith.ac.in>

