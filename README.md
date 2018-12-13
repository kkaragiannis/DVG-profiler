# DVG profiler
___

### Description
The defective viral genomes (DVG) profiler is a post sequence alignment processing algorithm that creates the
best discontinuous alignments of each read. The algorithm is relying on the fact that sensitive Smith
Waterman alignment is performed and that all the possible hits of each short read are returned. DVG
profiler will search all possible combinations of all reported hits and will select those that better explain the reads. 
The algorithm, also resolves conflicts of multiple solutions by maximizing the score of the total alignment. 
It is considering all possible combination of direction and relative position of split reads and as such can detect 
insertions, deletion, 5' and 3' copy backs.


### Requirements

The algorithm has been implemented as part of [HIVE][hive_git-repo-url] platform, a robust infrastructure for next-generation sequence (NGS) data analysis co-developed by Food and Drug Administration and George Washington University
[![HIVE][hive_logo]][hive_git-repo-url] 

DVG profiler itself is open source with a [public repository][dvgp_git-repo-url] on GitHub.

### Installation

DVG profiler includes the required code base for HIVE and is installed as an application during HIVE installation. 
HIVE installation instructions can be found [here][hive_readme].


   [dvg_pipeline]: <https://raw.githubusercontent.com/kkaragiannis/DVG-profiler/master/doc/images/Resized_overview_final_fig.png>
   [hive_readme]: <https://github.com/kkaragiannis/hexahedron/blob/master/doc/HIVE_README.md>
   [hive_logo]: <https://raw.githubusercontent.com/FDA/fda-hive/master/doc/images/hive_logo.png>
   [hive_git-repo-url]: <https://github.com/FDA/fda-hive>
   [hh_logo]: <https://raw.githubusercontent.com/kkaragiannis/hexahedron/master/doc/images/hexahedron_logo.png>
   [dvgp_git-repo-url]: <https://github.com/kkaragiannis/DVG-profiler>
