TwoAdvSelfSims.c 

Simulation calculating fixation probability of second beneficial allele, given existing
ben allele at initial frequency p. For use in the study "Linkage and the limits to natural
selection in partially selfing populations".

Simulation uses routines found with the GNU Scientific Library (GSL)
(http://www.gnu.org/software/gsl/)
Since GSL is distributed under the GNU General Public License 
(http://www.gnu.org/copyleft/gpl.html), you must download it 
separately from this file.

This program can be compiled with e.g. GCC using a command like:
gcc TwoAdvSelfSims -lm -lgsl -lgslcblas -I/usr/local/include -L/usr/local/lib TwoAdvSelfSims.c

Then run by executing:
./TwoAdvSelfSims N self rec ha sa hb sb p reps
Where:
- N is the population size
- self is the rate of self-fertilisation
- rec is recombination rate
- ha, hb is dominance at the original, introduced beneficial allele
- sa, sb is selection coefficient of the original, introduced beneficial allele
- p is the initial frequency of the first beneficial allele (when the second is introduced)
- reps is how many times the second allele should FIX before simulation stops 
(the number of actual runs is greater due to stochastic loss of second allele)

Note that haplotypes are defined as:
x1 = ab
x2 = Ab
x3 = aB
x4 = AB

Genotypes defined as:
g11 = g1 = ab/ab
g12 = g2 = Ab/ab
g13 = g3 = aB/ab
g14 = g4 = AB/ab
g22 = g5 = Ab/Ab
g23 = g6 = Ab/aB
g24 = g7 = Ab/AB
g33 = g8 = aB/aB
g34 = g9 = aB/AB
g44 = g10 = AB/AB

Output files are the parameters;
followed by number of times each haplotype fixed;
followed by average total generations elapsed in each case;
Then total number of simulations ran;
Then fixation prob of allele, both unscaled and scaled to unlinked case, 
along with 95% CI intervals for the latter case;
then number of allele fixations.

Note that 'fixation' DIFFERS depending on the inputs of sa, sb.
If sa >= sb (interference case) then 'fixation' counts as fixation of second allele on any genetic background.
If sa < sb (replacement case) then 'fixation' only considers fixation of second allele with neutral haplotype
(I.e. where the 'less fit' neutral allele at locus A fixes, instead of selected allele).
