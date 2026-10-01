# MinimalResolution

This is a p=3 fork of [Guozhen Wang's code](https://github.com/pouiyter/MinimalResolution)
for computing the Adams-Novikov E<sub>2</sub> page for the sphere (originally
at p=2), using the algebraic Novikov spectral sequence. See
[this repository](https://github.com/ebelmont/ANSS_data) for sample data and an
explanation of how to interpret the data files output by the program. The
rest of this README is Guozhen's original instructions.

**This branch computes Ext(BP<sub>*</sub>/3<sup>40</sup>), not
Ext(BP<sub>*</sub>).** Coefficients are integers mod 3<sup>40</sup>
(`Z3.cpp`). On the `prime3` branch that is an approximation, which breaks
down once the computation reaches algebraic Novikov filtration 40 and limits
that branch to stem + ANSS filtration < 160. Here it is the intended
coefficient ring, so the computation is exact in every filtration.

The long exact sequence of
0 → BP<sub>*</sub> —3<sup>40</sup>→ BP<sub>*</sub> → BP<sub>*</sub>/3<sup>40</sup> → 0
shows what changes: almost every class of Ext(BP<sub>*</sub>) appears twice
in Ext(BP<sub>*</sub>/3<sup>40</sup>), once as its image and once as a
"top-cell" class one Ext degree lower. So the output is roughly twice the size
of Ext(BP<sub>*</sub>), with v<sub>0</sub>-towers of height up to 40, and that
is the correct answer for this target.

The second argument L to `mr_BP` caps only the ANSS filtration s, so the
output is complete for s ≤ L − 1 in every stem up to the maximum degree.

******************************************************************************************************

The algorithm is explained in the pdf file MinimalResolution.pdf

The codes can be compiled with GCC. The GNU Multiple Precision Arithmetic Library should be installed.

******************************************************************************************************

To compile, run the following batch files:

sh st_compiling

sh BPtable_complile

sh BP_compile

*******************************************************************************************************

To get the minimal resolution for BP/I, for t<=50, s<=21 (say), run

./mr_st 25 21

To get the structure maps of the BP Hopf algebroid for t<=50, run

./BPtab 25

To get the minimal resolution for BP, for t<=50, s<=20, run

./mr_BP 25 20

*******************************************************************************************************

Warning:

The first input parameter is half of t.

The three executalbes are dependent, and should be run in the above order. 

The s for the minimal resolution for BP/I should be at least one larger than the s for that of BP.

Any mistake of the input could result in unpredictible behaviour, usually a break-down of the program such as a segmentation error.

[EB: the p=3 version may have different restrictions on the degrees, but I haven't thought about what the rule should be.]
