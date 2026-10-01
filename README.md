# MinimalResolution

This is a p=3 fork of [Guozhen Wang's code](https://github.com/pouiyter/MinimalResolution)
for computing the Adams-Novikov E<sub>2</sub> page for the sphere (originally
at p=2), using the algebraic Novikov spectral sequence. See
[this repository](https://github.com/ebelmont/ANSS_data) for sample data and an
explanation of how to interpret the data files output by the program. The
rest of this README is Guozhen's original instructions.

**Range limit.** Coefficients are stored as integers mod 3<sup>40</sup>
(`Z3.cpp`). The computation is incorrect at algebraic Novikov filtration 40
and above. This happens even in degrees where no class of
Ext(BP<sub>*</sub>) is divisible by 3<sup>40</sup>.

The second argument L to `mr_BP` caps algebraic Novikov filtration n, and the
computation is only guaranteed to be correct in the range s + n ≤ L − 1 and
L ≤ 40.

**Consequently the output is only complete for stem + s < 4L, i.e.
stem + s < 160 at L = 40.** The reason: every generator in ANSS filtration s
has internal degree t ≥ 4s, and v<sub>1</sub> is the v<sub>n</sub> of
smallest degree per unit of algebraic Novikov filtration, so a class with no
v<sub>0</sub> in it has s + n ≤ t/4. Classes containing v<sub>0</sub> need
slightly more room, but at degree 185 every difference from an uncapped
computation was at stem + s ≥ 160.

The `Z_mod_3e40` branch removes this limit by computing
Ext(BP<sub>*</sub>/3<sup>40</sup>) instead. Arithmetic mod 3<sup>40</sup> is
exact for that target, and L caps only the ANSS filtration s.

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
