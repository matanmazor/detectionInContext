Research log
================
Matan Mazor
2023-01-04

## Rationale for occluded visual search task (2022-01-04)

In Exp. 3 (occluded detection), subjects did not adjust their
target-absent decision times based on expected search difficulty.
Although they took longer to detect the target when more of the image
was occluded, the proportion of occluded pixels had no effect on the
time it took them to infer absence. This finding is surprising, because
subjects should know that detecting a target is harder when more of the
visual field is occluded, and so they should be more conservative in
detecting absence in such circumstances.

One potential explanation is that people never adjust the timing of
their decisions about absence based on the expected visibility of the
target. For example, Gorea and Sagi (2000) found that subjects fail to
simultaneously hold in working memory two separate detection criteria
for concurrent detection tasks. Interestingly, in visual search tasks,
subjects flexibly adjust their search termination heuristic between
trials as a function of set size and the expected salience of the target
(Mazor and Fleming 2022). In visual search, manipulations that affect
search time in target-present trials only are relatively rare. Together,
it seems that subjects struggle to use metacognitive inference (“I would
have seen it by now”) in a detection setting, but not in a visual search
setting.

One possible explanation is that metacognitive knowledge about attention
and perception is not represented in units of time (e.g., “finding the
item would take me 2 seconds”) but in units of attention allocation
steps instead (e.g., “finding the item would take me 20 saccades”). As a
result, metacognitive stopping can only occur in tasks where subjects
actively allocate attention to different sources.

To test this possibility directly, we can repeat Exp. 3, but this time
the target can appear or not in one of several (e.g., 4) different
locations. Maybe it can look something like this (just a rough sketch):

<img src="figures/occludedSearchSketch.png" width="40%" />

If we now find that the proportion of occluded pixels does have an
effect on search time in target-absent trials, this lends some support
to the idea that people are tracking something other than time in making
decisions about target absence.

## An issue with the experiment code (2022-01-06)

In all versions of the experiment I ran so far, there were two
jatos.onLoad commands, and cryptographic pre-registration only happened
in the second one. I ran some tests and all subjects were run with the
correct version of the onLoad function, because we have 100 unique
subject sums in Exp. 1, 323 in Exp. 2, and 253 in Exp. 3: the same as
the number of participants. I also checked and the setting of the
randomization seed affects the contents of the experiment, even if this
setting only appears in the second function. This bug, which did not
affect the data, will be fixed in all following experiments.

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-gorea2000" class="csl-entry">

Gorea, Andrei, and Dov Sagi. 2000. “Failure to Handle More Than One
Internal Representation in Visual Detection Tasks.” *Proceedings of the
National Academy of Sciences* 97 (22): 12380–84.
<https://doi.org/10.1073/pnas.97.22.12380>.

</div>

<div id="ref-mazor2022" class="csl-entry">

Mazor, Matan, and Stephen M. Fleming. 2022. “Efficient Search
Termination Without Task Experience.” *Journal of Experimental
Psychology: General* 151: 2494–2510.
<https://doi.org/10.1037/xge0001188>.

</div>

</div>