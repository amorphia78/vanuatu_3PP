# Supplementary Material: Detailed Analysis Methods

Supplement to *The Development of Third-Party Punishment and Reward in ni-Vanuatu Children: Evidence of Cost-Benefit Calculations*

This document describes, in full, every analysis reported in the main manuscript. It is intended to be read alongside the analysis script (`analysis_vanuatu_3pp.R`) and the raw data file (`data_vanuatu_3pp.tsv`), both openly available in the project repository (Anon, 2026), together with a rendered HTML document (`analysis_vanuatu_3pp.html`) showing the complete console output of the script. Section headings below follow the order of the script.

---

## S1. Software

All analyses were conducted in R 4.5.0 (R Core Team, 2025). Model fitting used the `glm` function of the base **stats** package; all classical tests (`binom.test`, `prop.test`, `fisher.test`, `anova.glm`) are also from **stats**. Additional packages were:

- **splines** (bundled with R) for natural cubic spline bases (`ns`; Hastie, 1992).
- **sandwich** 3.1-1 (Zeileis, 2006; Zeileis et al., 2020) for cluster-robust covariance estimation (`vcovCL`).
- **lmtest** 0.9-40 (Zeileis & Hothorn, 2002) for coefficient confidence intervals under a supplied covariance matrix (`coefci`).
- **pwr** 1.3-0 (Champely, 2020) for the sensitivity (minimum detectable effect) analysis.

Figures were drawn with base R graphics. The analysis involves no random number generation, no simulation and no stochastic model fitting, so results are exactly reproducible without setting a seed.

---

## S2. Data file and derived variables

The data file contains one row per participant and no missing values (303 complete cases, 303 unique participant identifiers). The columns are:

| Column | Description |
| --- | --- |
| `Age` | Age in decimal years at test |
| `ID` | Participant identifier |
| `Sex` | `F` / `M` (referred to as gender in the manuscript) |
| `Location` | `Urban` (Port Vila, Efate) / `Rural` (*Niprainitata*, Tanna) |
| `Condition` | `Private` = anonymous / `Public` = in person |
| `Cost` | `N` = free of economic cost / `Y` = economically costly |
| `AntiGot` | Allocation to the antisocial actor: `Bad` / `Empty` / `Good` |
| `NeutGot` | Allocation to the neutral actor: `Bad` / `Empty` / `Good` |

Note that the `Condition` labels are historical: `Public` denotes the **in-person** condition (in which the child expected to hand the box over face to face) and `Private` denotes the **anonymous** condition (in which the box was to be posted).

From `AntiGot` and `NeutGot`, six binary outcome variables were derived, one for each combination of actor (antisocial, neutral) and allocation (bad sweets = punishment, empty box, good sweets = reward). Each is coded `Y` if that item was allocated to that actor and `N` otherwise. Because factor levels sort as `N` then `Y`, all logistic models predict the probability of the `Y` (item-allocated) outcome. The three allocations to a given actor are mutually exclusive and exhaustive, so for each actor the three binary indicators sum to the full sample (antisocial: 149 bad, 49 empty, 105 good; neutral: 35 bad, 43 empty, 225 good). Each binary outcome therefore contrasts one allocation against the other two pooled, which is the quantity of substantive interest ("how often was this actor punished / rewarded?") rather than a contrast among the three options (see S14).

A centred age variable, `AgeCentred`, was created by subtracting the sample grand mean (7.00 years) from `Age`. Centring does not change the coefficients for the experimental factors in models without interactions, but it makes the intercept and the main effects in the interaction models interpretable at the mean age rather than at age zero. The models used for Figure 2 and for the cross-cultural comparison use uncentred `Age`, because for those models predictions are required at a specified age on the raw scale.

---

## S3. Sample and design summaries

Participant descriptives reported in the Participants section — sample size, mean, standard deviation and range of age, and gender counts, for the whole sample and separately for the urban and rural samples — are printed by the script directly from the data file.

Table 2 of the manuscript gives the number of urban children by year of age within each of the four between-subjects cells (economic cost × anonymity). Age groups were formed by truncating decimal age to the integer year (`floor`), so that, for example, the "5" row contains all children aged 5.0 to 5.99 years. The rural column of Table 2 was formed with the same year bands. All 21 rural children were tested in a single cell: **anonymous** and **free of economic cost**.

---

## S4. Preferential punishment of the antisocial actor (RQ1, RQ2)

Each child allocated an item to *both* actors, so the two punishment proportions are paired and cannot be compared with a test for independent samples. Two quantities are reported.

**The p value** is McNemar's exact test: an exact two-sided binomial test, against a null probability of .5, on the participants who punished exactly one of the two actors (the discordant pairs). Over the whole sample, 139 children punished only the antisocial actor and 25 only the neutral actor. This test evaluates the null hypothesis of marginal homogeneity — that punishment is equally likely for the two actors — and is exact, so it remains valid in the small age subgroups.

**The odds ratio and its confidence interval** compare the two marginal proportions rather than the discordant pairs. They were obtained by restructuring the data into long form (two rows per child, one per actor, with an indicator for actor type), fitting an ordinary logistic regression of punishment on actor type, and estimating the sampling variance with a cluster-robust ("sandwich") covariance matrix clustered on participant (`sandwich::vcovCL`), from which the interval was computed with `lmtest::coefci`. Clustering accounts for the repeated measurement of each child. The resulting odds ratio is the ratio of the *marginal* odds of punishment for the two actors; it is a different, though related, quantity from the discordant-pair (conditional) odds ratio implicit in McNemar's test, and both are equal to 1 under the same null hypothesis.

The same procedure was applied to the whole sample and to three age subgroups, defined as 4.0 ≤ age < 5.0 (*n* = 32), 5.0 ≤ age < 6.0 (*n* = 65) and age ≥ 6.0 (*n* = 206). The subgroups were chosen a priori to address RQ2, which asked at what age preferential punishment is first evident; the oldest group was left unsplit because splitting it further would have produced cells too small to be informative. No adjustment was made for the three subgroup tests.

---

## S5. Logistic regression models of allocation (Tables 4 and 5; RQ3–RQ9)

For each of the six binary outcomes, two logistic regression models were fitted to the full sample of 303 children (binomial family, logit link, `glm`):

- a **main-effects model** with anonymity, economic cost, gender and centred age as predictors; and
- an **interaction model** adding anonymity × age and gender × age.

The interaction with economic cost was not fitted. The economically costly condition was introduced part-way through data collection and, as Table 2 shows, contains 45 children who are almost all under 7 years old, which leaves too little overlap with the free condition across the age range to estimate an age-varying cost effect (RQ9).

Reference categories are anonymous (for anonymity), free of economic cost (for cost) and girls (for gender), so the tabulated odds ratios describe the in-person condition relative to anonymous, the costly condition relative to free, and boys relative to girls. Because age is centred, the odds ratio for age is the multiplicative change in the odds per additional year, and in the interaction models the main effects of the experimental factors and gender are their effects at the sample mean age of 7.00 years.

Reporting in Tables 4 and 5 follows a two-stage convention that is worth stating explicitly:

- The **main-effect rows** (anonymity, economic cost, gender, age) are taken from the main-effects model, so they estimate average effects across the age range and are not conditioned on the interaction terms.
- The **interaction rows** are taken from the interaction model.
- Between them, the tables report a single **likelihood-ratio test** of the interaction block, obtained by an analysis of deviance comparing the two nested models on 2 degrees of freedom (`anova(..., test = "Chisq")`). This is the omnibus test of whether age moderates the effects of anonymity or gender, and is the appropriate gate on interpreting the individual interaction terms below it.

Odds ratios are exponentiated coefficients. The confidence intervals are **profile-likelihood** intervals (`stats::confint` applied to a `glm` object), which have better small-sample coverage than Wald intervals for logistic models; the *z* statistics and *p* values in the same rows are Wald statistics from the model summary. For these data the two interval types agree to within about 0.01 on the odds-ratio scale, so no inferential inconsistency arises, but the distinction is noted for completeness.

The models for allocation of the empty box were fitted (they are used for Figure 2 and for the mutually exclusive breakdown shown there) but are not tabulated, since neutral outcomes were not a research question.

---

## S6. Figure 2: model-based allocation probabilities

Figure 2 shows a fitted probability for each of the six allocations, in a panel per actor (columns) and per allocation (rows). Because the design is unbalanced (Table 2), raw cell proportions confound the experimental factors with age and with each other, so the figure shows model predictions instead.

A separate set of six logistic models was fitted for this figure, each with anonymity, economic cost and **uncentred** age as predictors. Gender is not included, since the figure does not display it and it has no material effect on the other estimates. Predictions were made at the **median age of the sample, 6.6 years**, which is printed to the console by the script.

Within each panel the three bars are **design cells rather than averages over a factor**: a common baseline in which neither manipulation applies (anonymous and free of economic cost), followed by each manipulation applied on its own (in person while still free of economic cost; economically costly while still anonymous). Every bar is therefore a prediction for a cell in which children were actually tested, and can be read against the corresponding row of Table 2. Three consequences of this arrangement are worth noting:

- Each manipulation is read against the *same* reference bar within the panel, so the two experimental effects are directly comparable and the baseline is drawn once rather than repeated.
- Because the models contain no interaction between the two experimental factors, the odds ratio between the baseline bar and either manipulation bar is exactly the odds ratio for that factor in the underlying model. Averaging over the other factor on the probability scale would not preserve this, since the odds ratio is not collapsible.
- The reference cell is the same one used for Figure 3 (free of economic cost), Figure 5 (free of economic cost) and Figure 6 (anonymous and free of economic cost), so all four figures describe children in a consistent condition.

The fourth design cell, in which both manipulations apply, is not drawn in Figure 2, since the joint effect of the two manipulations is not a research question here. Those 22 children nevertheless inform every bar in the figure: the models are fitted to all 303 participants, so the cell contributes to each coefficient, and because the models contain no interaction between the two factors there is nothing in it that the three plotted cells fail to capture. The cell is shown, with its observed proportion, in Supplementary Figure 1.

The alternative of averaging over the other factor, weighted by how participants happened to divide between its levels (a predictive margin; Graubard & Korn, 1999), was considered and not adopted. Those weights would carry no substantive meaning here, because the split between the cost conditions (258 free, 45 costly) is an artefact of the economically costly condition having been introduced part-way through data collection rather than a feature of any population of interest. The cost of the choice made is a modest loss of precision: holding the other factor fixed uses a narrower slice of the design than averaging over it, and widens the plotted intervals by roughly a tenth.

Intervals were formed on the link (log-odds) scale from the prediction standard error and back-transformed through the inverse logit, so that they cannot extend outside (0, 1). This is the same interval construction used for Figures 3, 4 and 5. The plotted values are printed to the console by the script.

The figure carries no significance markers; formal inference for these predictors is in Tables 4 and 5.

---

## S7. Supplementary Figure 1: model fit against raw data

Supplementary Figure 1 is the same figure with observed proportions plotted alongside each fitted bar, to allow the reader to check that the model-based summary does not distort the data. Each raw bar is the observed proportion in the same design cell that its fitted neighbour predicts, with a 95% interval from `prop.test` (Wilson score interval with continuity correction). The only difference within a pair of bars is therefore that the fit holds age at the sample median whereas the raw proportion pools across the age range.

It also adds the fourth design cell, in which children were both in person and in the economically costly condition, so that every cell of the design appears somewhere and no participant is absent from the figure. Across the resulting 24 pairs, fitted and observed values differ by 0.036 on average and by at most 0.11.

Because the raw bars are restricted to a single design cell, some rest on small samples — in particular the two economically costly cells, drawn from 23 and 22 children respectively. Their intervals are correspondingly wide, and the comparison with the fitted bars is weak in those cells: the four largest fit-to-data discrepancies all fall there. Both figures are produced by the same drawing code, differing only in whether the fourth cell and the raw bars are drawn.

<!-- INSERT supplementaryFigure.png HERE -->

***Supplementary Figure 1.*** *Probabilities of allocations to the antisocial actor (left panels) and the neutral actor (right panels), with the corresponding observed proportions. Dark bars are model predictions: fitted values (±95% CI) for a child of median age (6.6 years) from six logistic regression models of each type of allocation to each type of actor, with predictors age, economic cost and anonymity. Bars run as a baseline in which neither manipulation applies (anonymous and free of economic cost), then each manipulation on its own, then both together; the first three are those shown in Figure 2 of the main manuscript. Every bar is a prediction for a design cell in which children were tested. Pale bars are the observed proportions (±95% CI) in that same design cell, pooled across age. Fitted and observed values differ by 0.036 on average across the 24 pairs. Cell sizes for the observed proportions are 152 (anonymous, free of economic cost), 106 (in person, free of economic cost), 23 (anonymous, economically costly) and 22 (in person, economically costly).*

---

## S8. Age-flexible (spline) models: Figures 3 and 4a

The models in Tables 4 and 5 constrain the effect of age on the log-odds scale to be linear. To examine whether the anonymity effect varied with age in a more complex way, exploratory models were fitted in which age enters through a **natural cubic spline** basis (`splines::ns`; Durrleman & Simon, 1989; Hastie, 1992). Two interior knots were placed at the 33rd and 67th percentiles of age (6.12 and 7.67 years), dividing the sample into three approximately equal-sized age bands, with boundary knots at the observed age range (4.15 and 10.91 years). A natural spline constrains the fit to be linear beyond the boundary knots, which limits the instability that otherwise arises at the sparse extremes of the age distribution. Knot placement by quantile, rather than at substantively chosen ages, was decided before inspecting the fits.

- **Figure 3** (punishment of the antisocial actor) comes from a model with anonymity, economic cost, the age spline, and the anonymity × age-spline interaction. Fitted lines are drawn for the free-of-economic-cost condition, and the plotted points are observed proportions from that condition only, so lines and points describe the same subsample.
- **Figure 4a** (reward of the antisocial actor) comes from a model with anonymity, the age spline, and their interaction. This model contains no economic-cost term, and the plotted points are observed proportions from the whole sample, so lines and points here describe all children rather than one cost condition.

For both figures, confidence bands were obtained by predicting on the link (log-odds) scale with standard errors (`predict(..., type = "link", se.fit = TRUE)`), forming a symmetric 95% interval there, and back-transforming through the inverse logit. Constructing the interval on the link scale keeps it inside (0, 1).

The plotted points are observed proportions in one-year age bands centred on 4.5, 5.5, …, 9.5 years, each band covering the year up to and including its upper bound. The eight children older than 10.0 years fall outside the plotted bands; they are included in all model fits.

These spline analyses are exploratory and are used to characterise the *shape* of the age relationship; the confirmatory tests of the anonymity × age interaction are the likelihood-ratio tests reported in Tables 4 and 5.

### S8.1 Age-specific tests of the anonymity contrast

Because the spline fits in Figures 3 and 4a suggest that the two anonymity conditions separate during middle childhood and converge at either end of the age range, that impression was examined formally. Two different questions have to be kept apart, and the script reports both.

Supplementary Tables 1 and 2 report, for each year of age, the difference between the two fitted anonymity conditions tested against zero. These contrasts are Wald tests of a linear combination of the model coefficients, with standard errors taken from the model covariance matrix, and they follow the pattern the figures suggest: for punishment the effect is clear at ages 6 and 7 and not detectable at 4, 5, 8, 9 or 10; for reward it is clear at ages 5, 6 and 7 and not detectable outside that range.

Two points bear on how these contrasts should be read. First, their significance reflects precision as well as effect size, and precision is poorest at the ends of the age range where the cells are smallest — the punishment contrast at age 4 has a confidence interval running from 0.07 to 17.79, so the youngest children are better described as uninformative than as showing no effect. Second, the omnibus likelihood-ratio test of the whole anonymity × spline-age block is not significant for either behaviour (punishment, Δdeviance = 5.78, 3 *df*, *p* = .123; reward, Δdeviance = 5.83, 3 *df*, *p* = .120), although these tests spread their power across three degrees of freedom and are correspondingly conservative. For reward, the linear anonymity × age interaction reported in Table 5 is significant (*p* = .042), so an age-varying effect is supported there even if its precise shape is not.

The position taken in the manuscript is therefore as follows. For both behaviours the anonymity effect is clear during middle childhood — at ages 6 and 7 for punishment, and at ages 5 to 7 for reward — and is undetectable in the oldest children. The fitted curves are consistent with the two conditions following a similar developmental path offset in time, the in-person group reaching comparable rates a year or two later. This characterisation of the age dependence is exploratory and is offered as a description of the fitted curves rather than as a demonstrated interaction.

***Supplementary Table 1.*** *Contrasts between anonymity conditions for punishment of the antisocial actor, from the spline model shown in Figure 3, in the free-of-economic-cost condition. Odds ratios above 1 indicate more punishment when in person than when anonymous.*

| Contrast | OR [95% CI] | Wald *z* | *p* |
| --- | --- | --- | --- |
| Age 4 | 1.07 [0.07, 17.79] | 0.05 | .962 |
| Age 5 | 0.36 [0.12, 1.12] | −1.77 | .077 |
| Age 6 | **0.22 [0.10, 0.52]** | **−3.45** | **.001** |
| Age 7 | **0.42 [0.20, 0.88]** | **−2.31** | **.021** |
| Age 8 | 0.95 [0.36, 2.48] | −0.11 | .910 |
| Age 9 | 1.12 [0.48, 2.61] | 0.25 | .801 |
| Age 10 | 0.82 [0.16, 4.10] | −0.25 | .804 |

***Supplementary Table 2.*** *Contrasts between anonymity conditions for reward of the antisocial actor, from the spline model shown in Figure 4a. Odds ratios above 1 indicate more reward when in person than when anonymous.*

| Contrast | OR [95% CI] | Wald *z* | *p* |
| --- | --- | --- | --- |
| Age 4 | 1.08 [0.13, 9.02] | 0.08 | .940 |
| Age 5 | **2.72 [1.14, 6.50]** | **2.25** | **.024** |
| Age 6 | **4.01 [1.76, 9.13]** | **3.31** | **.001** |
| Age 7 | **2.28 [1.07, 4.86]** | **2.13** | **.033** |
| Age 8 | 1.04 [0.36, 2.97] | 0.07 | .942 |
| Age 9 | 0.74 [0.27, 1.99] | −0.60 | .546 |
| Age 10 | 0.73 [0.09, 5.92] | −0.30 | .766 |

**Note.** Statistical significance is marked by bold type. All values are printed by the analysis script.

---

## S9. Gender × age model: Figure 4b

Figure 4b shows reward of the antisocial actor by gender and age, from a logistic model with gender, centred age and their interaction. Age is entered **linearly** here (no spline), and the model contains no anonymity or economic-cost terms; it is the two-predictor model that the figure displays. Points are observed proportions in the same one-year age bands as in S8, computed separately for girls and boys across the whole sample. Confidence bands were constructed as in S8.

---

## S10. Rewarding the antisocial actor and punishing the neutral actor

The association reported in the Results was first described as a simple contrast: the proportion of children rewarding the antisocial actor among those who punished the neutral actor (*n* = 35), against the proportion among those who did not (*n* = 268).

It was then tested by adding punishment of the neutral actor as a predictor to the interaction model for reward of the antisocial actor described in S5 (anonymity, economic cost, gender, centred age, anonymity × age, gender × age). Model improvement was assessed by a likelihood-ratio test on 1 degree of freedom against the model without that term. The odds ratio for the added term is the exponentiated coefficient, with a profile-likelihood confidence interval and a Wald *z* and *p*.

This analysis is post hoc, prompted by the unexpected frequency of reward of the antisocial actor, and is reported as such. Because punishment of the neutral actor is a behaviour measured at the same moment as the outcome rather than a prior or manipulated variable, the coefficient describes an association between two simultaneous choices and carries no causal interpretation in either direction.

---

## S11. Cross-cultural comparison with Kenward and Östh (2015): Figure 5

The Swedish data are not available as a raw file; the counts were read from Figure 2 of Kenward and Östh (2015) and entered into the script as literal values. In that study 24 children were tested per anonymity condition. The counts used are: anonymous condition, 14 of 24 punished the antisocial actor and 1 of 24 the neutral actor; in-person condition, 5 of 24 and 1 of 24 respectively. Swedish proportions and their 95% intervals were computed with `prop.test`.

Three features of the Swedish sample constrain the comparison. It was tested only under the anonymity manipulation, with no economic cost manipulation; it covers a narrow age range with a mean of 5.2 years; and only summary data are available. The comparison was therefore restricted to punishment, restricted to ni-Vanuatu children in the free-of-economic-cost condition, and made at the Swedish mean age.

The Vanuatu values are model predictions rather than raw proportions, because no ni-Vanuatu age band matches the Swedish sample closely enough to yield a usable raw estimate. Two models were fitted, for punishment of the antisocial and of the neutral actor, each with anonymity, economic cost, uncentred age, and anonymity × age. Predictions were made at age 5.167 years (5 years 2 months) in the free-of-economic-cost condition, with intervals formed on the link scale and back-transformed as in S8.

No *p* values are reported. The comparison is between a model-based prediction for one sample and a raw proportion from another, the two samples differ on multiple dimensions besides culture, and the analysis was framed as exploratory from the outset; the manuscript therefore restricts itself to describing the overlap of the two sets of intervals. The plotted values are printed to the console by the script.

---

## S12. Urban–rural comparison: age matching and Figure 6

The rural sample (*n* = 21) is too small for multivariate modelling and is not age-matched to the urban sample, so the comparison was made against a **matched urban subsample** rather than by regression (following the approach of Robbins & Rochat, 2011, Study 6).

The donor pool was restricted to urban children in the same experimental cell as the rural children — anonymous and free of economic cost (*n* = 131, of whom 94 fall within the rural age range). Matching was greedy, deterministic and without replacement: taking the rural children in data order, each was matched to the not-yet-selected urban child with the smallest absolute age difference, with ties broken by data order. Because no randomisation is involved, the procedure returns the same subsample on every run. Four complete passes were made, so each rural child contributed four urban matches, giving an urban comparison subsample of 84 children.

Four passes were chosen as the largest number that preserved close age matching. The script prints, after each pass, the running mean absolute age difference over all matches made so far: 0.05 years after the first pass, rising to 0.22 years (about two and a half months) after the fourth. Matching therefore remains close even at the fourth pass, when each rural child is being paired with its fourth-nearest available urban counterpart.

The matching succeeded to one decimal place on all three summary statistics: rural mean 7.3 years, range 5.5 to 9.9; matched urban mean 7.3 years, range 5.5 to 9.9.

Figure 6 shows observed proportions for the four behaviours of interest in the two samples, with 95% intervals from `prop.test`. Fisher's exact tests comparing the two samples on each behaviour were also computed and are printed by the script; all were far from significance. Because the rural sample is small, this comparison has limited power, and the manuscript's Limitations section states that the intervals remain compatible with moderate undetected differences.

---

## S13. Sensitivity (minimum detectable effect) analysis

The study was not powered in advance: sample size was determined by the time available across three field visits. A **sensitivity analysis** was therefore conducted post hoc, reporting the smallest effect each predictor could detect with 80% power at α = .05 (two-sided) given the realised sample sizes. This is reported in place of an observed (post hoc) power calculation, which would be uninformative because it is a deterministic function of the obtained *p* value.

For the three binary predictors, `pwr::pwr.2p2n.test` was used with the actual group sizes (anonymity 175 / 128; economic cost 45 / 258; gender 148 / 155). The function returns Cohen's *h*, which was converted to a second proportion and thence to an odds ratio under an assumed baseline probability of .5.

For age, `pwr::pwr.r.test` was used with *n* = 303 to obtain the minimum detectable correlation (*r* = .160), read here as the point-biserial correlation between age and the binary behaviour. Because the logistic approximation that yields an odds ratio takes Cohen's *d* rather than *r*, the correlation was first converted to *d* using the point-biserial relation at the .5 outcome split assumed throughout this section, giving *d* = 0.324. That was converted to an odds ratio per standard deviation of age (1.80) and rescaled to an odds ratio per year (1.43) using the observed age standard deviation of 1.64 years.

These figures are approximations and should be read as such. Each treats its predictor in isolation, without the other covariates present in the fitted models, and the assumed baseline of .5 is a convention rather than the observed base rate for any particular behaviour — it is the value at which the conversion from Cohen's *h* to an odds ratio is exact. They are intended to convey the approximate resolution of the design, not to substitute for the confidence intervals in Tables 4 and 5, which are the appropriate basis for judging the precision of any individual estimate. It is worth noting in particular that the detected age effect on punishment of the antisocial actor (OR 1.41 per year) sits close to this threshold, so the design had approximately the resolution needed to detect an effect of the size observed.

---

## S14. Analytic decisions and their rationale

**Separate binary models rather than a multinomial model.** The three allocations available for each actor are mutually exclusive, so an alternative would be a multinomial model with the empty box as reference. Six separate binary models were preferred because the research questions concern the rate of punishment and the rate of reward — "how often was this actor punished?" — rather than the relative preference among three options, and because binary odds ratios are directly comparable with the effect sizes reported in the developmental literature this study engages with. A multinomial specification gives the same qualitative conclusions for every predictor. One consequence of the binary specification is worth flagging: the economic-cost manipulation acts largely by moving children toward the empty box (Figure 2, middle row), so its odds ratios for punishment and reward reflect that shift.

**Separate models for the two actors.** Allocations to the antisocial and to the neutral actor were modelled separately rather than jointly, because the manipulations were expected (and found) to act differently on the two, and because the substantive claims concern each behaviour in its own right. The one comparison that is directly between the two actors — RQ1 and RQ2 — is the paired analysis in S4, which does account for the within-participant structure.

**No adjustment for multiple comparisons.** The punishment analyses address the a priori research questions listed in Table 1 and were specified before analysis, though not preregistered (data collection began in 2013). The reward analyses, and the association between reward of the antisocial actor and punishment of the neutral actor, are exploratory and are labelled as such throughout the manuscript. No family-wise correction was applied; readers should weight the exploratory results accordingly.

**Retention of the rural children in the main models.** All 303 children are included in the models of Tables 4 and 5, and location is not entered as a covariate. The 21 rural children were all tested anonymously and free of economic cost, so they contribute to the intercept and to the age term but not to the contrast between cost conditions.

**Age range of the economic cost manipulation.** As Table 2 shows, the economically costly condition was administered almost entirely to children under 7 (mean age 5.7 years, against 7.2 in the free condition), because the condition was added late in data collection. The estimated cost effect therefore relies on the age term to adjust for this imbalance, and the two conditions have limited common support above age 7. This is the reason no cost × age interaction was fitted (RQ9), and it is a limitation of the cost estimates rather than of the anonymity estimates: assignment to the anonymity conditions was random and the two anonymity groups are well balanced on age (means 7.0 and 6.9 years) and gender.

---

## S15. Reproducibility

All analyses are contained in the single script `analysis_vanuatu_3pp.R`, which reads the single data file and runs top to bottom. The first line calls a wrapper function that sets the working directory; users should replace it with a `setwd()` call for their own copy of the repository. The script writes the following outputs:

| File | Corresponds to |
| --- | --- |
| `figure2.png` | Figure 2 |
| `supplementaryFigure.png` | Supplementary Figure 1 |
| `figure3.png` | Figure 3 |
| `figure4a.png`, `figure4b.png` | Figure 4a, 4b |
| `OR_table_punishment.csv` | Table 4 |
| `OR_table_reward.csv` | Table 5 |

Figures 5 and 6 were drawn from the values printed to the console by the cross-cultural and urban–rural sections of the script; those printed tables contain every plotted proportion and interval. Figure 1 is a photograph of the apparatus and involves no analysis.

A rendered HTML document (`analysis_vanuatu_3pp.html`), produced from the script with `rmarkdown::render`, reproduces the complete console output including all model summaries, so that every number in the manuscript can be located without re-running the code. The full experimental script is available as `procedure_details.pdf`.

---

## References

Champely, S. (2020). *pwr: Basic functions for power analysis* (R package version 1.3-0). https://CRAN.R-project.org/package=pwr

Durrleman, S., & Simon, R. (1989). Flexible regression models with cubic splines. *Statistics in Medicine, 8*(5), 551–561. https://doi.org/10.1002/sim.4780080504

Graubard, B. I., & Korn, E. L. (1999). Predictive margins with survey data. *Biometrics, 55*(2), 652–659. https://doi.org/10.1111/j.0006-341x.1999.00652.x

Hastie, T. J. (1992). Generalized additive models. In J. M. Chambers & T. J. Hastie (Eds.), *Statistical models in S* (pp. 249–307). Wadsworth & Brooks/Cole.

Kenward, B., & Östh, T. (2015). Five-year-olds punish antisocial adults. *Aggressive Behavior, 41*(5), 413–420. https://doi.org/10.1002/ab.21568

R Core Team. (2025). *R: A language and environment for statistical computing* (Version 4.5.0). R Foundation for Statistical Computing. https://www.R-project.org/

Robbins, E., & Rochat, P. (2011). Emerging signs of strong reciprocity in human ontogeny. *Frontiers in Psychology, 2*, 353. https://doi.org/10.3389/fpsyg.2011.00353

Zeileis, A. (2006). Object-oriented computation of sandwich estimators. *Journal of Statistical Software, 16*(9), 1–16. https://doi.org/10.18637/jss.v016.i09

Zeileis, A., & Hothorn, T. (2002). Diagnostic checking in regression relationships. *R News, 2*(3), 7–10. https://CRAN.R-project.org/doc/Rnews/

Zeileis, A., Köll, S., & Graham, N. (2020). Various versatile variances: An object-oriented implementation of clustered covariances in R. *Journal of Statistical Software, 95*(1), 1–36. https://doi.org/10.18637/jss.v095.i01
