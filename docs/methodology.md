# Measuring Feature Importance with the Survey Q-Sort

**Context, the problem, and the models we should try.**

> Companion methodology note for the estimator recovery study.

---

## The context

Teams constantly have to prioritize a list of things — product features to build, customer needs to address, messages to lead with. The underlying question is simple to state: of these items, which matter most to people, and by how much?

A popular survey tool for answering this is the **Q-sort**. Respondents sort a list of items into a fixed shape — a few "most important," a few "least important," and the bulk in the middle. In its modern survey form it is run as a short sequence of tap-to-choose questions: pick the single most important item, then the next couple, then switch and pick the least important, then the next couple. Whatever is left lands in the middle. Each item then gets a score based on where it ended up — commonly `5, 4, 3, 2, 1` from top band to bottom — and we average those scores across everyone to get an importance ranking.

The method has a long pedigree (it began as a qualitative technique for studying people's viewpoints), but the market-research version keeps only the sorting step. Its credibility rests on a well-known 2006 comparison study that found the Q-sort separates items almost as sharply as the heavier "MaxDiff" method, in less than half the response time, and without specialized software. That combination — *fast, cheap, and decent* — is why it stays in use.

We are doing our own empirical work on this survey version. To do that responsibly, we first have to be precise about one thing that is easy to gloss over: **what quantity does this instrument actually measure well?**

> **What we are trying to estimate**
>
> A population-level importance result for the list — not just the order of the items, but the **spacing** between them: how much more important one item is than another. We want this expressed as **shares that add up to 100%**, where "twice as important" genuinely means twice. Spacing, not just rank, is what drives real decisions about where to invest.

---

## The problem

The standard way of scoring a Q-sort — average the `5-4-3-2-1` codes — quietly assumes two things. First, that the buckets are **evenly spaced**: that the gap between a "5" item and a "4" item is the same as the gap between a "2" and a "1." Second, that **every item belongs** somewhere on the forced shape, even items a respondent doesn't care about at all.

Neither assumption holds in practice. Real importance is lumpy: often one item dominates, a few matter a lot, and there is a long tail nobody cares much about. On top of that, the forced shape caps how far the leader can pull away — a runaway favorite can only be placed in the top bucket, so it cannot show that it is three times more important than the runner-up. The two effects combine to **compress the spacing**: the averaging shrinks the leaders and inflates the tail. The ranking usually still comes out right, but the spacing — the thing we actually care about — comes out wrong.

That has a second, quieter consequence: **uncertainty**. We want a 95% confidence interval that behaves like one — that actually contains the true value 95% of the time. The concern is that the naive average produces tight, confident-looking intervals that are **centered on the slightly-wrong (compressed) number**. When an interval is centered on the wrong value, even a correctly-sized interval misses the truth more often than its label claims. The figure below shows the pattern we expect to find: the choice-based estimate stays honest (near the diagonal), while the naive average sags below it.

![Predicted interval coverage](../figures/coverage_calibration.png)

*Predicted interval coverage. Illustrative pattern, not yet results from a run. A point on the dashed line is honest; below it means the interval is overconfident.*

> **The empirical question, stated precisely**
>
> Across realistic conditions, which estimator recovers the true spacing accurately, and which estimator's confidence intervals are actually honest (cover the truth at their stated rate)? We answer it by simulation: invent a world where the true importances are known, simulate respondents going through the four-pass exercise, then check which method recovers the planted truth and reports trustworthy intervals — repeated many times across different conditions.

---

## What we should try, and why

These are the estimators to put through the simulation. The first two are the status quo we are scrutinizing; the next two are the principled alternative in its two forms; the rest are supporting checks and optional extensions.

One subtlety drives the two Plackett-Luce rows. The exercise is not a single ranking gesture; it is **four separate choices** — two selecting the most important items, two selecting the least. We can either collapse all four into the one ranking they imply and fit a single likelihood, or model the four choices directly, with the two "least" passes treated as choices on the reversed scale (pick the worst). The two are **not the same likelihood** — they assign different probabilities to the very same response and agree only when a choice is between two items — so both are listed, and the simulation decides which fits real respondents. The [appendix](#appendix--why-we-model-the-least-passes-as-separate-trials) works the difference out on a small example.

| Model | What it does / assumes | Why try it, and what we expect |
|---|---|---|
| **Naive mean of bucket codes** (`5-4-3-2-1`) <br> *Baseline* | Averages the band codes across respondents and treats them as evenly-spaced numbers. The current standard practice. | It is the method we are testing. Expected to recover the ranking but compress the spacing, and to produce intervals that under-cover — i.e., it demonstrates the problem. |
| **Same mean, alternative coding vectors** <br> *Baseline stress-test* | Re-runs the naive mean using different numbers for the bands (e.g., a steeper `10-5-3-1-0`, or a flatter set) instead of `5-4-3-2-1`. | Isolates the equal-spacing assumption as the culprit. If the answer shifts just because we changed an arbitrary code, that is direct evidence the spacing claim is unearned. |
| **Plackett-Luce — single partial ranking** <br> *Proposed estimator* | Collapses the four passes into the one partial ranking they imply (top ≻ tied middle ≻ bottom) and fits one Plackett-Luce likelihood, built top-down on the forward worth scale. The convenient default — what the `PlackettLuce` package gives you directly. Runs in R, no Stan. | Recovers spacing far better than the naive mean, and it is the same engine behind MaxDiff. But it scores the "least" passes implicitly, as the tail of a top-down ranking — a specific assumption, not a neutral one (see next row). |
| **Plackett-Luce — four-trial factorization** <br> *Proposed — faithful variant* | Keeps the elicitation's real structure: four conditional choices — two "most" on the forward scale (`v`) and two "least" on the reversed scale (`1/v`) — each over the set that survived to that pass. This is exactly how MaxDiff treats its best/worst picks. | The behaviorally faithful version. It diverges from the single-ranking fit for any choice over three or more items (about 0.58 vs 0.55 on a 3-item least pick), and the gap grows over the real 7-item least pass. Fit both; let the simulation show which recovers the truth with honest coverage. |
| **Respondent-level (cluster) bootstrap** for the naive mean <br> *Required check* | Not a new model — the correct way to compute the naive method's intervals. Resamples whole respondents, never individual item-responses, because each person's scores are locked together (they sum to a fixed total). | Gives the naive estimator its fairest shot at honest intervals and prevents falsely narrow ones. If it still under-covers, that rules out "we just did the intervals wrong" as the explanation. |
| **Dirichlet-multinomial / compositional model** <br> *Optional comparator* | Models each respondent's allocation as a set of shares that sum to one, honoring the locked-together structure without the full choice machinery. | A middle option between the naive mean and Plackett-Luce. Useful to see how much of the fix comes simply from respecting the sum-to-one structure. |
| **Hierarchical Bayes version** (Stan) <br> *Later / only if needed* | Estimates a separate worth profile per respondent, drawn from a population distribution. Out of scope for the first pass. | Needed only if the audience is split into distinct segments with opposite priorities — the case where a single average describes no one. Add it if the first results show that risk is real. |

> **How we will judge them**
>
> - **Accuracy of the spacing:** how close each estimate's shares are to the truth, including the size of the gaps between items — not just whether the order is right.
> - **Honesty of the intervals:** whether a stated 95% interval actually contains the truth about 95% of the time (the figure above), reported alongside how wide the intervals are.
> - **Stratified** for the top/bottom items versus the muddy middle, and repeated across different list lengths, sample sizes, how lumpy the truth is, and how noisy respondents are.

> **Note.** The coverage figure uses illustrative numbers to show the pattern we are testing for; the real curves come out of the simulation. The whole plan runs in R — the hierarchical (Stan) row is optional and only comes into play if heterogeneity turns out to matter.

---

## Appendix — Why we model the "least" passes as separate trials

The survey Q-sort is not one ranking gesture; it is four separate choices — two selecting the most important items, two selecting the least. The two "least" passes are choices on the **reversed scale** (pick the worst), which is not the same as treating the bottom of the list as the tail of a single top-down ranking.

Let each item $i$ have worth $v_i = \exp(\beta_i)$.

**Treat all as one (single Plackett-Luce ranking).** Collapse the four passes into the partial ranking they imply and give it one likelihood, built top-down by repeatedly "choose the best of what remains," always on the forward scale:

$$P(\text{ranking}) = \prod_{t} \frac{v_{(t)}}{\sum_{j \in R_t} v_j}$$

Every factor is a *max*-type choice on $v$; the bottom item earns its probability only by being what is left at the end of a top-down peel. This is what `PlackettLuce` gives you by default.

**Model the trials (best + worst).** Keep the four conditional choices, two forward and two reversed:

$$L = P(\text{most}\mid S_1)\cdot P(\text{2 next-most}\mid S_2)\cdot P(\text{least}\mid S_3)\cdot P(\text{2 next-least}\mid S_4)$$

with the "most" factors forward and max-type, $\dfrac{v_i}{\sum v}$, and the "least" factors reversed and min-type, $\dfrac{1/v_i}{\sum 1/v}$, each over the set that survived to that pass. This is exactly how MaxDiff treats its worst choices.

**Where they split.** Take a single "least" pick over three items $\{B, C, D\}$ with worths $3, 2, 1$, where $D$ is chosen as least. Both models assign a probability to the *identical* event "$D$ is the worst of these three":

$$P(D\text{ last}) = P(B\succ C\succ D) + P(C\succ B\succ D) = \tfrac{3}{6}\cdot\tfrac{2}{3} + \tfrac{2}{6}\cdot\tfrac{3}{4} = 0.333 + 0.250 = 0.583$$

$$P(D\text{ worst}) = \frac{1/v_D}{\,1/v_B + 1/v_C + 1/v_D\,} = \frac{1/1}{\,1/3 + 1/2 + 1/1\,} = \frac{1}{1.833} = 0.545$$

$0.583 \neq 0.545$ — same data, two likelihoods. The single-ranking model only ever peels from the top (forward, $v$), so the bottom of the list is governed by $v$; the trial model peels from both ends, so the bottom is governed by $1/v$. They reconcile only when a choice is between two items (Bradley-Terry symmetry). In the real four-pass over 10 items, the least pass ranges over **seven** survivors, so the gap is larger and compounds across both least passes.

![Single ranking vs four trials](../figures/likelihood_split.png)

*Two ways to score the same Q-sort response. Left: collapse to one Plackett-Luce ranking (peel from the top, forward scale `v`). Right: model the four trials (peel from both ends; "least" uses the reversed scale `1/v`). Same event, two likelihoods — 0.583 vs 0.545.*

This is why the table lists both versions. The single-ranking collapse is the easy default; the four-trial factorization is the faithful one and matches how MaxDiff handles its worst choices. **Which one matches real respondents is an empirical question for the simulation**, not a modeling detail to assume away — so we fit both and compare which recovers the planted truth with honest interval coverage.

### Implementation note

| Model | Engine |
|---|---|
| Single partial ranking | `PlackettLuce::PlackettLuce()` on a tied-ranking object (top-1, tied-2, tied-4, tied-2, bottom-1) |
| Four-trial factorization | `survival::clogit()` on stacked choice tasks; negate item dummies on the two "least" tasks so one coefficient vector governs both directions (`support.BWS::bws.dataset()` shows the encoding pattern, incl. maxdiff / marginal / sequential variants) |
