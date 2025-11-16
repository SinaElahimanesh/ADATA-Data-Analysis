# ADATA Post-Task Questionnaire: Comprehensive Baseline vs Experimental Analysis

## Dataset overview
- Sources: merged three “OCJS Resubmission – ADATA Post-Task Questionnaire” CSV exports.
- Responses after filtering (excludes P13, P14, p15, P16): 53 total rows, 24 with explicit condition labels (12 Baseline, 12 Experimental System).
- Measures: time estimate (minutes), satisfaction, 5-point agreement items, NASA-TLX style workload questions, feature helpfulness ratings, rich open-ended comments.
- Comments: 22 text responses were detected across open-ended prompts; each is now labelled with themes & sentiment for downstream qualitative analysis (`outputs/qualitative/comment_labels.csv`).
- Outputs: see `outputs/` (tables), `outputs/stats/` (tests), `outputs/qualitative/` (themes, TF-IDF, sentiment, representative quotes).
- Figures: see `figs/` (`ratings_comparative.svg`, `time_satisfaction_comparative.svg`, `task_completion_comparative.svg`, condition-specific word clouds).

## Methods summary
- Ingest & normalization: all matching CSVs merged via union schema; excluded specified participants; whitespace-normalized categorical values; Likert text converted to numeric via leading integer.
- Quantitative:
  - Descriptive stats (mean/median/std/min/max) overall + per condition for every metric (`outputs/numeric_summary.csv`).
  - Difference table with effect sizes (`outputs/numeric_differences.csv`).
  - NASA-TLX composite (average of Mental, Temporal, Frustration) derived per row.
  - Permutation tests (10k iterations) for numeric and categorical outcomes, with bootstrap CIs and Holm-adjusted p-values (`outputs/stats/*.csv`).
- Qualitative:
  - Comment detection across all open-text fields; tokenization + TF-IDF (`outputs/qualitative/condition_tfidf.csv`).
  - Thematic labeling via keyword dictionaries (`outputs/qualitative/comment_labels.csv` & `.../theme_prevalence_by_condition.csv`).
  - Lexicon-based sentiment (positive/neutral/negative) per condition (`outputs/qualitative/sentiment_by_condition.csv`).
  - Representative quote extraction and SVG word cloud summaries.

---

## 1. Quantitative findings

### 1.1 Condition participation
| Metric | Baseline Count | Baseline % | Experimental System Count | Experimental System % |
|---|---:|---:|---:|---:|
| Responses with labeled condition | 12 | 50.0% | 12 | 50.0% |

**Implication:** Balanced samples enable fair Baseline vs Experimental comparisons.

### 1.2 Task completion
| Completion | Baseline Count | Baseline % | Experimental System Count | Experimental System % |
|---|---:|---:|---:|---:|
| Partially | 1 | 8.3% | 0 | 0.0% |
| Yes, Completely | 10 | 83.3% | 11 | 91.7% |
| Yes, Mostly (minor elements incomplete) | 1 | 8.3% | 1 | 8.3% |

**Implication:** Completion outcomes are comparable; differentiation must come from workload, ideation ease, and feature utility.

### 1.3 Time, satisfaction, and NASA-TLX composite
- Overall average time = **17.6 minutes** (Baseline 24.3 vs Experimental 10.9).
- Overall satisfaction = **3.92/5** (Baseline 3.58 vs Experimental 4.25).
- NASA-TLX composite (mean of Mental, Temporal, Frustration):
  - Baseline ≈ **3.36**
  - Experimental ≈ **1.69**
- See `outputs/numeric_summary.csv` for means/medians/std per metric and `outputs/numeric_differences.csv` for effect sizes.

**Implication:** Experimental users work significantly faster, report higher satisfaction, and roughly halve subjective workload.

### 1.4 Ratings (overall descriptive)
| Question | Average | N |
|---|---:|---:|
| [Frustration Level…] | 2.67 | 24 |
| [Mental Demand…] | 2.54 | 24 |
| [Temporal Demand…] | 2.38 | 24 |
| Agreement [Confident aligns with topic] | 3.75 | 24 |
| Agreement [Innovation/creativity] | 3.21 | 24 |
| Agreement [Pedagogical quality is low] | 2.71 | 24 |
| Agreement [It was easy to come up with an idea] | 2.83 | 24 |
| Agreement [It was hard to write/format] | 2.58 | 24 |

### 1.5 Ratings (Baseline vs Experimental)
| Question | Baseline Avg | Baseline N | Experimental Avg | Experimental N |
|---|---:|---:|---:|---:|
| [Frustration Level…] | 3.00 | 12 | 2.33 | 12 |
| [Mental Demand…] | 3.67 | 12 | 1.42 | 12 |
| [Temporal Demand…] | 3.42 | 12 | 1.33 | 12 |
| Agreement [Confident aligns with topic] | 3.50 | 12 | 4.00 | 12 |
| Agreement [Innovation/creativity] | 2.83 | 12 | 3.58 | 12 |
| Agreement [Pedagogical quality is low] | 2.83 | 12 | 2.58 | 12 |
| Agreement [It was easy to come up with an idea] | 1.92 | 12 | 3.75 | 12 |
| Agreement [It was hard to write/format] | 3.17 | 12 | 2.00 | 12 |

**Implication:** Experimental users experience dramatically lower workload and find ideation far easier; formatting stress drops by ≥1 point.

### 1.6 Feature helpfulness (Experimental respondents)
| Feature | % rating 4–5 |
|---|---:|
| Brainstorm module helpfulness | 90.9% |
| Design module helpfulness | 83.3% |
| Generate module helpfulness | 75.0% |
| Document generation time-saving | 100% |

**Implication:** The integrated workflow is consistently perceived as helpful, with PDF generation unanimously saving time.

### 1.7 Permutation tests (10k iterations, Holm-adjusted p)
| Measure | Baseline Mean | (CI) | Experimental Mean | (CI) | Mean Diff (B–E) | Cohen d | Cliff Δ | p | Holm p | N B | N E |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Time (minutes) | 24.33 | [17.67, 31.42] | 10.92 | [9.33, 12.92] | 13.42 | 1.45 | 0.68 | 0.0024 | 0.0192 | 12 | 12 |
| Satisfaction (1–5) | 3.58 | [3.08, 4.08] | 4.25 | [3.92, 4.58] | -0.67 | -0.86 | -0.44 | 0.0844 | 0.4482 | 12 | 12 |
| Easy to come up with idea | 1.92 | [1.42, 2.50] | 3.75 | [2.92, 4.50] | -1.83 | -1.45 | -0.67 | 0.0047 | 0.0329 | 12 | 12 |
| Mental Demand | 3.67 | [3.00, 4.25] | 1.42 | [1.17, 1.67] | 2.25 | 2.38 | 0.85 | 0.0002 | 0.0020 | 12 | 12 |
| Temporal Demand | 3.42 | [2.58, 4.17] | 1.33 | [1.08, 1.58] | 2.08 | 1.86 | 0.75 | 0.0005 | 0.0045 | 12 | 12 |
| Other measures | (see `outputs/stats/numeric_tests.csv` for full table + CI/adjusted p) |

**Significant improvements (Holm p < 0.05):**
- Time to complete task
- Ease of coming up with an idea
- Mental demand
- Temporal demand

**Implication:** The Experimental workflow tangibly reduces cognitive/temporal load and accelerates ideation; remaining measures trend positive but need larger samples for definitive claims.

Categorical tests (task completion, module helpfulness) show no significant differences because Baseline respondents rarely rated those modules (most ratings come from Experimental users).

### 1.8 Figures
- Ratings comparison: `figs/ratings_comparative.svg`
- Time and satisfaction comparison: `figs/time_satisfaction_comparative.svg`
- Task completion comparison: `figs/task_completion_comparative.svg`

---

## 2. Qualitative findings
Data-driven outputs: `outputs/qualitative/` (comment labels, theme prevalence, TF-IDF, sentiment, representative quotes); SVG word clouds in `figs/`. Every detected comment is auto-labeled with themes and sentiment (`outputs/qualitative/comment_labels.csv`) for traceability.

### 2.1 Theme prevalence (mentions per comment)
| Condition | Brainstorm support | Design refinement | LaTeX/PDF | Formatting effort | Helpfulness of AI | Time saving |
|---|---:|---:|---:|---:|---:|---:|
| Baseline (10 comments) | 0.40 | 0.30 | 0.20 | 0.40 | 0.40 | 0.10 |
| Experimental (12 comments) | 0.58 | 0.50 | 0.58 | 0.33 | 0.67 | 0.25 |

**Implication:** Experimental comments focus on AI assistance, document generation, and idea support; Baseline comments emphasize manual formatting effort and constraints.

### 2.2 Sentiment distribution
| Condition | Positive | Neutral | Negative |
|---|---:|---:|---:|
| Baseline | 20.0% | 70.0% | 10.0% |
| Experimental System | 25.0% | 66.7% | 8.3% |

**Implication:** Tone is mostly neutral-positive in both groups, with slightly more positivity under Experimental; negatives in Baseline highlight manual overhead.

### 2.3 TF-IDF highlights & representative quotes
- Baseline top terms: *formatting*, *manual*, *prompt engineering*, *references* — pointing to effortful authoring.
- Experimental top terms: *idea*, *brainstorm*, *pdf*, *latex*, *helpful*, *generate* — mirroring the strong helpfulness ratings.
- Representative quotes (see `outputs/qualitative/representative_quotes.csv`):
  - Baseline: “It was well defined… but challenging to ensure the question stayed on topic without tooling.”
  - Experimental: “The brainstorming part + developing the idea meaningfully… the Latex generator solved a big problem for me.”

### 2.4 Comment labels
- Each comment is tagged with assigned themes, sentiment, and token counts in `outputs/qualitative/comment_labels.csv` for deeper review or supervised modeling.

### 2.5 Labeled comment examples
| Participant | Condition | Themes | Sentiment | Excerpt |
|---|---|---|---|---|
| 10 | Baseline | design_refinement; helpfulness_ai | negative | “GPT do everything… write a roadmap to design the task for LLMs.” |
| P10 | Experimental System | helpfulness_ai; latex_pdf | positive | “The Brain storm queries are great… would like active feedback; noted a LaTeX compiling error.” |
| P5 | Experimental System | brainstorm_support; helpfulness_ai | neutral | “The brainstorming part + developing the idea meaningfully… details tailor it conveniently.” |

**Implication:** Baseline users describe manual or unclear workflows, whereas Experimental users highlight AI brainstorming/LaTeX benefits alongside actionable feedback for iteration.

### 2.5 Word clouds
- `figs/wordcloud_baseline.svg`: dominated by “formatting”, “time”, “manual”, “prompt”.
- `figs/wordcloud_experimental_system.svg`: dominated by “idea”, “pdf”, “latex”, “helpful”, “generate”.

---

## 3. Practical implications
- **Adopt Experimental workflow:** Expect ~13-minute time savings and large reductions in cognitive load.
- **Ideation support:** Brainstorm + design modules effectively boost creativity and ease-of-idea generation; emphasize them in onboarding.
- **Formatting automation:** PDF/LaTeX generation eliminates a frequent Baseline pain point; ensure the feature is reliable (monitor compile errors).
- **Guidance for future improvements:** Address Baseline concerns (prompt engineering, manual formatting) and monitor Experimental UI/UX issues (occasional bugs, desire for editable brainstorm history).

---

## Notes
- Numeric stats reference `outputs/numeric_summary.csv` and `outputs/numeric_differences.csv`.
- Statistical tests with bootstrap CIs and Holm correction are in `outputs/stats/numeric_tests.csv` and `outputs/stats/categorical_tests.csv`.
- Qualitative artifacts (themes, sentiment, TF-IDF, quotes) live in `outputs/qualitative/`.
- All figures are SVG for easy embedding (`figs/`).

## Condition breakdown (comparative)
| Metric | Baseline Count | Baseline % | Experimental System Count | Experimental System % |
|---|---:|---:|---:|---:|
| Responses | 12 | 50.0% | 12 | 50.0% |

## Task completion (comparative)
| Completion | Baseline Count | Baseline % | Experimental System Count | Experimental System % |
|---|---:|---:|---:|---:|
| Partially | 1 | 8.3% | 0 | 0.0% |
| Yes, Completely | 10 | 83.3% | 11 | 91.7% |
| Yes, Mostly (minor elements incomplete) | 1 | 8.3% | 1 | 8.3% |

## Ratings (1-5) averages
| Question | Average | N |
|---|---:|---:|
| [Frustration Level: How insecure, discouraged, irritated, and stressed versus secure, gratified, content, and relaxed did you feel during the task?] | 2.67 | 24 |
| [Mental Demand: How much mental and perceptual activity was required (e.g., thinking, deciding, calculating)?] | 2.54 | 24 |
| [Temporal Demand (Pace): How much time pressure did you feel due to the rate or pace at which the task needed to be performed?] | 2.38 | 24 |
| Agreement [I am confident that my designed question aligns well with the assignments topic.] | 3.75 | 24 |
| Agreement [I am happy with the level of innovation and creativity of my designed question.] | 3.21 | 24 |
| Agreement [I feel like my designed question doesn't have high pedagogical quality (i.e., it is clear, fair, and effectively tests the intended concepts).] | 2.71 | 24 |
| Agreement [It was easy to come up with an idea for the new question.] | 2.83 | 24 |
| Agreement [It was hard to write and format the new question.] | 2.58 | 24 |

## Per-condition summary (Baseline vs Experimental)
| Condition | N | Avg Minutes | Avg Satisfaction |
|---|---:|---:|---:|
| Baseline | 12 | 24.3 | 3.58 |
| Experimental System | 12 | 10.9 | 4.25 |

### Ratings (1-5) averages by condition (comparative)
| Question | Baseline Avg | Baseline N | Experimental System Avg | Experimental System N |
|---|---:|---:|---:|---:|
| [Frustration Level: How insecure, discouraged, irritated, and stressed versus secure, gratified, content, and relaxed did you feel during the task?] | 3.00 | 12 | 2.33 | 12 |
| [Mental Demand: How much mental and perceptual activity was required (e.g., thinking, deciding, calculating)?] | 3.67 | 12 | 1.42 | 12 |
| [Temporal Demand (Pace): How much time pressure did you feel due to the rate or pace at which the task needed to be performed?] | 3.42 | 12 | 1.33 | 12 |
| Agreement [I am confident that my designed question aligns well with the assignments topic.] | 3.50 | 12 | 4.00 | 12 |
| Agreement [I am happy with the level of innovation and creativity of my designed question.] | 2.83 | 12 | 3.58 | 12 |
| Agreement [I feel like my designed question doesn't have high pedagogical quality (i.e., it is clear, fair, and effectively tests the intended concepts).] | 2.83 | 12 | 2.58 | 12 |
| Agreement [It was easy to come up with an idea for the new question.] | 1.92 | 12 | 3.75 | 12 |
| Agreement [It was hard to write and format the new question.] | 3.17 | 12 | 2.00 | 12 |

### Feature helpfulness distributions (comparative)

#### How helpful was the Brainstorm module in generating ideas for your question?
| Response | Baseline Count | Baseline % | Experimental System Count | Experimental System % |
|---|---:|---:|---:|---:|
| 3 | 0 | 0.0% | 1 | 9.1% |
| 4 | 0 | 0.0% | 4 | 36.4% |
| 5 | 0 | 0.0% | 6 | 54.5% |

#### How helpful was the Design module in refining and structuring your question?
| Response | Baseline Count | Baseline % | Experimental System Count | Experimental System % |
|---|---:|---:|---:|---:|
| 2 | 0 | 0.0% | 1 | 8.3% |
| 3 | 0 | 0.0% | 1 | 8.3% |
| 4 | 0 | 0.0% | 4 | 33.3% |
| 5 | 0 | 0.0% | 6 | 50.0% |

#### How helpful was the Generate module in finalizing your question?
| Response | Baseline Count | Baseline % | Experimental System Count | Experimental System % |
|---|---:|---:|---:|---:|
| 3 | 0 | 0.0% | 3 | 25.0% |
| 4 | 0 | 0.0% | 3 | 25.0% |
| 5 | 0 | 0.0% | 6 | 50.0% |

#### The final document generation feature saved me time (or can save me time in my actual TA process):
| Response | Baseline Count | Baseline % | Experimental System Count | Experimental System % |
|---|---:|---:|---:|---:|
| 5 | 0 | 0.0% | 12 | 100.0% |

## Statistical tests (Baseline vs Experimental)

### Numeric outcomes (permutation tests on mean differences)
| Measure | Baseline Mean | (CI) | Experimental Mean | (CI) | Mean Diff (B - E) | Cohen d | Cliff Delta | Permutation p | Holm p | N B | N E |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Time (minutes) | 24.33 | [17.67, 31.42] | 10.92 | [9.33, 12.92] | 13.42 | 1.45 | 0.68 | 0.0024 | 0.019198080191980802 | 12 | 12 |
| Satisfaction (1-5) | 3.58 | [3.08, 4.08] | 4.25 | [3.92, 4.58] | -0.67 | -0.86 | -0.44 | 0.0844 | 0.44815518448155184 | 12 | 12 |
| Rate your agreement with the following statements [I am confident that my designed question aligns well with the assignments topic.] | 3.50 | [2.67, 4.33] | 4.00 | [3.08, 4.67] | -0.50 | -0.34 | -0.17 | 0.5010 | 1.0 | 12 | 12 |
| Rate your agreement with the following statements [I am happy with the level of innovation and creativity of my designed question.] | 2.83 | [2.25, 3.50] | 3.58 | [2.75, 4.33] | -0.75 | -0.57 | -0.32 | 0.2306 | 0.8807119288071192 | 12 | 12 |
| Rate your agreement with the following statements [I feel like my designed question doesn't have high pedagogical quality (i.e., it is clear, fair, and effectively tests the intended concepts).] | 2.83 | [2.25, 3.33] | 2.58 | [1.92, 3.25] | 0.25 | 0.21 | 0.14 | 0.7226 | 1.0 | 12 | 12 |
| Rate your agreement with the following statements [It was easy to come up with an idea for the new question.] | 1.92 | [1.42, 2.50] | 3.75 | [2.92, 4.50] | -1.83 | -1.45 | -0.67 | 0.0047 | 0.0328967103289671 | 12 | 12 |
| Rate your agreement with the following statements [It was hard to write and format the new question.] | 3.17 | [2.25, 4.08] | 2.00 | [1.42, 2.58] | 1.17 | 0.83 | 0.42 | 0.0747 | 0.44815518448155184 | 12 | 12 |
| Please rate the following aspects of the task you just completed. [Mental Demand: How much mental and perceptual activity was required (e.g., thinking, deciding, calculating)?] | 3.67 | [3.00, 4.25] | 1.42 | [1.17, 1.67] | 2.25 | 2.38 | 0.85 | 0.0002 | 0.001999800019998 | 12 | 12 |
| Please rate the following aspects of the task you just completed. [Temporal Demand (Pace): How much time pressure did you feel due to the rate or pace at which the task needed to be performed?] | 3.42 | [2.58, 4.17] | 1.33 | [1.08, 1.58] | 2.08 | 1.86 | 0.75 | 0.0005 | 0.0044995500449955 | 12 | 12 |
| Please rate the following aspects of the task you just completed. [Frustration Level: How insecure, discouraged, irritated, and stressed versus secure, gratified, content, and relaxed did you feel during the task?] | 3.00 | [2.33, 3.67] | 2.33 | [1.83, 2.92] | 0.67 | 0.58 | 0.35 | 0.2202 | 0.8807119288071192 | 12 | 12 |

### Categorical outcomes (permutation chi-square tests)
| Measure | Chi-square Stat | Permutation p | Holm p | N B | N E |
|---|---:|---:|---:|---:|---:|
| Task Completion | 1.05 | 1.0000 | 1.0 | 12 | 12 |
| How helpful was the Brainstorm module in generating ideas for your question? | 0.00 | 1.0000 | 1.0 | 0 | 11 |
| How helpful was the Design module in refining and structuring your question? | 0.00 | 1.0000 | 1.0 | 0 | 12 |
| How helpful was the Generate module in finalizing your question? | 0.00 | 1.0000 | 1.0 | 0 | 12 |
| The final document generation feature saved me time (or can save me time in my actual TA process): | 0.00 | 1.0000 | 1.0 | 0 | 12 |

Notes: p-values are permutation-based (two-sided), robust for small samples and non-normality.

## Notes
- Likert-style answers were converted by extracting the leading integer (1-5).
- Time fields were normalized by extracting the first number as minutes.
- See the `outputs/` folder for CSV summaries you can load into Excel or Python.
