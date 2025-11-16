#!/usr/bin/env python3
import csv
import os
import re
from collections import Counter, defaultdict
from statistics import mean, median, pstdev
from typing import Dict, List, Tuple, Any, Set
import random
import math

INPUT_CSV = "OCJS Resubmission - ADATA Post-Task Questionnaire (Responses) - Form Responses 1.csv"
OUTPUT_DIR = "outputs"
REPORT_MD = "report.md"
FIGS_DIR = "figs"


def ensure_output_dir(path: str) -> None:
    if not os.path.exists(path):
        os.makedirs(path)


def read_rows(csv_path: str) -> Tuple[List[str], List[Dict[str, str]]]:
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        rows = [row for row in reader]
    return fieldnames, rows

def read_all_questionnaires() -> Tuple[List[str], List[Dict[str, str]]]:
    """
    Load all CSV files that match the OCJS questionnaire naming pattern,
    align them to a union schema, and concatenate the rows.
    """
    files = sorted([
        fn for fn in os.listdir(".")
        if fn.startswith("OCJS Resubmission") and fn.endswith(".csv")
    ])
    if not files:
        return read_rows(INPUT_CSV)
    union_fields: List[str] = []
    seen: Set[str] = set()
    per_file_data: List[Tuple[List[str], List[Dict[str, str]]]] = []
    for fn in files:
        flds, rows = read_rows(fn)
        per_file_data.append((flds, rows))
        for f in flds:
            if f not in seen:
                seen.add(f)
                union_fields.append(f)
    merged_rows: List[Dict[str, str]] = []
    for flds, rows in per_file_data:
        for r in rows:
            merged_row = {field: r.get(field, "") for field in union_fields}
            merged_rows.append(merged_row)
    return union_fields, merged_rows

def parse_leading_int(value: str) -> Any:
    """
    Extract leading integer from values like:
    - '5 (Strongly Agree)'
    - '1 (Very Low)'
    - '3'
    Returns int if found, else None.
    """
    if value is None:
        return None
    s = str(value).strip()
    m = re.match(r"^\s*(\d+)", s)
    if m:
        try:
            return int(m.group(1))
        except ValueError:
            return None
    return None


def parse_minutes(value: str) -> Any:
    """
    Extract integer minutes from a variety of formats:
    - '15'
    - '15min'
    - '10 minutes'
    - '45'
    Returns int if found, else None.
    """
    if value is None:
        return None
    s = str(value).strip().lower()
    # Replace common tokens and keep digits
    m = re.search(r"(\d+)", s)
    if m:
        try:
            return int(m.group(1))
        except ValueError:
            return None
    return None


def write_csv(path: str, header: List[str], rows: List[List[Any]]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for r in rows:
            writer.writerow(r)

def percent(count: int, total: int) -> str:
    return f"{(count / total * 100.0):.1f}%" if total else "0.0%"


def summarize_counter(counter: Counter) -> List[List[Any]]:
    total = sum(counter.values())
    rows = []
    for key, count in counter.most_common():
        pct = (count / total * 100.0) if total else 0.0
        rows.append([key, count, f"{pct:.1f}%"])
    rows.append(["TOTAL", total, "100.0%" if total else "0.0%"])
    return rows

def compute_cohen_d(a: List[float], b: List[float]) -> float:
    if not a or not b:
        return float("nan")
    ma = mean(a)
    mb = mean(b)
    def var(xs: List[float]) -> float:
        m = mean(xs)
        return sum((x - m) ** 2 for x in xs) / (len(xs) - 1) if len(xs) > 1 else 0.0
    va = var(a)
    vb = var(b)
    sp = ((len(a) - 1) * va + (len(b) - 1) * vb) / (len(a) + len(b) - 2) if (len(a) + len(b) - 2) > 0 else 0.0
    sd = (sp ** 0.5) if sp > 0 else 0.0
    return (ma - mb) / sd if sd > 0 else (float("inf") if ma != mb else 0.0)

def compute_cliffs_delta(a: List[float], b: List[float]) -> float:
    if not a or not b:
        return float("nan")
    n_greater = 0
    n_less = 0
    for x in a:
        for y in b:
            if x > y:
                n_greater += 1
            elif x < y:
                n_less += 1
    denom = len(a) * len(b)
    return (n_greater - n_less) / denom if denom else float("nan")

def permutation_pvalue_diff_means(a: List[float], b: List[float], n_perm: int = 10000, seed: int = 42) -> float:
    if not a or not b:
        return float("nan")
    rnd = random.Random(seed)
    observed = abs(mean(a) - mean(b))
    pooled = list(a + b)
    n_a = len(a)
    count = 0
    for _ in range(n_perm):
        rnd.shuffle(pooled)
        a_s = pooled[:n_a]
        b_s = pooled[n_a:]
        stat = abs(mean(a_s) - mean(b_s))
        if stat >= observed - 1e-12:
            count += 1
    return (count + 1) / (n_perm + 1)

def chi_square_stat_from_counts(labels: List[str], a_counts: Counter, b_counts: Counter) -> float:
    cats = labels
    a_total = sum(a_counts.values())
    b_total = sum(b_counts.values())
    grand_total = a_total + b_total
    if grand_total == 0:
        return 0.0
    chi2 = 0.0
    for c in cats:
        o_a = a_counts.get(c, 0)
        o_b = b_counts.get(c, 0)
        col_total = o_a + o_b
        if col_total == 0:
            continue
        e_a = col_total * (a_total / grand_total)
        e_b = col_total * (b_total / grand_total)
        if e_a > 0:
            chi2 += (o_a - e_a) ** 2 / e_a
        if e_b > 0:
            chi2 += (o_b - e_b) ** 2 / e_b
    return chi2

def permutation_pvalue_contingency(a_labels: List[str], b_labels: List[str], categories: List[str], n_perm: int = 10000, seed: int = 42) -> float:
    rnd = random.Random(seed)
    observed_counts_a = Counter(a_labels)
    observed_counts_b = Counter(b_labels)
    observed = chi_square_stat_from_counts(categories, observed_counts_a, observed_counts_b)
    pooled = list(a_labels + b_labels)
    n_a = len(a_labels)
    count = 0
    for _ in range(n_perm):
        rnd.shuffle(pooled)
        a_s = pooled[:n_a]
        b_s = pooled[n_a:]
        cs_a = Counter(a_s)
        cs_b = Counter(b_s)
        stat = chi_square_stat_from_counts(categories, cs_a, cs_b)
        if stat >= observed - 1e-12:
            count += 1
    return (count + 1) / (n_perm + 1)

def bootstrap_mean_ci(values: List[float], iters: int = 5000, alpha: float = 0.05, seed: int = 42) -> Tuple[float, float]:
    if not values:
        return (float("nan"), float("nan"))
    rnd = random.Random(seed)
    n = len(values)
    samples = []
    for _ in range(iters):
        draw = [values[rnd.randrange(0, n)] for _ in range(n)]
        samples.append(mean(draw))
    samples.sort()
    low_idx = max(0, int((alpha / 2) * iters) - 1)
    high_idx = min(iters - 1, int((1 - alpha / 2) * iters) - 1)
    return samples[low_idx], samples[high_idx]

def holm_adjust(p_values: List[float]) -> List[float]:
    if not p_values:
        return []
    indexed = sorted([(p, i) for i, p in enumerate(p_values)])
    m = len(p_values)
    adjusted = [0.0] * m
    running = 0.0
    for rank, (p, idx) in enumerate(indexed, start=1):
        adj = (m - rank + 1) * p
        adj = max(adj, running)
        adjusted[rank - 1] = min(1.0, adj)
        running = adjusted[rank - 1]
    result = [0.0] * m
    for (p, original_idx), adj in zip(indexed, adjusted):
        result[original_idx] = adj
    return result


def main() -> None:
    ensure_output_dir(OUTPUT_DIR)
    ensure_output_dir(FIGS_DIR)
    fields, rows = read_all_questionnaires()

    # Exclude specific participants from all analyses
    participant_col = "Participant ID (P#)"
    exclude_ids = {"p13", "p14", "p15", "p16"}
    if participant_col in fields:
        rows = [
            r for r in rows
            if str(r.get(participant_col, "")).strip().lower() not in exclude_ids
        ]

    # Identify key columns by prefix matching to be robust to exact wording
    col_condition = "Condition"
    col_task_completed = "Did you successfully complete the task of designing a new question?"
    col_time_spent = "How much time do you think you spent on this task? (in minutes)"
    col_satisfaction = "How satisfied are you with the quality of the question you designed?"

    # Agreement/ratings columns (parse leading int 1-5)
    agreement_prefixes = [
        "Rate your agreement with the following statements [I am confident",
        "Rate your agreement with the following statements [I am happy",
        "Rate your agreement with the following statements [I feel like my designed question",
        "Rate your agreement with the following statements [It was easy",
        "Rate your agreement with the following statements [It was hard",
    ]
    tlx_cols = [
        "Please rate the following aspects of the task you just completed. [Mental Demand",
        "Please rate the following aspects of the task you just completed. [Temporal Demand",
        "Please rate the following aspects of the task you just completed. [Frustration Level",
    ]

    helpfulness_cols = [
        "How helpful was the Brainstorm module in generating ideas for your question?",
        "How helpful was the Design module in refining and structuring your question?",
        "How helpful was the Generate module in finalizing your question?",
        "The final document generation feature saved me time (or can save me time in my actual TA process):",
    ]

    # Map actual column names present
    def find_cols_by_prefix(prefix_list: List[str]) -> List[str]:
        found = []
        for p in prefix_list:
            for f in fields:
                if f.startswith(p):
                    found.append(f)
        return found

    agreement_cols = find_cols_by_prefix(agreement_prefixes)
    tlx_actual_cols = find_cols_by_prefix(tlx_cols)
    helpfulness_actual_cols = [c for c in helpfulness_cols if c in fields]

    # Basic counts
    n_responses = len(rows)
    by_condition = Counter([r.get(col_condition, "").strip() for r in rows if r.get(col_condition)])
    by_task_completed = Counter([r.get(col_task_completed, "").strip() for r in rows if r.get(col_task_completed)])

    # Time
    minutes_list: List[int] = []
    for r in rows:
        mins = parse_minutes(r.get(col_time_spent, ""))
        if mins is not None:
            minutes_list.append(mins)
    avg_minutes = mean(minutes_list) if minutes_list else None

    # Satisfaction (1-5)
    satisfaction_list: List[int] = []
    for r in rows:
        v = parse_leading_int(r.get(col_satisfaction, ""))
        if v is not None:
            satisfaction_list.append(v)
    avg_satisfaction = mean(satisfaction_list) if satisfaction_list else None

    # Agreement/ratings
    per_col_numeric: Dict[str, List[int]] = defaultdict(list)
    for c in agreement_cols + tlx_actual_cols:
        for r in rows:
            v = parse_leading_int(r.get(c, ""))
            if v is not None:
                per_col_numeric[c].append(v)

    per_col_avg: List[Tuple[str, float, int]] = []
    for c, vals in per_col_numeric.items():
        if vals:
            per_col_avg.append((c, float(mean(vals)), len(vals)))
    per_col_avg.sort(key=lambda x: x[0])

    # Helpfulness distributions
    helpfulness_dist: Dict[str, Counter] = {}
    for c in helpfulness_actual_cols:
        helpfulness_dist[c] = Counter([r.get(c, "").strip() for r in rows if r.get(c)])

    # Grouped by Condition
    conditions = sorted([k for k in by_condition.keys() if k])
    rows_by_cond: Dict[str, List[Dict[str, str]]] = {c: [] for c in conditions}
    for r in rows:
        cond_val = r.get(col_condition, "").strip()
        if cond_val in rows_by_cond:
            rows_by_cond[cond_val].append(r)

    # Per-condition summary metrics (N, avg minutes, avg satisfaction)
    per_cond_summary: List[List[Any]] = []
    for cond in conditions:
        rs = rows_by_cond[cond]
        n = len(rs)
        mins_vals: List[int] = []
        sat_vals: List[int] = []
        for r in rs:
            mm = parse_minutes(r.get(col_time_spent, ""))
            if mm is not None:
                mins_vals.append(mm)
            sv = parse_leading_int(r.get(col_satisfaction, ""))
            if sv is not None:
                sat_vals.append(sv)
        avg_m = f"{mean(mins_vals):.1f}" if mins_vals else ""
        avg_s = f"{mean(sat_vals):.2f}" if sat_vals else ""
        per_cond_summary.append([cond, n, avg_m, avg_s])
    write_csv(
        os.path.join(OUTPUT_DIR, "condition_summary.csv"),
        ["Condition", "N", "Avg Minutes", "Avg Satisfaction (1-5)"],
        per_cond_summary,
    )

    # Task completion by condition
    task_completion_by_cond_rows: List[List[Any]] = []
    for cond in conditions:
        rs = rows_by_cond[cond]
        total = len([r for r in rs if r.get(col_task_completed)])
        cnt = Counter([r.get(col_task_completed, "").strip() for r in rs if r.get(col_task_completed)])
        for label, count in cnt.most_common():
            task_completion_by_cond_rows.append([cond, label, count, percent(count, total)])
        task_completion_by_cond_rows.append([cond, "TOTAL", total, "100.0%" if total else "0.0%"])
    write_csv(
        os.path.join(OUTPUT_DIR, "task_completion_by_condition.csv"),
        ["Condition", "Completion", "Count", "Percent"],
        task_completion_by_cond_rows,
    )

    # Ratings (agreement + TLX) averages by condition
    rating_by_cond_rows: List[List[Any]] = []
    for cond in conditions:
        rs = rows_by_cond[cond]
        for c in agreement_cols + tlx_actual_cols:
            vals: List[int] = []
            for r in rs:
                v = parse_leading_int(r.get(c, ""))
                if v is not None:
                    vals.append(v)
            if vals:
                rating_by_cond_rows.append([cond, c, f"{mean(vals):.2f}", len(vals)])
    write_csv(
        os.path.join(OUTPUT_DIR, "rating_averages_by_condition.csv"),
        ["Condition", "Question", "Average (1-5)", "N"],
        rating_by_cond_rows,
    )

    # Build comparative (wide) ratings table: one row per question with per-condition avg and N
    # Determine canonical question order from overall per_col_avg
    questions_order = [q for (q, _, _) in per_col_avg]
    # Map (cond, question) -> (avg, n)
    cond_q_to_vals: Dict[Tuple[str, str], Tuple[float, int]] = {}
    for cond in conditions:
        rs = rows_by_cond[cond]
        for q in questions_order:
            vals: List[int] = []
            for r in rs:
                v = parse_leading_int(r.get(q, ""))
                if v is not None:
                    vals.append(v)
            if vals:
                cond_q_to_vals[(cond, q)] = (float(mean(vals)), len(vals))
    # Prepare CSV header dynamically for conditions
    wide_header = ["Question"]
    for cond in conditions:
        wide_header += [f"{cond} Avg", f"{cond} N"]
    wide_rows: List[List[Any]] = []
    for q in questions_order:
        row: List[Any] = [q]
        for cond in conditions:
            avg_n = cond_q_to_vals.get((cond, q))
            if avg_n:
                row += [f"{avg_n[0]:.2f}", avg_n[1]]
            else:
                row += ["", ""]
        wide_rows.append(row)
    write_csv(
        os.path.join(OUTPUT_DIR, "rating_averages_by_condition_comparative.csv"),
        wide_header,
        wide_rows,
    )

    # Helpfulness distributions by condition
    for c in helpfulness_actual_cols:
        rows_out: List[List[Any]] = []
        for cond in conditions:
            rs = rows_by_cond[cond]
            answers = [r.get(c, "").strip() for r in rs if r.get(c)]
            cnt = Counter(answers)
            total = sum(cnt.values())
            for label, count in cnt.most_common():
                rows_out.append([cond, label, count, percent(count, total)])
            rows_out.append([cond, "TOTAL", total, "100.0%" if total else "0.0%"])
        fname = re.sub(r"[^a-z0-9]+", "_", c.lower()).strip("_") + "_by_condition.csv"
        write_csv(
            os.path.join(OUTPUT_DIR, fname),
            ["Condition", "Response", "Count", "Percent"],
            rows_out,
        )

    # -------------------- Numeric descriptive summaries --------------------
    def pretty_label(text: str) -> str:
        return text.replace("Rate your agreement with the following statements ", "Agreement ").replace("Please rate the following aspects of the task you just completed. ", "").strip()

    metric_specs: List[Tuple[str, Any]] = []
    metric_specs.append(("Time (minutes)", lambda r: parse_minutes(r.get(col_time_spent, ""))))
    metric_specs.append(("Satisfaction (1-5)", lambda r: parse_leading_int(r.get(col_satisfaction, ""))))
    for col in agreement_cols + tlx_actual_cols:
        metric_specs.append((pretty_label(col), lambda r, c=col: parse_leading_int(r.get(c, ""))))
    if tlx_actual_cols:
        def tlx_composite(row: Dict[str, str]) -> Any:
            vals = [parse_leading_int(row.get(c, "")) for c in tlx_actual_cols]
            vals = [v for v in vals if v is not None]
            if not vals:
                return None
            return sum(vals) / len(vals)
        metric_specs.append(("NASA TLX Composite", tlx_composite))

    def calc_stats(values: List[float]) -> Tuple[str, str, str, str, str]:
        if not values:
            return ("", "", "", "", "")
        mu = mean(values)
        med = median(values)
        sd = pstdev(values) if len(values) > 1 else 0.0
        return (f"{mu:.2f}", f"{med:.2f}", f"{sd:.2f}", f"{min(values):.2f}", f"{max(values):.2f}")

    numeric_summary_rows: List[List[Any]] = []
    numeric_diff_rows: List[List[Any]] = []

    for label, getter in metric_specs:
        overall_vals: List[float] = []
        per_cond_vals: Dict[str, List[float]] = {c: [] for c in conditions}
        for row in rows:
            val = getter(row)
            if val is None:
                continue
            overall_vals.append(float(val))
            cond = row.get(col_condition, "").strip()
            if cond in per_cond_vals:
                per_cond_vals[cond].append(float(val))
        # overall stats
        mean_s, median_s, sd_s, min_s, max_s = calc_stats(overall_vals)
        numeric_summary_rows.append([label, "Overall", len(overall_vals), mean_s, median_s, sd_s, min_s, max_s])
        for cond in conditions:
            vals = per_cond_vals.get(cond, [])
            mean_s, median_s, sd_s, min_s, max_s = calc_stats(vals)
            numeric_summary_rows.append([label, cond or "(empty)", len(vals), mean_s, median_s, sd_s, min_s, max_s])
        if len(per_cond_vals.get("Baseline", [])) and len(per_cond_vals.get("Experimental System", [])):
            b_vals = per_cond_vals["Baseline"]
            e_vals = per_cond_vals["Experimental System"]
            diff = mean(b_vals) - mean(e_vals)
            cd = compute_cohen_d(b_vals, e_vals)
            cliffs = compute_cliffs_delta(b_vals, e_vals)
            numeric_diff_rows.append([label, f"{mean(b_vals):.2f}", f"{mean(e_vals):.2f}", f"{diff:.2f}", f"{cd:.2f}", f"{cliffs:.2f}", len(b_vals), len(e_vals)])

    write_csv(
        os.path.join(OUTPUT_DIR, "numeric_summary.csv"),
        ["Metric", "Condition", "N", "Mean", "Median", "Std Dev", "Min", "Max"],
        numeric_summary_rows,
    )
    write_csv(
        os.path.join(OUTPUT_DIR, "numeric_differences.csv"),
        ["Metric", "Baseline Mean", "Experimental Mean", "Mean Diff (B - E)", "Cohen d", "Cliff Delta", "N Baseline", "N Experimental"],
        numeric_diff_rows,
    )

    # -------------------- Qualitative analysis --------------------
    QUAL_DIR = os.path.join(OUTPUT_DIR, "qualitative")
    ensure_output_dir(QUAL_DIR)

    COMMENT_HINTS = [
        "What did you like most",
        "What was the most challenging",
        "If you could've improved",
        "Any additional comments",
    ]

    def detect_comment_cols(all_fields: List[str]) -> List[str]:
        cols: List[str] = []
        for f in all_fields:
            if any(f.strip().startswith(h) for h in COMMENT_HINTS):
                cols.append(f)
        if cols:
            return cols
        # fallback: take last few columns that are free-text-ish
        for f in all_fields[-8:]:
            if "Rank the" in f or "Rate the" in f:
                continue
            cols.append(f)
        return cols

    comment_cols = detect_comment_cols(fields)

    def extract_comment_text(row: Dict[str, str]) -> str:
        texts = []
        for col in comment_cols:
            val = row.get(col, "")
            if isinstance(val, str) and val.strip():
                texts.append(val.strip())
        return " ".join(texts).strip()

    STOPWORDS = set("""
    a an and are as at be by for from has have i if in into is it its of on or our so than that the their them they this to was we were what when where which who why will with you your
    """.split())
    WORD_RE = re.compile(r"[A-Za-z]{2,}")

    def tokenize(text: str) -> List[str]:
        return [t for t in (w.lower() for w in WORD_RE.findall(text)) if t not in STOPWORDS]

    THEME_KEYWORDS: Dict[str, List[str]] = {
        "brainstorm_support": ["brainstorm", "idea", "generate", "generation", "inspiration"],
        "design_refinement": ["design", "refine", "refinement", "edit", "enhance", "structure"],
        "latex_pdf": ["latex", "pdf", "document", "compile", "format"],
        "usability_ui": ["ui", "ux", "interface", "bug", "bugs", "screen", "responsive"],
        "innovation_quality": ["innovation", "creative", "creativity", "quality", "pedagogical"],
        "formatting_effort": ["formatting", "writing", "hard", "difficult", "manual"],
        "helpfulness_ai": ["helpful", "hallucination", "rag", "llm", "ai"],
        "time_saving": ["time", "save", "saving", "faster", "speed", "quick"],
    }
    THEME_NAMES = list(THEME_KEYWORDS.keys())

    POSITIVE_WORDS = {"good", "great", "helpful", "love", "easy", "fast", "convenient", "positive", "awesome", "smooth"}
    NEGATIVE_WORDS = {"hard", "difficult", "bad", "slow", "bug", "issue", "problem", "confusing", "annoying", "frustrating"}

    def assign_themes(text: str) -> List[str]:
        low = text.lower()
        tags = [theme for theme, kws in THEME_KEYWORDS.items() if any(kw in low for kw in kws)]
        return tags if tags else ["general"]

    def simple_sentiment(text: str) -> Tuple[str, int, int]:
        tokens = tokenize(text)
        pos = sum(1 for t in tokens if t in POSITIVE_WORDS)
        neg = sum(1 for t in tokens if t in NEGATIVE_WORDS)
        if pos > neg:
            return "positive", pos, neg
        if neg > pos:
            return "negative", pos, neg
        return "neutral", pos, neg

    comment_records: List[List[Any]] = []
    theme_counts_by_cond: Dict[str, Counter] = {c: Counter() for c in conditions}
    sentiment_counts_by_cond: Dict[str, Counter] = {c: Counter() for c in conditions}
    comments_per_cond: Counter = Counter()
    cond_docs: Dict[str, List[str]] = {c: [] for c in conditions}

    for idx, row in enumerate(rows):
        cond = row.get(col_condition, "").strip()
        text = extract_comment_text(row)
        if not text:
            continue
        comments_per_cond[cond] += 1
        tokens = tokenize(text)
        cond_docs.setdefault(cond, []).extend(tokens)
        themes = assign_themes(text)
        for th in themes:
            theme_counts_by_cond.setdefault(cond, Counter())[th] += 1
        sentiment_label, pos_count, neg_count = simple_sentiment(text)
        sentiment_counts_by_cond.setdefault(cond, Counter())[sentiment_label] += 1
        participant = row.get(participant_col, "").strip() if participant_col in row else f"resp_{idx+1}"
        comment_records.append([
            participant,
            cond,
            text,
            ", ".join(sorted(set(themes))),
            sentiment_label,
            pos_count,
            neg_count,
        ])

    write_csv(
        os.path.join(QUAL_DIR, "comment_labels.csv"),
        ["Participant", "Condition", "Comment", "Themes", "Sentiment", "Positive Words", "Negative Words"],
        comment_records,
    )

    theme_rows: List[List[Any]] = []
    for cond in conditions:
        total_comments = max(1, comments_per_cond.get(cond, 0))
        for theme in THEME_NAMES + ["general"]:
            count = theme_counts_by_cond.get(cond, Counter()).get(theme, 0)
            theme_rows.append([
                cond,
                theme,
                count,
                f"{(count / total_comments):.2f}",
            ])
    write_csv(
        os.path.join(QUAL_DIR, "theme_prevalence_by_condition.csv"),
        ["Condition", "Theme", "Mentions", "Mentions per comment"],
        theme_rows,
    )

    sentiment_rows: List[List[Any]] = []
    for cond in conditions:
        total = max(1, sum(sentiment_counts_by_cond.get(cond, Counter()).values()))
        for label in ["positive", "neutral", "negative"]:
            count = sentiment_counts_by_cond.get(cond, Counter()).get(label, 0)
            sentiment_rows.append([
                cond,
                label,
                count,
                f"{(count / total * 100):.1f}%",
            ])
    write_csv(
        os.path.join(QUAL_DIR, "sentiment_by_condition.csv"),
        ["Condition", "Sentiment", "Count", "Percent"],
        sentiment_rows,
    )

    # TF-IDF per condition
    tfidf_rows: List[List[Any]] = []
    all_terms = set()
    term_freq_by_cond: Dict[str, Counter] = {}
    for cond in conditions:
        cnt = Counter(cond_docs.get(cond, []))
        term_freq_by_cond[cond] = cnt
        all_terms.update(cnt.keys())
    df: Dict[str, int] = {term: sum(1 for cond in conditions if term_freq_by_cond[cond].get(term, 0) > 0) for term in all_terms}
    for cond in conditions:
        cnt = term_freq_by_cond[cond]
        total_terms = sum(cnt.values()) or 1
        scored = []
        for term, freq in cnt.items():
            idf = math.log((1 + len(conditions)) / (1 + df.get(term, 0))) + 1
            tfidf = (freq / total_terms) * idf
            scored.append((term, tfidf, freq))
        scored.sort(key=lambda x: x[1], reverse=True)
        for term, score, freq in scored[:50]:
            tfidf_rows.append([cond, term, f"{score:.6f}", freq])
    write_csv(
        os.path.join(QUAL_DIR, "condition_tfidf.csv"),
        ["Condition", "Term", "TF-IDF", "TF"],
        tfidf_rows,
    )

    # Representative quotes per condition
    quotes_rows: List[List[Any]] = []
    quotes_by_cond: Dict[str, List[str]] = defaultdict(list)
    for record in comment_records:
        _, cond, text, _, _, _, _ = record
        quotes_by_cond[cond].append(text)
    for cond in conditions:
        quotes = sorted(quotes_by_cond.get(cond, []), key=lambda s: len(s), reverse=True)[:5]
        for q in quotes:
            quotes_rows.append([cond, q])
    write_csv(
        os.path.join(QUAL_DIR, "representative_quotes.csv"),
        ["Condition", "Quote"],
        quotes_rows,
    )

    # Statistical tests (Baseline vs Experimental)
    stats_dir = os.path.join(OUTPUT_DIR, "stats")
    ensure_output_dir(stats_dir)
    base_rows = rows_by_cond.get("Baseline", [])
    exp_rows = rows_by_cond.get("Experimental System", [])
    numeric_tests_rows: List[List[Any]] = []
    numeric_pvalues: List[float] = []
    # Time
    base_time = [parse_minutes(r.get(col_time_spent, "")) for r in base_rows]
    base_time = [v for v in base_time if v is not None]
    exp_time = [parse_minutes(r.get(col_time_spent, "")) for r in exp_rows]
    exp_time = [v for v in exp_time if v is not None]
    if base_time and exp_time:
        p = permutation_pvalue_diff_means(base_time, exp_time)
        d = compute_cohen_d(base_time, exp_time)
        cd = compute_cliffs_delta(base_time, exp_time)
        ci_b = bootstrap_mean_ci(base_time)
        ci_e = bootstrap_mean_ci(exp_time)
        numeric_tests_rows.append(["Time (minutes)", mean(base_time), mean(exp_time), mean(base_time) - mean(exp_time), d, cd, p, len(base_time), len(exp_time), ci_b[0], ci_b[1], ci_e[0], ci_e[1]])
        numeric_pvalues.append(p)
    # Satisfaction
    base_sat = [parse_leading_int(r.get(col_satisfaction, "")) for r in base_rows]
    base_sat = [v for v in base_sat if v is not None]
    exp_sat = [parse_leading_int(r.get(col_satisfaction, "")) for r in exp_rows]
    exp_sat = [v for v in exp_sat if v is not None]
    if base_sat and exp_sat:
        p = permutation_pvalue_diff_means(base_sat, exp_sat)
        d = compute_cohen_d(base_sat, exp_sat)
        cd = compute_cliffs_delta(base_sat, exp_sat)
        ci_b = bootstrap_mean_ci(base_sat)
        ci_e = bootstrap_mean_ci(exp_sat)
        numeric_tests_rows.append(["Satisfaction (1-5)", mean(base_sat), mean(exp_sat), mean(base_sat) - mean(exp_sat), d, cd, p, len(base_sat), len(exp_sat), ci_b[0], ci_b[1], ci_e[0], ci_e[1]])
        numeric_pvalues.append(p)
    # Per-question numeric (Likert)
    for q in agreement_cols + tlx_actual_cols:
        base_vals = [parse_leading_int(r.get(q, "")) for r in base_rows]
        base_vals = [v for v in base_vals if v is not None]
        exp_vals = [parse_leading_int(r.get(q, "")) for r in exp_rows]
        exp_vals = [v for v in exp_vals if v is not None]
        if base_vals and exp_vals:
            p = permutation_pvalue_diff_means(base_vals, exp_vals)
            d = compute_cohen_d(base_vals, exp_vals)
            cd = compute_cliffs_delta(base_vals, exp_vals)
            ci_b = bootstrap_mean_ci(base_vals)
            ci_e = bootstrap_mean_ci(exp_vals)
            numeric_tests_rows.append([q, mean(base_vals), mean(exp_vals), mean(base_vals) - mean(exp_vals), d, cd, p, len(base_vals), len(exp_vals), ci_b[0], ci_b[1], ci_e[0], ci_e[1]])
            numeric_pvalues.append(p)
    numeric_holm = holm_adjust(numeric_pvalues)
    for idx, row in enumerate(numeric_tests_rows):
        row.append(numeric_holm[idx] if idx < len(numeric_holm) else "")
    write_csv(
        os.path.join(stats_dir, "numeric_tests.csv"),
        ["Measure", "Baseline Mean", "Experimental Mean", "Mean Diff (B - E)", "Cohen d", "Cliff Delta", "Permutation p", "N Baseline", "N Experimental", "Baseline CI Low", "Baseline CI High", "Experimental CI Low", "Experimental CI High", "p_holm"],
        numeric_tests_rows,
    )
    # Categorical tests (task completion + helpfulness)
    cat_tests_rows: List[List[Any]] = []
    cat_pvalues: List[float] = []
    base_comp = [r.get(col_task_completed, "").strip() for r in base_rows if r.get(col_task_completed)]
    exp_comp = [r.get(col_task_completed, "").strip() for r in exp_rows if r.get(col_task_completed)]
    if base_comp and exp_comp:
        cats = sorted(set(base_comp + exp_comp))
        p = permutation_pvalue_contingency(base_comp, exp_comp, cats)
        chi2 = chi_square_stat_from_counts(cats, Counter(base_comp), Counter(exp_comp))
        cat_tests_rows.append(["Task Completion", "; ".join(cats), chi2, p, len(base_comp), len(exp_comp), ""])
        cat_pvalues.append(p)
    for c in helpfulness_actual_cols:
        base_lbls = [r.get(c, "").strip() for r in base_rows if r.get(c)]
        exp_lbls = [r.get(c, "").strip() for r in exp_rows if r.get(c)]
        if base_lbls or exp_lbls:
            cats = sorted(set(base_lbls + exp_lbls))
            if cats:
                p = permutation_pvalue_contingency(base_lbls, exp_lbls, cats)
                chi2 = chi_square_stat_from_counts(cats, Counter(base_lbls), Counter(exp_lbls))
                cat_tests_rows.append([c, "; ".join(cats), chi2, p, len(base_lbls), len(exp_lbls), ""])
                cat_pvalues.append(p)
    cat_holm = holm_adjust(cat_pvalues)
    idx = 0
    for i in range(len(cat_tests_rows)):
        if idx < len(cat_holm):
            cat_tests_rows[i][6] = cat_holm[idx]
            idx += 1
    write_csv(
        os.path.join(stats_dir, "categorical_tests.csv"),
        ["Measure", "Categories", "Chi-square Stat", "Permutation p", "N Baseline", "N Experimental", "p_holm"],
        cat_tests_rows,
    )
    # Save summaries
    write_csv(
        os.path.join(OUTPUT_DIR, "condition_counts.csv"),
        ["Condition", "Count", "Percent"],
        summarize_counter(by_condition),
    )
    write_csv(
        os.path.join(OUTPUT_DIR, "task_completion_counts.csv"),
        ["Task Completion", "Count", "Percent"],
        summarize_counter(by_task_completed),
    )
    # Comparative condition counts CSV
    baseline_n = by_condition.get("Baseline", 0)
    experimental_n = by_condition.get("Experimental System", 0)
    total_n = baseline_n + experimental_n
    write_csv(
        os.path.join(OUTPUT_DIR, "condition_counts_comparative.csv"),
        [
            "Metric",
            "Baseline Count",
            "Baseline Percent",
            "Experimental System Count",
            "Experimental System Percent",
        ],
        [["Responses", baseline_n, percent(baseline_n, total_n), experimental_n, percent(experimental_n, total_n)]],
    )
    # Comparative task completion CSV
    completion_labels = sorted(list(set([k for k in by_task_completed.keys() if k])))
    completion_by_cond: Dict[str, Counter] = {}
    for cond in conditions:
        rs = rows_by_cond[cond]
        completion_by_cond[cond] = Counter([r.get(col_task_completed, "").strip() for r in rs if r.get(col_task_completed)])
    comp_baseline_total = sum(completion_by_cond.get("Baseline", Counter()).values())
    comp_experimental_total = sum(completion_by_cond.get("Experimental System", Counter()).values())
    task_comp_rows: List[List[Any]] = []
    for label in completion_labels:
        bcnt = completion_by_cond.get("Baseline", Counter()).get(label, 0)
        ecnt = completion_by_cond.get("Experimental System", Counter()).get(label, 0)
        task_comp_rows.append([label, bcnt, percent(bcnt, comp_baseline_total), ecnt, percent(ecnt, comp_experimental_total)])
    write_csv(
        os.path.join(OUTPUT_DIR, "task_completion_comparative.csv"),
        ["Completion", "Baseline Count", "Baseline Percent", "Experimental System Count", "Experimental System Percent"],
        task_comp_rows,
    )
    write_csv(
        os.path.join(OUTPUT_DIR, "time_minutes_samples.csv"),
        ["Minutes"],
        [[m] for m in minutes_list],
    )
    write_csv(
        os.path.join(OUTPUT_DIR, "rating_averages.csv"),
        ["Question", "Average (1-5)", "N"],
        [[q, f"{avg:.2f}", n] for (q, avg, n) in per_col_avg],
    )
    for c, cnt in helpfulness_dist.items():
        fname = re.sub(r"[^a-z0-9]+", "_", c.lower()).strip("_") + "_distribution.csv"
        write_csv(
            os.path.join(OUTPUT_DIR, fname),
            ["Response", "Count", "Percent"],
            summarize_counter(cnt),
        )

    # Simple SVG figure utilities (no external libs)
    def save_svg_bar_dual(path: str, title: str, categories: List[str], baseline_vals: List[float], experimental_vals: List[float], max_val: float, xlabel: str) -> None:
        width = 900
        height = 40 + len(categories) * 28 + 60
        bar_h = 10
        gap = 18
        left = 160
        right = 40
        baseline_color = "#4E79A7"
        experimental_color = "#E15759"
        with open(path, "w", encoding="utf-8") as f:
            f.write(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">\n')
            f.write('<style>text{font-family:Arial,sans-serif;font-size:12px}</style>\n')
            f.write(f'<text x="{width/2}" y="20" text-anchor="middle" font-size="16" font-weight="bold">{title}</text>\n')
            axis_y = 30
            axis_w = width - left - right
            for i in range(0, 6):
                x = left + axis_w * i / 5
                val = max_val * i / 5
                f.write(f'<line x1="{x}" y1="{axis_y}" x2="{x}" y2="{height-30}" stroke="#ddd"/>\n')
                f.write(f'<text x="{x}" y="{height-10}" text-anchor="middle">{val:.1f}</text>\n')
            f.write(f'<text x="{width/2}" y="{height-32}" text-anchor="middle">{xlabel}</text>\n')
            f.write(f'<rect x="{left}" y="{axis_y}" width="12" height="12" fill="{baseline_color}"/><text x="{left+18}" y="{axis_y+11}">Baseline</text>\n')
            f.write(f'<rect x="{left+100}" y="{axis_y}" width="12" height="12" fill="{experimental_color}"/><text x="{left+118}" y="{axis_y+11}">Experimental</text>\n')
            y = axis_y + 24
            for idx, cat in enumerate(categories):
                f.write(f'<text x="{left-8}" y="{y+bar_h}" text-anchor="end">{cat}</text>\n')
                bw = 0 if max_val <= 0 else (baseline_vals[idx] / max_val) * axis_w
                ew = 0 if max_val <= 0 else (experimental_vals[idx] / max_val) * axis_w
                f.write(f'<rect x="{left}" y="{y}" width="{bw}" height="{bar_h}" fill="{baseline_color}"/>\n')
                f.write(f'<rect x="{left}" y="{y+bar_h+2}" width="{ew}" height="{bar_h}" fill="{experimental_color}"/>\n')
                f.write(f'<text x="{left+bw+4}" y="{y+bar_h-1}" fill="#333">{baseline_vals[idx]:.2f}</text>\n')
                f.write(f'<text x="{left+ew+4}" y="{y+2*bar_h+1}" fill="#333">{experimental_vals[idx]:.2f}</text>\n')
                y += (bar_h * 2 + gap)
            f.write("</svg>\n")

    # Figures: ratings comparative (averages)
    if 'questions_order' in locals() and questions_order:
        cats = []
        base_vals = []
        exp_vals = []
        for q in questions_order:
            cats.append(q.replace("Rate your agreement with the following statements ", "Agreement ").replace("Please rate the following aspects of the task you just completed. ", ""))
            bv = cond_q_to_vals.get(("Baseline", q), (0.0, 0))[0] if ("Baseline", q) in cond_q_to_vals else 0.0
            ev = cond_q_to_vals.get(("Experimental System", q), (0.0, 0))[0] if ("Experimental System", q) in cond_q_to_vals else 0.0
            base_vals.append(bv)
            exp_vals.append(ev)
        save_svg_bar_dual(os.path.join(FIGS_DIR, "ratings_comparative.svg"), "Ratings (1-5) Averages by Condition", cats, base_vals, exp_vals, 5.0, "Average Rating (1–5)")

    # Figure: time and satisfaction comparative
    base_time_mean = mean(base_time) if 'base_time' in locals() and base_time else 0.0
    exp_time_mean = mean(exp_time) if 'exp_time' in locals() and exp_time else 0.0
    base_sat_mean = mean(base_sat) if 'base_sat' in locals() and base_sat else 0.0
    exp_sat_mean = mean(exp_sat) if 'exp_sat' in locals() and exp_sat else 0.0
    save_svg_bar_dual(
        os.path.join(FIGS_DIR, "time_satisfaction_comparative.svg"),
        "Time and Satisfaction by Condition",
        ["Time (minutes)", "Satisfaction (1–5)"],
        [base_time_mean, base_sat_mean],
        [exp_time_mean, exp_sat_mean],
        max_val=max(max(base_time_mean, exp_time_mean), 5.0),
        xlabel="Value"
    )

    # Figure: task completion comparative (percent)
    if 'task_comp_rows' in locals() and task_comp_rows:
        labels_tc = [r[0] for r in task_comp_rows]
        bperc_tc = [float(r[2].strip('%')) for r in task_comp_rows]
        eperc_tc = [float(r[4].strip('%')) for r in task_comp_rows]
        save_svg_bar_dual(
            os.path.join(FIGS_DIR, "task_completion_comparative.svg"),
            "Task Completion by Condition",
            labels_tc,
            bperc_tc,
            eperc_tc,
            max_val=100.0,
            xlabel="Percent"
        )

    def save_svg_wordcloud(path: str, title: str, word_freq: Counter) -> None:
        items = word_freq.most_common(60)
        if not items:
            return
        width = 900
        height = 500
        cols = 3
        col_w = width // cols
        x_positions = [int(col_w * 0.2), int(col_w * 1.2), int(col_w * 2.2)]
        y_positions = [70, 70, 70]
        max_tf = max(freq for _, freq in items)
        with open(path, "w", encoding="utf-8") as f:
            f.write(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">\n')
            f.write('<style>text{font-family:Arial,sans-serif}</style>\n')
            f.write(f'<text x="{width/2}" y="30" text-anchor="middle" font-size="18" font-weight="bold">{title}</text>\n')
            for idx, (word, freq) in enumerate(items):
                col = idx % cols
                font_size = 10 + int(20 * (freq / max_tf))
                x = x_positions[col]
                y = y_positions[col]
                f.write(f'<text x="{x}" y="{y}" font-size="{font_size}">{word} ({freq})</text>\n')
                y_positions[col] += font_size + 10
            f.write("</svg>\n")

    for cond in conditions:
        counter = Counter(cond_docs.get(cond, []))
        save_svg_wordcloud(
            os.path.join(FIGS_DIR, f"wordcloud_{re.sub(r'[^a-z0-9]+', '_', cond.lower()).strip('_')}.svg"),
            f"Top terms - {cond}",
            counter,
        )

    # Write report.md
    lines: List[str] = []
    lines.append("# ADATA Post-Task Questionnaire: Basic Analysis\n")
    lines.append(f"- Total responses: {n_responses}\n")
    if minutes_list:
        lines.append(f"- Average time spent: {avg_minutes:.1f} minutes (N={len(minutes_list)})\n")
    if satisfaction_list:
        lines.append(f"- Average satisfaction with question: {avg_satisfaction:.2f} / 5 (N={len(satisfaction_list)})\n")

    lines.append("\n## Condition breakdown (comparative)\n")
    lines.append("| Metric | Baseline Count | Baseline % | Experimental System Count | Experimental System % |\n")
    lines.append("|---|---:|---:|---:|---:|\n")
    lines.append(f"| Responses | {baseline_n} | {percent(baseline_n, total_n)} | {experimental_n} | {percent(experimental_n, total_n)} |\n")

    lines.append("\n## Task completion (comparative)\n")
    lines.append("| Completion | Baseline Count | Baseline % | Experimental System Count | Experimental System % |\n")
    lines.append("|---|---:|---:|---:|---:|\n")
    for label, bcnt, bpct, ecnt, epct in task_comp_rows:
        lines.append(f"| {label} | {bcnt} | {bpct} | {ecnt} | {epct} |\n")

    if per_col_avg:
        lines.append("\n## Ratings (1-5) averages\n")
        lines.append("| Question | Average | N |\n")
        lines.append("|---|---:|---:|\n")
        for q, avg, n in per_col_avg:
            short_q = q.replace("Rate your agreement with the following statements ", "Agreement ").replace("Please rate the following aspects of the task you just completed. ", "")
            lines.append(f"| {short_q} | {avg:.2f} | {n} |\n")

    # Per-condition summary section
    if conditions:
        lines.append("\n## Per-condition summary (Baseline vs Experimental)\n")
        lines.append("| Condition | N | Avg Minutes | Avg Satisfaction |\n")
        lines.append("|---|---:|---:|---:|\n")
        for cond, n, avg_m, avg_s in per_cond_summary:
            lines.append(f"| {cond} | {n} | {avg_m} | {avg_s} |\n")

        # (Per-condition task completion subtables removed; see comparative section above)

        # Ratings averages by condition
        # Comparative ratings table (wide)
        lines.append("\n### Ratings (1-5) averages by condition (comparative)\n")
        # header
        header_cells = ["Question"]
        for cond in conditions:
            header_cells += [f"{cond} Avg", f"{cond} N"]
        lines.append("| " + " | ".join(header_cells) + " |\n")
        sep_cells = ["---"] + ["---:" for _ in header_cells[1:]]
        lines.append("|" + "|".join(sep_cells) + "|\n")
        # rows
        for q in questions_order:
            short_q = q.replace("Rate your agreement with the following statements ", "Agreement ").replace("Please rate the following aspects of the task you just completed. ", "")
            row_cells = [short_q]
            for cond in conditions:
                avg_n = cond_q_to_vals.get((cond, q))
                if avg_n:
                    row_cells += [f"{avg_n[0]:.2f}", str(avg_n[1])]
                else:
                    row_cells += ["", ""]
            lines.append("| " + " | ".join(row_cells) + " |\n")

        # Helpfulness distributions (comparative)
        if helpfulness_actual_cols:
            lines.append("\n### Feature helpfulness distributions (comparative)\n")
            for c in helpfulness_actual_cols:
                lines.append(f"\n#### {c}\n")
                base_cnt = Counter([r.get(c, "").strip() for r in rows_by_cond.get("Baseline", []) if r.get(c)])
                exp_cnt = Counter([r.get(c, "").strip() for r in rows_by_cond.get("Experimental System", []) if r.get(c)])
                base_total = sum(base_cnt.values())
                exp_total = sum(exp_cnt.values())
                all_labels = sorted(set(list(base_cnt.keys()) + list(exp_cnt.keys())))
                # CSV for this question
                helpful_comp_csv = re.sub(r"[^a-z0-9]+", "_", c.lower()).strip("_") + "_comparative.csv"
                write_csv(
                    os.path.join(OUTPUT_DIR, helpful_comp_csv),
                    ["Response", "Baseline Count", "Baseline Percent", "Experimental System Count", "Experimental System Percent"],
                    [[label, base_cnt.get(label, 0), percent(base_cnt.get(label, 0), base_total), exp_cnt.get(label, 0), percent(exp_cnt.get(label, 0), exp_total)] for label in all_labels],
                )
                # Report table
                lines.append("| Response | Baseline Count | Baseline % | Experimental System Count | Experimental System % |\n")
                lines.append("|---|---:|---:|---:|---:|\n")
                for label in all_labels:
                    bcnt = base_cnt.get(label, 0)
                    ecnt = exp_cnt.get(label, 0)
                    lines.append(f"| {label or '(empty)'} | {bcnt} | {percent(bcnt, base_total)} | {ecnt} | {percent(ecnt, exp_total)} |\n")

    # Statistical tests summary
    lines.append("\n## Statistical tests (Baseline vs Experimental)\n")
    if 'numeric_tests_rows' in locals() and numeric_tests_rows:
        lines.append("\n### Numeric outcomes (permutation tests on mean differences)\n")
        lines.append("| Measure | Baseline Mean | (CI) | Experimental Mean | (CI) | Mean Diff (B - E) | Cohen d | Cliff Delta | Permutation p | Holm p | N B | N E |\n")
        lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for row in numeric_tests_rows:
            measure, mb, me, diff, d, cd, p, nb, ne, b_lo, b_hi, e_lo, e_hi, p_holm = row
            lines.append(f"| {measure} | {mb:.2f} | [{b_lo:.2f}, {b_hi:.2f}] | {me:.2f} | [{e_lo:.2f}, {e_hi:.2f}] | {diff:.2f} | {d:.2f} | {cd:.2f} | {p:.4f} | {p_holm if p_holm != '' else ''} | {nb} | {ne} |\n")
    if 'cat_tests_rows' in locals() and cat_tests_rows:
        lines.append("\n### Categorical outcomes (permutation chi-square tests)\n")
        lines.append("| Measure | Chi-square Stat | Permutation p | Holm p | N B | N E |\n")
        lines.append("|---|---:|---:|---:|---:|---:|\n")
        for measure, cats, chi2, p, nb, ne, p_holm in cat_tests_rows:
            lines.append(f"| {measure} | {chi2:.2f} | {p:.4f} | {p_holm if p_holm != '' else ''} | {nb} | {ne} |\n")
        lines.append("\nNotes: p-values are permutation-based (two-sided), robust for small samples and non-normality.\n")

    # Omit separate overall helpfulness section; comparative tables included above

    lines.append("\n## Notes\n")
    lines.append("- Likert-style answers were converted by extracting the leading integer (1-5).\n")
    lines.append("- Time fields were normalized by extracting the first number as minutes.\n")
    lines.append("- See the `outputs/` folder for CSV summaries you can load into Excel or Python.\n")

    with open(REPORT_MD, "w", encoding="utf-8") as f:
        f.writelines(lines)

    print(f"Wrote {REPORT_MD} and CSV summaries in {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()


